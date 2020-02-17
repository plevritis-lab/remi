system.file("extdata", "pathwaygenelist.RData", package = "cabernet")
system.file("extdata", "gbiogrid.RData", package = "cabernet")
system.file("extdata", "LR_pairs_filt.txt", package="cabernet")

#' Clean the data
#'
#' This function loads a file as a matrix, 
#' removes genes with low expression, 
#' and creates cell type-specific gene expression matrices.
#'
#' @param dat Data matrix of normalized log2 transformed TPM counts
#' @param cellmarkers list of cell markers as labeled in the columns of the dataset.
#' @param var Number of standard deviations user would like to cut the data by
#' @return A list of gene expression per cell type
#' @export
#' 
cleanData <- function(dat, cellmarkers, var) {
  cellexps <- list()
  dat.cols <- unique(unlist(lapply(colnames(dat), 
                                   function(x) {strsplit(x, "_")[[1]][1]})))
  for(i in 1:length(cellmarkers)) {
    curr.name <- names(cellmarkers)[i]
    celltype.filt <- dat[,grep(cellmarkers[i], colnames(dat))]
    
    # Removing genes with low expression
    genemean <- rowMeans(celltype.filt)
    b <- as.numeric(bin(genemean, var))
    print(length(b))
    names(b) <- rownames(celltype.filt)
    removegenes <- names(b[b == 1])
    
    # Finalized cleaned data matrix
    cleaned.filt <- celltype.filt[-which(rownames(celltype.filt) %in% removegenes),]
    rownames(cleaned.filt) <- paste0(curr.name, "_", rownames(cleaned.filt))
    colnames(cleaned.filt) <- dat.cols
    
    cellexps[[curr.name]] <- cleaned.filt
  }
  return(cellexps)
}


#' Generating cell type-specific ligand receptor pair network
#'
#' Network used for cabernet algorithm. Nodes represent ligand and receptor genes from
#' all cell types. Edges represent a literature-supported ligand receptor pair. 
#'
#' @param lr.table List of all ligand receptor pairs and their pair information
#' @param celltypes List of all cell types in the dataset
#' @param datlist List of all cell type-specific gene expression matrices
#' @param pgenes Will be removed. List of pathway genes downloaded from KEGG
#' @return List including igraph ligand receptor network and data table of 
#' ligand receptor pairs across all cell types
#' @export
#' 
generateLRnet <- function(lr.table, celltypes, datlist, pgenes) {
  # Creating vectors of cell type-specific ligand, receptor, and pathway genes
  all.rs <- c()
  all.ls <- c()
  all.paths <- c()
  for(i in celltypes) {
    all.rs <- c(all.rs, unique(paste0(i, "_", lr.table$Receptor)))
    all.ls <- c(all.ls, unique(paste0(i, "_", lr.table$Ligand)))
    all.paths <- c(all.paths, paste0(i, "_", pgenes))
  }
  
  all.genes <- unlist(lapply(datlist, function(x) {rownames(x)}))
  
  allcombolr <- expand.grid(all.ls, all.rs, stringsAsFactors=F)
  first.ind <- seq(from=1, by=2, length.out = nrow(allcombolr))
  second.ind <- seq(from=2, by=2, length.out = nrow(allcombolr))
  
  lr.network.pairs <- data.table(allcombolr) %>%
    dplyr::rename(L = Var1) %>%
    dplyr::rename(R = Var2) %>%
    mutate(n1cell = unlist(strsplit(L, split="_"))[first.ind]) %>%
    mutate(n2cell = unlist(strsplit(R, split="_"))[first.ind]) %>%
    mutate(ligand = unlist(strsplit(L, split="_"))[second.ind]) %>%
    mutate(receptor = unlist(strsplit(R, split="_"))[second.ind]) %>%
    mutate(pairname = paste0(ligand, "_", receptor)) %>%
    filter(pairname %in% lr.table$Pair.Name) %>%
    mutate(combo1 = paste0(L, "_", R)) %>%
    mutate(combo2 = paste0(R, "_", L)) %>%
    filter(L %in% all.genes) %>%
    filter(R %in% all.genes)
  
  pairwise.cor <- calculateCor(lr.network.pairs, datlist)
  
  lr.network <- graph_from_edgelist(as.matrix(lr.network.pairs[,c(1,2)]), 
                                    directed=F)
  E(lr.network)$weight <- abs(as.numeric(paste0(pairwise.cor$cor)))
  
  return(list(net=lr.network, mat=lr.network.pairs))
}

#' calculateCor - calculating pairwise correlation for all the LR pairs
#'
#' This is an internal function used in generateLRnet. It calculates 
#' Pearson correlation of all ligand receptor pairs.
#'
#' @param  lr.table List of cell type-specific LR pairs
#' @param dat.list List of gene expression for each cell type
#' @return A table of correlations for each LR pair
#' @export
calculateCor <- function(lr.table, dat.list) {

  allcellexp <- dat.list[[names(dat.list)[1]]]
  for(c in names(dat.list)[2:length(dat.list)]) {
    allcellexp <- rbind(allcellexp, dat.list[[c]])
  }
  
  pairwise.lr <- matrix(NA, nrow=nrow(lr.table), ncol=5)
  
  for(i in 1:nrow(lr.table)) {
    ligand <- lr.table$L[i]
    receptor <- lr.table$R[i]
    
    l.exp <- as.numeric(allcellexp[which(rownames(allcellexp) == ligand),])
    r.exp <- as.numeric(allcellexp[which(rownames(allcellexp) ==
                                           receptor),])
    
    pairwise.lr[i,3] <- mean(l.exp)
    pairwise.lr[i,4] <- mean(r.exp)
    cor.res <- cor.test(l.exp, r.exp)
    pairwise.lr[i,5] <- cor.res$estimate
  } 
  pairwise.lr[,1] <- lr.table$L
  pairwise.lr[,2] <- lr.table$R
  pairwise.lr <- data.frame(pairwise.lr)
  colnames(pairwise.lr) <- c("L", "R", "Lexp", "Rexp", "cor")
  pairwise.lr$cor <- as.numeric(paste0(pairwise.lr$cor))
  
  return(pairwise.lr)
}

#' Calculating eigenvector centrality for each receptor given downstream
#' protein protein-interaction signaling network.
#'
#'
#' @param  lr.table List of cell type-specific LR pairs
#' @param dat.list List of gene expression for each cell type
#' @param pgenelist TO BE REMOVED
#' @param numgenes Number of genes user wants to seed clustering algorithm with
#' @return List of receptors with high scoring eigenvector centrality measurements
#' @export
#' 
pickECgenes <- function(lr.table, dat.list, pgenelist, numgenes, cutoff, seed=10) {
  set.seed(seed)
  pgenes <- unlist(pgenelist)
  
  cell.nets <- list()
  for(c in names(dat.list)) {
    receive.exp <- dat.list[[c]]
    receive.net <- makeNetwork(receive.exp, lr.table$L, pgenes, c)
    cell.nets[[c]] <- receive.net
  }

  cell.ec.list <- list()
  all.ec <- c()
  for(c in names(dat.list)) {
    cell.ec <- eigen_centrality(cell.nets[[c]])$vector
    cell.ec.r <- cell.ec[which(names(cell.ec) %in% lr.table$R)]
    cell.ec.list[[c]] <- cell.ec.r
    all.ec <- c(all.ec, cell.ec.r)
  }
  
  sig.top.r.path <- c()
  for(c in names(dat.list)) {
    plength <- c()
    for(p1 in 1:length(pathway.genelist$genesets)) {
      cell.path.all <- paste0(c, "_", pathway.genelist$genesets[[p1]])
      singler.path <- unique(intersect(cell.path.all, lr.table$R))
      top.genes <- head(sort(cell.ec.list[[c]][singler.path], decreasing=T), 
                        numgenes)
      sig.top.r.path <- c(sig.top.r.path, top.genes)
    }
  }
  
  allcellexp <- dat.list[[names(dat.list)[1]]]
  for(c in names(dat.list)[2:length(dat.list)]) {
    allcellexp <- rbind(allcellexp, dat.list[[c]])
  }

  sig.path.exp <- allcellexp[which(rownames(allcellexp) %in% 
                                     unique(names(sig.top.r.path))),]
  sig.path.high <- sig.path.exp[which(rowMeans(sig.path.exp) > cutoff),]
  
  return(list(eclist=sig.top.r.path[rownames(sig.path.high)], 
              net=cell.nets))
  #return(rownames(sig.path.high))
  #return(list(eclist=sig.path.high, net=cell.nets))
}

#' cabernetCommunities
#' 
#' Community detection ft. label propagation and louvain clustering on
#' the LR network
#'
#' @param  net List of LR network and LR table
#' @param dat.list List of gene expression for each cell type
#' @param labelednodes List of high scoring eigenvector centrality receptors to seed
#' the label propagation clustering
#' @return Communities in LR network 
#' @export
#' 
cabernetCommunities <- function(net, dat.list, lnodes, seed) {
  set.seed(seed)
  # Identifying what components have a labeled receptor
  clu <- igraph::membership(components(net))
  
  labeled.comm.all <- table(clu[names(lnodes)])
  
  #labeled.comm.names <- names(labeled.comm.all[labeled.comm.all>1])
  labeled.comm.names <- names(labeled.comm.all)
  
  # Final list of communities
  lr.communities <- c()
  
  # Clustering using label propagation to seed out important receptors
  labelprop.comms <- clusterLabelProp(net, clu, labeled.comm.names, lnodes)
  
  # Calculate degree of each network and check how many are 
  num.oversized <- calculateOversizedComms(net, labelprop.comms, dat.list[[1]])
  
  # Loop through all the communities until their degrees match the sample size
  old.comms <- labelprop.comms
  old.oversized <- num.oversized
  if(length(num.oversized) == 0) {
    return(labelprop.comms)
  }
  while(length(num.oversized) > 0) {
    louvain.comms <- clusterLouvain(net, num.oversized, old.comms)
    num.oversized <- calculateOversizedComms(net, louvain.comms, dat.list[[1]])
    
    if(length(setdiff(num.oversized, old.oversized)) == 0) {
      next 
    }
    old.comms <- louvain.comms
    old.oversized <- num.oversized
  }
  return(louvain.comms)
}

#' Calculate graphical lasso on each community to sparsify the network
#'
#' @param  netlist List of LR network and table
#' @param communities LR network communities
#' @param datlist List of cell type-specific gene expression matrices
#' @param seednum Number of random seed to ensure reproducible results
#' @return List of edges remaining in LR network
#' @export
#' 
cabernetGlasso <- function(netlist, communities, dat.list, seednum) {
  set.seed(seednum)
  predicted.edges <- list()
  
  allcellexp <- dat.list[[names(dat.list)[1]]]
  for(c in names(dat.list)[2:length(dat.list)]) {
    allcellexp <- rbind(allcellexp, dat.list[[c]])
  }
  
  gedges_w <- matrix("", nrow=0, ncol=8)
  colnames(gedges_w) <- c("node1", "node2", "weight", "cor", 
                          "commnum", "commsize", "deg", "numedges")
  
  size.communities <- table(communities)
  size.communities <- names(size.communities[size.communities > 1])
  
  # Graphical Lasso within edges
  for(i in unique(size.communities)) {
    commgenes <- names(communities[communities == i])
    lr.mat <- t(allcellexp[which(rownames(allcellexp) %in% commgenes),])
    lr.comm.net <- induced_subgraph(netlist$net, colnames(lr.mat))
    
    W <- calculateCommGlasso(lr.mat, netlist)
    lr.W <- data.table(W) %>% 
      mutate(commnum = i) %>%
      mutate(commsize = length(commgenes)) %>%
      dplyr::select(node1, node2, weight, cor, commnum, commsize, deg, numedges)
     gedges_w <- rbind(gedges_w, lr.W)
  }
  
  # Graphical Lasso between edges
  for(i in 1:max(communities)) {
    for(j in seq(i+1, max(communities))) {
      if(i != j) {
        comm1genes <- names(communities)[which(communities == i)]
        comm2genes <- names(communities)[which(communities == j)]
        commgenes <- unique(c(comm1genes, comm2genes))
        all_edges <- E(netlist$net)[inc(commgenes)]
        all_edges_m <- data.table(get.edges(netlist$net, all_edges)) %>%
          mutate(name1 = V(netlist$net)$name[V1]) %>%
          mutate(name2 = V(netlist$net)$name[V2]) %>%
          mutate(comm1 = communities[name1]) %>%
          mutate(comm2 = communities[name2]) %>%
          mutate(sum = comm1+comm2) %>%
          filter(sum == i+j) 
        
        betweengenes <- unique(c(all_edges_m$name1, all_edges_m$name2))
        
        if(length(betweengenes) > 1) {
          
          lr.mat <- t(allcellexp[which(rownames(allcellexp) %in% betweengenes),])
          W <- calculateCommGlasso(lr.mat, netlist)
          betweenlr.W <- data.table(W) %>% 
            mutate(commnum = paste0(i, "_", j)) %>%
            mutate(commsize = length(betweengenes)) %>%
            dplyr::select(node1, node2, weight, cor, commnum, 
                          commsize, deg, numedges) 
          
          gedges_w <- rbind(gedges_w, betweenlr.W)
        }
      }
    }
  }
  
  return(gedges_w)
}


#' Cabernet algorithm
#'
#' Identifies communities in LR network and then performs graphical lasso
#' on the communites in the network to predict significant ligand
#' receptor interactions in the dataset of interest
#'
#' @param  netlist List of LR network and table
#' @param dat.list List of gene expression for each cell type
#' @param labelednodes List of high eigenvector centrality receptors
#' @param seed Random seed number to ensure reproducible results
#' @return List of cabernet predictions
#' @export
#' 
cabernet <- function(netlist, dat.list, labelednodes, seed) {
  # Cluster 

  commdetect.output <- cabernetCommunities(netlist$net, dat.list, labelednodes, seed)
  # Graphical Lasso

  glasso.output <- cabernetGlasso(netlist, commdetect.output, dat.list, seed)
  
  # Return final results
  strsplit.ind <- seq(from=1, by=2, length.out = nrow(glasso.output))
  strsplit.ind2 <- seq(from=2, by=2, length.out = nrow(glasso.output))
  net.edges <- data.table(glasso.output) %>% 
    mutate(n1cell = unlist(strsplit(node1, split="_"))[strsplit.ind]) %>%
    mutate(n2cell = unlist(strsplit(node2, split="_"))[strsplit.ind]) %>%
    mutate(ligand = unlist(strsplit(node1, split="_"))[strsplit.ind2]) %>%
    mutate(receptor = unlist(strsplit(node2, split="_"))[strsplit.ind2]) %>%
    dplyr::select(node1, node2, ligand, receptor, n1cell, n2cell,
                  weight, cor, commsize, commnum, deg, numedges)
  
  # need to switch the cell type
  for(i in 1:nrow(net.edges)) {
    if(net.edges[i,"node1"] %in% lr.net$mat$R) {
      temp.name <- net.edges[i,"node2"]
      temp.cell <- net.edges[i,"n2cell"]
      temp.solo <- net.edges[i,"receptor"]
      net.edges[i,"node2"] <- net.edges[i,"node1"]
      net.edges[i,"node1"] <- temp.name
      net.edges[i,"n2cell"] <- net.edges[i,"n1cell"]
      net.edges[i,"n1cell"] <- temp.cell
      net.edges[i,"receptor"] <- net.edges[i,"ligand"]
      net.edges[i,"ligand"] <- temp.solo
    }
  }
  
  predicted.edges <- net.edges %>%
    mutate(pairname = paste0(node1, "_", node2)) %>%
    filter(pairname %in% netlist$mat$combo1) %>%
    group_by(node1, node2) %>%
    filter(deg == max(deg)) %>%
    filter(numedges == max(as.numeric(numedges))) %>%
    filter(weight == min(as.numeric(weight))) %>%
    ungroup() %>%
    group_by_at(vars(-commnum)) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    filter(as.numeric(weight) > 0) %>%
    mutate(score=as.numeric(weight)*as.numeric(numedges)) %>%
    mutate(softmax=DMwR::SoftMax(as.numeric(score))) %>%
    mutate(scaledweight=scales::rescale(as.numeric(score))) %>%
    mutate(cc = paste0(n1cell, "_", n2cell)) 

  return(list(edges=predicted.edges, mat=commdetect.output))
}
