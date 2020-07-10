system.file("extdata", "pathwaygenelist.RData", package = "cabernet")
system.file("extdata", "gbiogrid.RData", package = "cabernet")
system.file("extdata", "LR_pairs_filt.txt", package="cabernet")

lr.full.pairs <- read_tsv("/Users/ayu1/Documents/cabernet/Data/PairsLigRec.txt")
curr.lr.filt <- lr.full.pairs %>%
  filter(Pair.Evidence == "literature supported") %>%
  dplyr::select(Pair.Name, Ligand.ApprovedSymbol, Receptor.ApprovedSymbol) %>%
  rename(Ligand = Ligand.ApprovedSymbol) %>%
  rename(Receptor = Receptor.ApprovedSymbol) %>%
  add_row(Pair.Name="CD274_PDCD1", Ligand="CD274", Receptor="PDCD1") %>%
  add_row(Pair.Name="CD80_CTLA4", Ligand="CD80", Receptor="CTLA4") %>%
  add_row(Pair.Name="CD80_CD28", Ligand="CD80", Receptor="CD28") %>%
  add_row(Pair.Name="CD86_CD28", Ligand="CD86", Receptor="CD28") %>%
  add_row(Pair.Name="GREM1_KDR", Ligand="GREM1", Receptor="KDR") %>%
  add_row(Pair.Name="PDCD1LG2_PDCD1", Ligand="PDCD1LG2", Receptor="PDCD1") %>%
  add_row(Pair.Name="NECTIN2_CD226", Ligand="NECTIN2", Receptor="CD226") %>%
  add_row(Pair.Name="NECTIN2-TIGHT", Ligand="NECTIN2", Receptor="TIGHT") %>%
  add_row(Pair.Name="PVR-TIGHT", Ligand="PVR", Receptor="TIGHT") %>%
  add_row(Pair.Name="SIGLEC1-SPN", Ligand="SIGLEC1", Receptor="SPN")

load('/Users/ayu1/Documents/cabernet/Data/gbiogrid.RData')
load('/Users/ayu1/Documents/cabernet/Data/pathwaygenelist.RData')
pathway.genes <- unique(unlist(pathway.genelist$genesets))

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
setupData <- function(dat, cellmarkers, filter=T, var = 3) {
  filtered.cellexps <- list()
  notfiltered.cellexps <- list()

  colsvec <- c()
  for(i in 1:length(cellmarkers)) {
    curr.name <- names(cellmarkers)[i]
    celltype.filt <- dat[,grep(cellmarkers[i], colnames(dat))]

    dat.cols <- unique(unlist(lapply(colnames(celltype.filt), 
                                     function(x) {strsplit(x, "_")[[1]][1]})))
    
    rownames(celltype.filt) <- paste0(curr.name, "_", rownames(celltype.filt))
    colnames(celltype.filt) <- dat.cols
    
    colsvec <- c(colsvec, dat.cols)
    
    # Removing genes with low expression
    genemean <- rowMeans(celltype.filt)
    
    b <- as.numeric(OneR::bin(genemean, var))
    names(b) <- rownames(celltype.filt)
    removegenes <- names(b[b == 1])
    
    # Finalized cleaned data matrix
    cleaned.filt <- celltype.filt[-which(rownames(celltype.filt) %in% removegenes),]
    
    filtered.cellexps[[curr.name]] <- cleaned.filt
    notfiltered.cellexps[[curr.name]] <- celltype.filt
  }

  allpresent.cols <- names(which(table(colsvec) == length(cellmarkers))) 

  #Filtering non-filtered data
  for(c in names(notfiltered.cellexps)) {
    allpresent.filt.data <- notfiltered.cellexps[[c]][,allpresent.cols]
    duplicate.cols <- which(apply(allpresent.filt.data, 1,
                                  function(x) length(unique(x))==1) == TRUE)
    if(length(duplicate.cols) > 0) {
      notfiltered.cellexps[[c]] <- allpresent.filt.data[-duplicate.cols,]
    } else {
      notfiltered.cellexps[[c]] <- allpresent.filt.data
    }
  }
  
  for(c in names(filtered.cellexps)) {
    allpresent.filt.data <- filtered.cellexps[[c]][,allpresent.cols]
    duplicate.cols <- which(apply(allpresent.filt.data, 1, 
                                  function(x) length(unique(x))==1) == TRUE)
    if(length(duplicate.cols) > 0) {
      filtered.cellexps[[c]] <- allpresent.filt.data[-duplicate.cols,]
    } else {
      filtered.cellexps[[c]] <- allpresent.filt.data
    }
  }
  
  return(list(filtered=filtered.cellexps, unfiltered=notfiltered.cellexps))
}

#' Generating cell type-specific ligand receptor pair network
#'
#' Network used for cabernet algorithm. Nodes represent ligand and receptor genes from
#' all cell types. Edges represent a literature-supported ligand receptor pair. 
#'
#' @param lr.table List of all ligand receptor pairs and their pair information
#' @param celltypes List of all cell types in the dataset
#' @param datlist List of all cell type-specific gene expression matrices
#' @return List including igraph ligand receptor network and data table of 
#' ligand receptor pairs across all cell types
#' @export
#' 
expandLRpairs <- function(lr.table, datlist, celltypes) {
  # Creating vectors of cell type-specific ligand, receptor, and pathway genes
  all.rs <- c()
  all.ls <- c()
  for(i in celltypes) {
    all.rs <- c(all.rs, unique(paste0(i, "_", lr.table$Receptor)))
    all.ls <- c(all.ls, unique(paste0(i, "_", lr.table$Ligand)))
  }
  
  all.genes <- as.vector(unlist(lapply(datlist, function(x) {rownames(x)})))
  
  allcellexp <- datlist[[names(datlist)[1]]]
  for(c in names(datlist)[2:length(datlist)]) {
    allcellexp <- rbind(allcellexp, datlist[[c]])
  }
  
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
  
  if(nrow(lr.network.pairs) == 0) {
    print("No ligand and receptor pairs in dataset. Try relaxing the gene expression cutoff filter.")
  }
  
  return(lr.network.pairs)
}

#' Generating cell type-specific ligand receptor pair network
#'
#' Network used for cabernet algorithm. Nodes represent ligand and receptor genes from
#' all cell types. Edges represent a literature-supported ligand receptor pair. 
#'
#' @param lr.table List of all ligand receptor pairs and their pair information
#' @param celltypes List of all cell types in the dataset
#' @param datlist List of all cell type-specific gene expression matrices
#' @return List including igraph ligand receptor network and data table of 
#' ligand receptor pairs across all cell types
#' @export
#' 
generateLRnet <- function(lr.table, celltypes, datlist, cor=T) {
  
  print("step1")
  lr.network.pairs <- expandLRpairs(lr.table, datlist, celltypes)
  
  print("step2")
  pairwise.cor <- calculateCor(lr.network.pairs, datlist)
  print("step3")
  lr.network <- graph_from_edgelist(as.matrix(lr.network.pairs[,c(1,2)]), 
                                    directed=F)
  E(lr.network)$weight <- abs(as.numeric(paste0(pairwise.cor$cor)))
  #cor.edges <- as.numeric(paste0(pairwise.cor$cor))
  #E(lr.network)$weight <- cor.edges*cor.edges
  #E(lr.network)$weight <- as.numeric(paste0(pairwise.cor$cor))
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
  
  lr.genes <- unique(c(lr.table$L, lr.table$R))
  lr.allcellexp <- allcellexp[lr.genes,]
  
  pairwise.lr <- matrix(NA, nrow=nrow(lr.table), ncol=3)
  
  lr.cor <- cor(t(lr.allcellexp))
  
  for(i in 1:nrow(lr.table)) {
    ligand <- lr.table$L[i]
    receptor <- lr.table$R[i]
    cor.res <- lr.cor[ligand, receptor]
    pairwise.lr[i,3] <- cor.res
  } 
  pairwise.lr[,1] <- lr.table$L
  pairwise.lr[,2] <- lr.table$R
  pairwise.lr <- data.frame(pairwise.lr)
  colnames(pairwise.lr) <- c("L", "R", "cor")
  pairwise.lr$cor <- as.numeric(paste0(pairwise.lr$cor))
  
  return(pairwise.lr)
}

#' Calculating eigenvector centrality for each receptor given downstream
#' protein protein-interaction signaling network.
#'
#'
#' @param  lrnet List of cell type-specific LR pairs
#' @param dat.list List of gene expression for each cell type
#' @param pgenelist TO BE REMOVED
#' @param numgenes Number of genes user wants to seed clustering algorithm with
#' @return List of receptors with high scoring eigenvector centrality measurements
#' @export
#' 
pickECgenes <- function(lrnet, dat.list, pgenelist, numgenes, cutoff, ppi, seed=10) {
  
  lr.table <- lrnet$mat
  
  set.seed(seed)
  pgenes <- unlist(pgenelist)
  
  cell.nets <- list()
  for(c in names(dat.list)) {
    receive.exp <- dat.list[[c]]
    receive.net <- makePPINetwork(receive.exp, unique(lr.table$L), pgenes, c, ppi)
    if(!is.na(receive.net[1])) {
      cell.nets[[c]] <- receive.net
    }   
  }

  cell.ec.list <- list()
  cell.ec.all.list <- list()
  all.ec <- c()
  for(c in names(cell.nets)) {
    cell.ec <- eigen_centrality(cell.nets[[c]])$vector
    cell.ec.r <- cell.ec[which(names(cell.ec) %in% lr.table$R)]
    cell.ec.list[[c]] <- cell.ec.r
    cell.ec.all.list[[c]] <- cell.ec
    all.ec <- c(all.ec, cell.ec.r)
  }
  
  sig.top.r.path <- c()
  for(c in names(cell.nets)) {
    plength <- c()
    for(p1 in 1:length(pgenelist)) {
      cell.path.all <- paste0(c, "_", pgenelist[[p1]])
      singler.path <- unique(intersect(cell.path.all, lr.table$R))
      
      top.genes <- head(sort(cell.ec.list[[c]][singler.path], decreasing=T), 
                        numgenes)
      top.genes <- top.genes[top.genes > 0]
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
              net=cell.nets,
              allec = sig.top.r.path))
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
  
  labeled.comm.names <- names(labeled.comm.all[labeled.comm.all>1])
  labeled.comm.names <- names(labeled.comm.all)
  
  # Final list of communities
  lr.communities <- c()
  
  # Clustering using label propagation to seed out important receptors
  labelprop.comms <- clusterLabelProp(net, clu, labeled.comm.names, lnodes)
  
  if(length(unique(labelprop.comms)) == 1) {
    labelprop.comms <- membership(cluster_louvain(net))
  }
  
  # Calculate degree of each network and check how many are 
  num.oversized <- calculateOversizedComms(net, labelprop.comms, dat.list[[1]])
  
  # Iterate through all the communities until their degrees match the sample size
  old.comms <- labelprop.comms
  old.oversized <- num.oversized
  
  if(length(num.oversized) == 0) {
    size.communities <- table(labelprop.comms)
    commnames <- unique(names(size.communities[size.communities > 1]))
    return(list(names=commnames, membership=labelprop.comms))
  }
  
  while(length(num.oversized) > 0) {
    louvain.comms <- clusterLouvain(net, num.oversized, old.comms)
    num.oversized <- calculateOversizedComms(net, louvain.comms, dat.list[[1]])
    
    if(length(setdiff(num.oversized, old.oversized)) == 0) {
      break 
    }
    old.comms <- louvain.comms
    old.oversized <- num.oversized
  }
  
  size.communities <- table(louvain.comms)
  commnames <- unique(names(size.communities[size.communities > 1]))
  
  print(paste0(length(commnames), " communities identified"))

  return(list(names=commnames, membership=louvain.comms))
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
cabernetGlasso <- function(netlist, communities, dat.list, seednum, scale=F) {
  set.seed(seednum)
  predicted.edges <- list()
  
  allcellexp <- dat.list[[names(dat.list)[1]]]
  for(c in names(dat.list)[2:length(dat.list)]) {
    allcellexp <- rbind(allcellexp, dat.list[[c]])
  }
  
  gedges_w <- matrix("", nrow=0, ncol=8)
  colnames(gedges_w) <- c("node1", "node2", "weight", "cor", 
                          "commnum", "commsize", "deg", "numgenes")
  
  cat("Removing conditionally independent edges within communities\n")
  # Graphical Lasso within edges
  pbar <- dplyr::progress_estimated(length(communities$names), min_time = 0)
  
  for(i in 1:length(communities$names)) {

    communitynum <- communities$names[i]
    commgenes <- names(communities$membership)[communities$membership == communitynum]

    lr.mat <- t(allcellexp[which(rownames(allcellexp) %in% commgenes),])
    lr.comm.net <- induced_subgraph(netlist$net, colnames(lr.mat))
    
    W <- calculateCommGlasso(cor(lr.mat), lr.mat, netlist, lambda.max=0.9, scale=scale)
    
    lr.W <- data.table(W$W) %>% 
      mutate(commnum = communitynum) %>%
      mutate(commsize = length(commgenes)) %>%
      dplyr::select(node1, node2, weight, cor, commnum, commsize, deg, numgenes)
    gedges_w <- rbind(gedges_w, lr.W)
    
    pbar$pause(0.01)$tick()$print()
  }
  
  if(length(communities$names) <= 1) {
    return(gedges_w)
  }
  
  cat("\nRemoving conditionally independent edges between communities\n")
  pbar <- dplyr::progress_estimated(length(communities$names), min_time = 0)
  
  noadjcomms <- c()
  
  for(i in 1:length(communities$names)) {

    comm1.num <- communities$names[i]
    comm1genes <- names(communities$membership)[which(communities$membership == comm1.num)]
    
    comm1.adjedges <- findAdjacentCommEdges(netlist, comm1genes, communities)
    
    adjcomms <- setdiff(unique(comm1.adjedges$comm1, comm1.adjedges$comm2), comm1.num)
    
    if(length(adjcomms) == 0) {
      noadjcomms <- c(noadjcomms, comm1.num)
    } else {
      for(j in adjcomms) {
        comm2.num <- j
        
        combo1 <- paste0(gedges_w$node1, "_", gedges_w$node2)
        combo2 <- paste0(gedges_w$node2, "_", gedges_w$node1)
        combo <- c(combo1, combo2)
        
        comm1genes <- names(communities$membership)[which(communities$membership == comm1.num)]
        comm2genes <- names(communities$membership)[which(communities$membership == comm2.num)]
        commgenes <- unique(c(comm1genes, comm2genes))
        
        betweencommedges <- findAdjacentCommEdges(netlist, commgenes, communities) %>%
          filter((comm1 == comm1.num & comm2 == comm2.num) | 
                   (comm1 == comm2.num & comm2 == comm1.num)) %>% 
                   filter(comm1 != comm2)
        
        betweengenes <- unique(c(betweencommedges$name1, betweencommedges$name2))
        
        if(length(betweengenes) > 2) {
  
          lr.mat <- t(allcellexp[which(rownames(allcellexp) %in% betweengenes),])
          W <- calculateCommGlasso(cor(lr.mat), lr.mat, netlist, scale=scale, additional=combo)$W
          betweenlr.W <- data.table(W) %>% 
            mutate(commnum = paste0(comm1.num, "_", comm2.num)) %>%
            mutate(commsize = length(betweengenes)) %>%
            dplyr::select(node1, node2, weight, cor, commnum, 
                          commsize, deg, numgenes) 
          
          gedges_w <- rbind(gedges_w, betweenlr.W)
            
        }
      }
      pbar$pause(0.01)$tick()$print()
    }
  }
  
  return(gedges_w)
}

#' Cleaning Output
#'
#' Identifies communities in LR network and then performs graphical lasso
#' on the communites in the network to predict significant ligand
#' receptor interactions in the dataset of interest
#'
#' @param  input Glasso output
#' @param netlist Ligand receptor network
#' @return Cleaned up list of LR pairs with metadata displayed in columns
#' @export
#' 
cleaningOutput <- function(input, netlist) {
  # Swapping LR pairs if they are in the order of RL
  switched.W <- input
  switch.inds <- which(input$node1 %in% netlist$mat$R)
  switched.W$node1[switch.inds] <- input$node2[switch.inds]
  switched.W$node2[switch.inds] <- input$node1[switch.inds]
  
  # Breaking apart columns for more detailed sorting
  strsplit.ind <- seq(from=1, by=2, length.out = nrow(input))
  strsplit.ind2 <- seq(from=2, by=2, length.out = nrow(input))
  
  net.edges <- data.table(switched.W) %>% 
    mutate(n1cell = unlist(strsplit(node1, split="_"))[strsplit.ind]) %>%
    mutate(n2cell = unlist(strsplit(node2, split="_"))[strsplit.ind]) %>%
    mutate(ligand = unlist(strsplit(node1, split="_"))[strsplit.ind2]) %>%
    mutate(receptor = unlist(strsplit(node2, split="_"))[strsplit.ind2]) %>%
    dplyr::select(node1, node2, ligand, receptor, n1cell, n2cell,
                  weight, cor, commsize, commnum, deg, numgenes) %>%
    mutate(pairname = paste0(node1, "_", node2)) %>%
    filter(pairname %in% netlist$mat$combo1) %>%
    group_by(node1, node2) %>%
    #filter(deg == max(deg)) %>%
    #filter(numedges == max(as.numeric(numedges))) %>%
    #filter(weight == min(as.numeric(weight))) %>%
    #ungroup() %>%
    unique() %>%
    #group_by_at(vars(-commnum)) %>%
    group_by(pairname) %>%
    mutate(n = n()) %>%
    #filter(n() == 1) %>%
    ungroup() 
  
  return(net.edges)
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
cabernet <- function(cellmarkers, dat.list, seed, numgenes=2, cutoff=1) {
  
  # Creating LR network
  cat("Building LR network\n")
  netlist <- generateLRnet(curr.lr.filt, names(cellmarkers), 
                          dat.list$filtered, pathway.genes)
  
  #Identify influential receptors
  cat("Identifying influential receptors\n")
  labelednodes <- pickECgenes(netlist, 
                              dat.list$unfiltered, 
                              pathway.genelist$genesets, 
                              numgenes = numgenes,
                              cutoff = cutoff,
                              ppi = g.biogrid,
                              seed=seed)

  # Cluster 
  cat("Detecting communities")
  commdetect.output <- cabernetCommunities(netlist$net, dat.list$filtered, labelednodes$eclist, seed)

  # Graphical Lasso
  cat("\n Estimating activated LR pairs")
  glasso.output <- cabernetGlasso(netlist, commdetect.output, dat.list$filtered, seed)

  predicted.edges <- cleaningOutput(glasso.output, netlist)
  
  numLRs <- length(unique(predicted.edges$node1, predicted.edges$node2))
  
  if(nrow(predicted.edges) == 0) {
    cat("\nNo interactome found\n")
    return(predicted.edges)
  }
  
  #scaled.edges <- predicted.edges %>%
  #  filter(as.numeric(weight) > 0) %>%
  #  mutate(commsizeprob = commsize/numLRs) %>%
  #  mutate(score = as.numeric(weight)*commsizeprob) %>%
  #  mutate(scaledweight=scales::rescale(as.numeric(score))) %>%
  #  mutate(cc = paste0(n1cell, "_", n2cell)) 

  return(list(lrnet=netlist, ir=labelednodes, 
              #predinteractome=scaled.edges, 
              wholeinteractome=predicted.edges,
              communities=commdetect.output))
  
}

#' Making chord diagram for caberNET results
#'
#'
#' @param  interactome Predicted interactome
#' @param grid.col Optional: color option for cell type
#' @return Chord Diagram highlighting proportion of cell-types in interactome
#' @export
#' 
ChordPlot <- function(interactome, grid.col) {
  chord.format <- mod.cab.filt %>%
    mutate(node1 = paste0(n1cell, "_", ligand)) %>%
    mutate(node2 = paste0(n2cell, "_", receptor))
  
  chord.sigmaxedges <- chord.format %>%
    dplyr::select(node1, node2, n1cell, n2cell) %>%
    unique()
  
  strsplit.ind <- seq(from=1, by=2, length.out = nrow(chord.sigmaxedges))
  adeno.lr <- chord.sigmaxedges %>% 
    mutate(celltypes = paste0(n1cell, "_", n2cell)) %>%
    group_by(celltypes) %>%
    mutate(count = n()) %>%
    ungroup() 
  
  df2 <- adeno.lr %>%
    dplyr::select(n1cell, n2cell, count) %>%
    unique()
  
  # Colors for the different cell types
  #grid.col <- 
  lwd_mat = matrix(1, nrow = nrow(df2), ncol = ncol(df2))
  
  #tiff("adenochord.tiff", width = 5, height = 5, units = 'in', res = 600)
  chordDiagram(df2, 
               grid.col=grid.col, 
               directional=1, 
               direction.type = c("diffHeight", "arrows"), 
               annotationTrack = c("grid"), link.lwd = lwd_mat) 
  
  legend("right", 
         legend = names(grid.col), 
         fill = grid.col,
         border=NA,
         bty="n")
  # dev.off()
}




