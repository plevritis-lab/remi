if (!(requireNamespace("org.Hs.eg.db", quietly = TRUE))) {
  BiocManager::install("org.Hs.eg.db")
}
if (!(requireNamespace("clusterProfiler", quietly = TRUE))) {
  BiocManager::install("clusterProfiler")
}
require(tidyverse)
require(glasso)
require(igraph)
require(data.table)
require(OneR)
require(preprocessCore)
require(viridis)
require(clusterProfiler)
require(networkD3)
require(org.Hs.eg.db)
require(msigdbr)
require(dplyr)
require(Seurat)


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
  for(curr.name in cellmarkers) {

    celltype.filt <- dat[,grep(paste0("_", curr.name), colnames(dat))]

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
  
  #filtering non-filtered data
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

  if(filter == F) {
    return(list(filtered = notfiltered.cellexps,
                unfiltered = notfiltered.cellexps))
  } else {
    return(list(filtered=filtered.cellexps,
                unfiltered=notfiltered.cellexps))
  }
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
  if(length(names(datlist)) > 1) {
    for(c in names(datlist)[2:length(datlist)]) {
      if(nrow(datlist[[c]]) == 0) {cat("No genes in ", c)}
      allcellexp <- rbind(allcellexp, datlist[[c]])
    }
  }

  allcellexp <- t(scale(t(allcellexp)))

  allcombolr <- expand.grid(all.ls, all.rs, stringsAsFactors=F)
  first.ind <- seq(from=1, by=2, length.out = nrow(allcombolr))
  second.ind <- seq(from=2, by=2, length.out = nrow(allcombolr))

  lr.network.pairs <- data.table::data.table(allcombolr) %>%
    dplyr::rename(L = Var1) %>%
    dplyr::rename(R = Var2) %>%
    dplyr::mutate(n1cell = unlist(strsplit(L, split="_"))[first.ind]) %>%
    dplyr::mutate(n2cell = unlist(strsplit(R, split="_"))[first.ind]) %>%
    dplyr::mutate(ligand = unlist(strsplit(L, split="_"))[second.ind]) %>%
    dplyr::mutate(receptor = unlist(strsplit(R, split="_"))[second.ind]) %>%
    dplyr::mutate(pairname = paste0(ligand, "_", receptor)) %>%
    dplyr::filter(pairname %in% lr.table$Pair.Name) %>%
    dplyr::mutate(combo1 = paste0(L, "_", R)) %>%
    dplyr::mutate(combo2 = paste0(R, "_", L)) %>%
    dplyr::filter(L %in% all.genes) %>%
    dplyr::filter(R %in% all.genes)

  if(nrow(lr.network.pairs) == 0) {
    cat("No ligand and receptor pairs in dataset. Try relaxing the gene expression cutoff\n")
  }

  return(list(lr.network.pairs=lr.network.pairs, expanddata = allcellexp))
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
generateLRnet <- function(lr.table, celltypes, datlist, cor=T, verbose=T) {

  if(verbose == T) cat("\nStep 1/3: Expanding LR pairs across all cell types\n")
  lr.network.pairs <- expandLRpairs(lr.table, datlist, celltypes)

  if(verbose == T) cat("Step 2/3: Calculating edge weight\n")
  pairwise.cor <- calculateCor(lr.network.pairs)

  if(verbose == T) cat("Step 3/3: Creating graph object\n")
  lr.network <- igraph::graph_from_edgelist(as.matrix(lr.network.pairs$lr.network.pairs[,c(1,2)]),
                                            directed=F)
  igraph::E(lr.network)$weight <- abs(as.numeric(paste0(pairwise.cor$pairwiseLR$cor)))

  return(list(net=lr.network, mat=lr.network.pairs$lr.network.pairs,
              pairwisecor=pairwise.cor, expanddata=lr.network.pairs$expanddata))
}



#' calculateCor - calculating pairwise correlation for all the LR pairs
#'
#' This is an internal function used in generateLRnet. It calculates
#' Pearson correlation of all ligand receptor pairs.
#'
#' @param  lr.table List of cell type-specific LR pairs
#' @return A table of correlations for each LR pair
#' @export
calculateCor <- function(lr.obj, randomize=F, l.chosen=NULL, r.chosen=NULL, T_star=0) {

  allcellexp <- lr.obj$expanddata
  lr.table <- lr.obj$lr.network.pairs

  lr.genes <- unique(c(lr.table$L, lr.table$R))

  lr.allcellexp <- allcellexp[which(rownames(allcellexp) %in% lr.genes),]

  pairwise.lr <- matrix(NA, nrow=nrow(lr.table), ncol=3)

  lr.cor <- cor(t(lr.allcellexp))

  if(randomize == T) {
    lr.cor[l.chosen, r.chosen] <- T_star
    lr.cor[r.chosen, l.chosen] <- T_star
  }

  for(i in 1:nrow(lr.table)) {
    ligand <- lr.table$L[i]
    receptor <- lr.table$R[i]
    cor.res <- lr.cor[which(rownames(lr.cor) == ligand),
                      which(colnames(lr.cor) == receptor)]
    pairwise.lr[i,3] <- cor.res
  }

  pairwise.lr[,1] <- lr.table$L
  pairwise.lr[,2] <- lr.table$R
  pairwise.lr <- data.frame(pairwise.lr)
  colnames(pairwise.lr) <- c("L", "R", "cor")
  pairwise.lr$cor <- as.numeric(paste0(pairwise.lr$cor))

  return(list(pairwiseLR=pairwise.lr, expanddata=allcellexp, lr.cor = lr.cor))
}




#' Calculating eigenvector centrality for each receptor given downstream
#' protein protein-interaction signaling network.
#'
#' @param  lrnet List of cell type-specific LR pairs
#' @param dat.list List of gene expression for each cell type
#' @param pgenelist TO BE REMOVED
#' @param numgenes Number of genes user wants to seed clustering algorithm with
#' @return List of receptors with high scoring eigenvector centrality measurements
#' @export
#'
pickECgenes <- function(lrnet, dat.list, pgenelist, numgenes, cutoff, ppi, seed=10, verbose=T) {

  lr.table <- lrnet$mat
  allcellexp <- lrnet$expanddata

  set.seed(seed)
  pgenes <- unlist(pgenelist)

  if(verbose == T) cat("Building downstream PPI networks\n")
  cell.nets <- list()
  for(c in names(dat.list)) {
    receive.exp <- t(scale(t(dat.list[[c]])))
    receive.net <- makePPINetwork(receive.exp, unique(lr.table$L), pgenes, c, lr.table, ppi)
    cell.nets[[c]] <- receive.net
  }

  if(verbose == T) cat("Calculating importance score\n")

  cell.ec.list <- list()
  cell.ec.all.list <- list()
  all.ec <- c()
  for(c in names(cell.nets)) {
    cell.ec <- igraph::eigen_centrality(cell.nets[[c]])$vector
    cell.ec.r <- cell.ec[which(names(cell.ec) %in% lr.table$R)]
    cell.ec.list[[c]] <- cell.ec.r
    cell.ec.all.list[[c]] <- cell.ec
    all.ec <- c(all.ec, cell.ec)
  }

  if(verbose == T) cat("Selecting important receptors\n")

  sig.top.r.path <- c()
  for(c in names(cell.nets)) {
    plength <- c()
    for(p1 in 1:length(pgenelist)) {
      cell.path.all <- paste0(c, "_", pgenelist[[p1]])
      singler.path <- unique(intersect(cell.path.all, lr.table$R))

      top.genes <- head(sort(cell.ec.list[[c]][singler.path], decreasing=T),
                        numgenes)
      sig.top.r.path <- c(sig.top.r.path, top.genes)
    }
  }

  if(length(sig.top.r.path) == 0) {
    return(NULL)
  }

  sig.path.exp <- allcellexp[which(rownames(allcellexp) %in%
                                     unique(names(sig.top.r.path))),]

  sig.path.high <- sig.path.exp[which(rowMeans(sig.path.exp) > cutoff),]

  return(list(eclist=sig.top.r.path[rownames(sig.path.high)],
              net=cell.nets,
              allec = all.ec))
}




#' remifiedCommunities
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
remifiedCommunities <- function(net, dat.list,
                                seed, verbose=T,
                                cd = "Louvain", maxNum) {
  set.seed(seed)

  # Final list of communities
  lr.communities <- c()

  if(cd == "Louvain") {
      density.comms <- igraph::membership(cluster_louvain(net))
  }

  # Calculate degree of each network and check how many are
  num.oversized <- calculateOversizedComms(net,
                                           density.comms,
                                           dat.list[[1]],
                                           maxNum = maxNum)
  
  # Iterate through all the communities until their
  # degrees match the sample size
  old.comms <- density.comms
  old.oversized <- num.oversized

  if(length(num.oversized) == 0) {
    size.communities <- table(density.comms)
    commnames <- unique(names(size.communities[size.communities > 1]))
    return(list(names=commnames, membership=density.comms))
  }

  while(length(num.oversized) > 0) {
    louvain.comms <- clusterLouvain(net, num.oversized, old.comms)
    num.oversized <- calculateOversizedComms(net,
                                             louvain.comms,
                                             dat.list[[1]],
                                             maxNum = maxNum)

    if(length(setdiff(num.oversized, old.oversized)) == 0) {
      break
    }
    old.comms <- louvain.comms
    old.oversized <- num.oversized
  }

  size.communities <- table(louvain.comms)
  commnames <- unique(names(size.communities[size.communities > 1]))

  if(verbose == T) cat(paste0(length(commnames),
                              " communities identified\n"))

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
remifiedGlasso <- function(netlist, communities, dat.list, seednum, lambda, scale=F) {
  set.seed(seednum)
  predicted.edges <- list()

  allcellexp <- netlist$expanddata

  gedges_w <- matrix("", nrow=0, ncol=9)
  colnames(gedges_w) <- c("node1", "node2", "weight", "cor",
                          "commnum", "commsize", "deg", "numgenes", "lambda")

  cat("Removing conditionally independent edges within communities\n")
  # Graphical Lasso within edges
  pbar <- dplyr::progress_estimated(length(communities$names), min_time = 0)

  for(i in 1:length(communities$names)) {

    communitynum <- communities$names[i]
    commgenes <- names(communities$membership)[communities$membership == communitynum]

    lr.mat <- t(allcellexp[which(rownames(allcellexp) %in% commgenes),])
    lr.comm.net <- igraph::induced_subgraph(netlist$net, colnames(lr.mat))

    W <- calculateCommGlasso(cor(lr.mat), lr.mat, netlist, lambda.max=0.9,
                             lambda=lambda, scale=scale)

    lr.W <- data.table::data.table(W$W) %>%
      dplyr::mutate(commnum = communitynum) %>%
      dplyr::mutate(commsize = length(commgenes)) %>%
      dplyr::select(node1, node2, weight, cor, commnum, commsize, deg, numgenes, lambda)
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

    adjcomms <- adjcomms[which(as.numeric(adjcomms) > as.numeric(comm1.num))]

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
          dplyr::filter((comm1 == comm1.num & comm2 == comm2.num) |
                          (comm1 == comm2.num & comm2 == comm1.num)) %>%
          dplyr::filter(comm1 != comm2)

        betweengenes <- unique(c(betweencommedges$name1, betweencommedges$name2))

        if(length(betweengenes) > 2) {

          lr.mat <- t(allcellexp[which(rownames(allcellexp) %in% betweengenes),])
          W <- calculateCommGlasso(cor(lr.mat), lr.mat, netlist, lambda=lambda, scale=scale, additional=combo)$W
          betweenlr.W <- data.table::data.table(W) %>%
            dplyr::mutate(commnum = paste0(comm1.num, "_", comm2.num)) %>%
            dplyr::mutate(commsize = length(betweengenes)) %>%
            dplyr::select(node1, node2, weight, cor, commnum,
                          commsize, deg, numgenes, lambda)

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
#' @param  obj REMI object
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

  net.edges <- data.table::data.table(switched.W) %>%
    dplyr::mutate(n1cell = unlist(strsplit(node1, split="_"))[strsplit.ind]) %>%
    dplyr::mutate(n2cell = unlist(strsplit(node2, split="_"))[strsplit.ind]) %>%
    dplyr::mutate(ligand = unlist(strsplit(node1, split="_"))[strsplit.ind2]) %>%
    dplyr::mutate(receptor = unlist(strsplit(node2, split="_"))[strsplit.ind2]) %>%
    dplyr::select(node1, node2, ligand, receptor, n1cell, n2cell,
                  weight, cor, commsize, commnum, deg, numgenes, lambda) %>%
    dplyr::mutate(pairname = paste0(node1, "_", node2)) %>%
    dplyr::filter(pairname %in% netlist$mat$combo1) %>%
    dplyr::group_by(node1, node2) %>%
    unique() %>%
    dplyr::group_by(pairname) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup()

  commIDs <- unique(net.edges$commnum)

  between.comms <- commIDs[grep("_", commIDs)]

  if(length(between.comms) > 0) {
    within.comms <- commIDs[-grep("_", commIDs)]
  } else {
    within.comms <- commIDs
  }

  within.edges <- net.edges %>%
    dplyr::filter(commnum %in% within.comms)

  between.edges <- net.edges %>%
    dplyr::filter(commnum %in% between.comms)

  between.edges <- between.edges %>%
    dplyr::filter(!(pairname %in% within.edges$pairname))

  net.edges.filt <- bind_rows(within.edges, between.edges)

  filtered.net.edges <- net.edges.filt %>%
    dplyr::filter(abs(as.numeric(weight)) > 0)

  return(list(net.edges=net.edges.filt,
              filtered.net.edges=filtered.net.edges))
}


#' REMI algorithm
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
remi <- function(cellmarkers, dat.list, seed=30,
                 lambda=NULL, lr.database=NULL,
                 downstreamgenes=NULL, ppi.net=NULL,
                 cd = "Louvain", maxNum = NULL) {

  if(is.null(lr.database)) {lr.database = curr.lr.filt}
  if(is.null(downstreamgenes)) {downstreamgenes = pathway.genelist}
  if(is.null(ppi.net)) {ppi.net = g.biogrid}

  # Creating LR network
  cat("Building LR network\n")
  pathwaygenes <- unique(unlist(pathway.genelist$genesets))
  netlist <- generateLRnet(lr.database, cellmarkers,
                           dat.list$filtered, downstreamgenes)

  # Cluster
  cat("Detecting communities")
  commdetect.output <- remifiedCommunities(netlist$net,
                                           dat.list$filtered,
                                           seed=seed,
                                           cd=cd,
                                           maxNum = maxNum)

  # Graphical Lasso
  cat("\n Estimating activated LR pairs \n")
  glasso.output <- remifiedGlasso(netlist, commdetect.output,
                                  dat.list$filtered, seed,
                                  lambda, scale=F)

  predicted.edges <- cleaningOutput(glasso.output, netlist)

  numLRs <- length(unique(predicted.edges$filtered.net.edges$node1,
                          predicted.edges$filtered.net.edges$node2))

  if(nrow(predicted.edges$filtered.net.edges) == 0) {
    cat("\nNo interactome found\n")
    return(NULL)
  }

  params <- list()
  params[["seed"]] <- seed

  return(list(lrnet=netlist,
              unfilteredinteractome=predicted.edges$net.edges,
              interactome=predicted.edges$filtered.net.edges,
              communities=commdetect.output,
              params=params))

}






#' Making chord diagram for REMI results
#'
#'
#' @param interactome Predicted interactome
#' @param grid.col Optional: color option for cell type
#' @return Chord Diagram highlighting proportion of cell-types in interactome
#' @export
#'
REMIPlot <- function(interactome, type="chord",
                     grid.col=NULL, size=10, thres=0,
                     selectcell = NULL,
                     legend = FALSE) {
  
  chord.format <- interactome$interactome %>%
    filter(weight > thres) %>%
    dplyr::mutate(node1 = paste0(n1cell, "_", ligand)) %>%
    dplyr::mutate(node2 = paste0(n2cell, "_", receptor))
  
  if(!(is.null(selectcell))) {
    chord.format <- chord.format %>%
      filter(n1cell == selectcell | n2cell == selectcell)
  }
  
  chord.sigmaxedges <- chord.format %>%
    dplyr::select(node1, node2, n1cell, n2cell) %>%
    unique()
  
  strsplit.ind <- seq(from=1, by=2, length.out = nrow(chord.sigmaxedges))
  adeno.lr <- chord.sigmaxedges %>%
    dplyr::mutate(celltypes = paste0(n1cell, "_", n2cell)) %>%
    dplyr::group_by(celltypes) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup()

  df2 <- adeno.lr %>%
    dplyr::select(n1cell, n2cell, count) %>%
    unique()

  # Colors for the different cell types
  lwd_mat = matrix(1, nrow = nrow(df2), ncol = ncol(df2))

  if(type == "chord") {
    circlize::chordDiagram(df2,
                           grid.col=grid.col,
                           directional=1,
                           direction.type = c("diffHeight"),
                           annotationTrack = c("grid"), link.lwd = lwd_mat)

    if(legend == TRUE) {
      legend("right",
             legend = names(grid.col),
             fill = grid.col,
             border=NA,
             bty="n")
    }
  }

  if(type == "alluvial") {
    p <- ggplot(df2,
           aes(y = count, axis1 = n1cell, axis2 = n2cell)) +
      geom_alluvium(aes(fill = n1cell), width = 1/12) +
      geom_stratum(width = 1/12) +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=size) +
      scale_x_discrete(limits = c("Signaling", "Receiving"), expand = c(.05, .05)) +
      scale_fill_manual(values=grid.col) +
      coord_flip() +
      theme_void()

    return(p)
  }
}


#' Calculating p-value for a given edge in a community
#'
#'
#' @param obj REMI object
#' @param l.chosen Ligand in edge
#' @param r.chosen Receptor in edge
#' @param comm.chosen Community that edge is present in
#' @param iterNum = Number of permutations for p-value calculation
#' @param seednum = Seed
#' @param lambda = Manually setting lambda if of interest
#' @return Chord Diagram highlighting proportion of cell-types in interactome
#' @export
#'
calculateSignificance <- function(obj,
                                  l.chosen,
                                  r.chosen,
                                  comm.chosen,
                                  maxNum,
                                  iterNum = 1000,
                                  seednum=30,
                                  lambda = NULL) {

  # Setting Varibales
  allcellexp <- obj$lrnet$expanddata
  node1 <- obj$interactome$node1
  node2 <- obj$interactome$node2
  commnums <- obj$interactome$commnum
  param <- obj$params
  between <- FALSE

  if(is.null(lambda)) {
    opt.lambda <- obj$interactome %>%
      dplyr::filter(commnum == comm.chosen) %>%
      dplyr::select(lambda) %>%
      pull() %>% unique()
  } else {
    opt.lambda <- lambda
  }

  if("_" %in% comm.chosen) {
    between = TRUE
  }

  # Obtaining community genes
  comm.num.index <- which(commnums == comm.chosen)

  orig.comm.genes <- unique(c(node1[comm.num.index],
                              node2[comm.num.index]))

  cat("Number of genes in community\t", length(orig.comm.genes),"\n")

  if(length(orig.comm.genes) > ncol(allcellexp)) {
    cat("Can not calculate significance due to community size.\n")
    return(NULL)
  }

  # Parameters for original community
  community.mat <- t(allcellexp[which(rownames(allcellexp) %in% orig.comm.genes),])
  R <- cor(community.mat)
  S <- cov(community.mat)
  l.ind <- which(colnames(R) == l.chosen)
  r.ind <- which(colnames(R) == r.chosen)

  if(length(l.ind) == 0) { cat("LR pair not in community") }

  comm.change <- 0

  cat("Calculating p-value\n")
  pbar <- dplyr::progress_estimated(iterNum, min_time = 0)

  pval.list <- c()
  T.list <- c()
  Y.list <- c()
  for(i in 1:iterNum) {
    set.seed(i*seednum)

    T_star = 2 * runif(1) - 1 # a uniform value in [-1,1]

    R_star_cor_df <- calculateCor(list(expanddata=allcellexp,
                                       lr.network.pairs=obj$lrnet$mat),
                                  randomize=T,
                                  l.chosen=l.chosen,
                                  r.chosen=r.chosen,
                                  T_star=T_star)

    boot.net <- graph_from_edgelist(as.matrix(R_star_cor_df$pairwiseLR[,c("L", "R")]), directed=F)
    E(boot.net)$weight <- abs(as.numeric(R_star_cor_df$pairwiseLR[,"cor"]))

    commboot.output <- remifiedCommunities(net=boot.net,
                                           dat.list=cellexp.list$filtered,
                                           seed=param$seed,
                                           verbose=F,
                                           maxNum = maxNum)

    commnums.boot <- commboot.output$membership

    l.boot.num <- commnums.boot[l.chosen]
    r.boot.num <- commnums.boot[r.chosen]

    if(between) {
      if(l.boot.num == r.boot.num) {
        T.list <- c(T.list, T_star)
        Y.list <- c(Y.list, 0)
        comm.change <- comm.change + 1
        next
      }

      comm1genes <- names(commnums.boot)[which(commnums.boot == l.boot.num)]
      comm2genes <- names(commnums.boot)[which(commnums.boot == r.boot.num)]
      commgenes <- unique(c(comm1genes, comm2genes))

      between.edges <- findAdjacentCommEdges(obj$lrnet, commgenes, commboot.output) %>%
        filter((comm1 == l.boot.num & comm2 == r.boot.num) |
                 (comm1 == r.boot.num & comm2 == l.boot.num)) %>%
        filter(comm1 != comm2)

      lr.boot.genes <- unique(c(between.edges$name1, between.edges$name2))
    } else { lr.boot.genes <- names(commnums.boot[commnums.boot == l.boot.num]) }

    comm.mat <- t(allcellexp[which(rownames(allcellexp) %in% lr.boot.genes),])
    comm.cor <- R_star_cor_df$lr.cor[which(rownames(R_star_cor_df$lr.cor) %in% lr.boot.genes),
                                     which(colnames(R_star_cor_df$lr.cor) %in% lr.boot.genes)]

    g <- calculateCommGlasso(comm.cor, comm.mat, obj$lrnet,
                             lambda=opt.lambda,
                             scale=F)$W %>%
                             dplyr::mutate(weight = if_else(abs(as.numeric(weight)) > 0, 1, 0))

    node.ind <- which(g$node1 == l.chosen & g$node2 == r.chosen)
    if(length(node.ind) == 0) {
      node.ind <- which(g$node2 == l.chosen & g$node1 == r.chosen)
    }

    Y_star <- g$weight[node.ind]

    T.list <- c(T.list, T_star)
    Y.list <- c(Y.list, Y_star)

    pbar$pause(0.01)$tick()$print()
  }

  D <- list(T=T.list, Y=Y.list)

  pval <- calculatePvalue(R, S, D, i_=l.ind, j_=r.ind,
                          n=ncol(allcellexp),
                          p=length(orig.comm.genes))

  return(list(pval=pval, D=D))
}

#' Creating single cell REMI object
#'
#'
#' @param obj REMI object
#' @param celltype.col what cell types to measure
#' @param remove.markers cell types that should not be used for calculations
#' @param gene.selct genes of interest
#' @param assay = Seurat assay
#' @param filter = remove low expressed ligand and receptors
#' @param thres = average assay expression filter
#' @return Single cell REMI object that can be used by the algorithm
#' @export
#'
setupSingleCell <- function(obj, sample.col,
                            celltype.col,
                            remove.markers = NULL,
                            gene.select = NULL,
                            assay="integrated",
                            filter=T,
                            thres=0,
                            expthres = 0.1) {

  cat("Calculating percent expressed for ligand and receptor genes\n")
  
  Idents(obj) <- celltype.col
  d <- Seurat::DotPlot(obj, features=rownames(obj))
  percexp <- d$data %>%
    filter(pct.exp > expthres) %>%
    dplyr::mutate(cell_gene = paste0(id, "_", features.plot))
  
  cat("Averaging expression\n")
  
  pseudobulk <- SingleToBulk(obj, assay, sample.col, celltype.col)

  num.markers <- length(pseudobulk$cellmarkers) - length(remove.markers)
  print(num.markers)

  filtered.cellexps <- list()
  notfiltered.cellexps <- list()

  colsvec <- c()

  for(i in 1:length(pseudobulk$cellmarkers)) {

    curr.name <- pseudobulk$cellmarkers[i]

    if(!(curr.name %in% remove.markers)) {

      nospace.name <- as.character(gsub(" ", "", curr.name))

      all.cols <- unlist(lapply(colnames(pseudobulk$dat),
                                function(x) {strsplit(x, "_")[[1]][2]}))
      all.cols <- gsub(" ", "", all.cols)

      cell.cols <- grep(paste0("^\\b", nospace.name, "\\b$"), all.cols)
      
      celltype.filt <- pseudobulk$dat[,cell.cols]

      # Match sample name
      dat.cols <- unlist(lapply(colnames(celltype.filt),
                                function(x) {strsplit(x, "_")[[1]][1]}))
      colsvec <- c(colsvec, dat.cols)

      # Removing genes with low expression
      removegenes <- rownames(celltype.filt)[which(rowMeans(celltype.filt) <= thres)]

      # Finalized cleaned data matrix
      if(length(removegenes) > 0) {
        cleaned.filt <- pseudobulk$scaleddat[-which(rownames(pseudobulk$scaleddat) %in%
                                                      removegenes), cell.cols]
      } else {
        cleaned.filt <- pseudobulk$scaleddat[,cell.cols]
      }

      # Adding in cell type into rownames
      rownames(celltype.filt) <- paste0(curr.name, "_", rownames(celltype.filt))
      colnames(celltype.filt) <- dat.cols

      rownames(cleaned.filt) <- paste0(curr.name, "_", rownames(cleaned.filt))
      colnames(cleaned.filt) <- dat.cols
      
      cleaned.filt <- cleaned.filt[which(rownames(cleaned.filt) %in% percexp$cell_gene),]

      if(is.null(gene.select)) {
        filtered.cellexps[[curr.name]] <- cleaned.filt
        notfiltered.cellexps[[curr.name]] <- celltype.filt
      } else {
        filtered.cellexps[[curr.name]] <- cleaned.filt[which(rownames(cleaned.filt) %in%
                                                               gene.select),]
        notfiltered.cellexps[[curr.name]] <- celltype.filt[which(rownames(celltype.filt) %in%
                                                                   gene.select),]

      }
    }
  }
  
  allpresent.cols <- names(which(table(colsvec) == num.markers))
  
  if(length(allpresent.cols) == 1) {
    cat("Please remove cell types. Only one patient with all cell types\n")
    stop
  }
  
  #Filtering non-filtered data
  for(c in names(notfiltered.cellexps)) {
    
    allpresent.filt.data <- notfiltered.cellexps[[c]][,which(colnames(notfiltered.cellexps[[c]]) %in%
                                                              allpresent.cols)]
    
    if(is.null(nrow(allpresent.filt.data))) {
      cat(c, " was omitted due to lack of uniform expression.\n")
    } else {
    
      duplicate.cols <- which(apply(allpresent.filt.data, 1,
                                    function(x) length(unique(x))==1) == TRUE)
  
      if(length(duplicate.cols) > 0) {
        notfiltered.cellexps[[c]] <- allpresent.filt.data[-duplicate.cols,]
      } else {
        notfiltered.cellexps[[c]] <- allpresent.filt.data
      }
    }
  }
  
  # Only using samples that have all the cell types of interest
  for(c in names(filtered.cellexps)) {
    allpresent.filt.data <- filtered.cellexps[[c]][,which(colnames(filtered.cellexps[[c]]) %in%
                                                            allpresent.cols)]
    
    if(is.null(nrow(allpresent.filt.data))) {
      temp <- "hi"
    } else {
      duplicate.cols <- which(apply(allpresent.filt.data, 1,
                                    function(x) length(unique(x)) == 1) == TRUE)
  
      if(length(duplicate.cols) > 0) {
        filtered.cellexps[[c]] <- allpresent.filt.data[-duplicate.cols,]
      } else {
        filtered.cellexps[[c]] <- allpresent.filt.data
      }
    }
  }

  if(filter == F) {
    return(list(filtered=notfiltered.cellexps,
                unfiltered=notfiltered.cellexps,
                cellmarkers=setdiff(pseudobulk$cellmarkers, remove.markers)))
  } else {
    return(list(filtered=filtered.cellexps,
                unfiltered=notfiltered.cellexps,
                cellmarkers=setdiff(pseudobulk$cellmarkers, remove.markers)))
  }
}

#' Averaging across expression levels
#'
#'
#' @param obj REMI object
#' @param assay what Seurat assay to use for calculations
#' @return Average expressiona cross patients
#' @export
#'
SingleToBulk <- function(obj, assay, samplecol, celltypecol) {
  temp.rownames <- rownames(obj@meta.data)
  obj@meta.data <- obj@meta.data %>%
  dplyr::mutate(samplesREMI = gsub("_", "", !!as.name(samplecol))) %>%
  dplyr::mutate(group.ctype = paste0(samplesREMI, "_" , !!as.name(celltypecol))) %>%
  dplyr::mutate(group.ctype = as.factor(group.ctype))
  rownames(obj@meta.data) <- temp.rownames

  Idents(obj) <- "group.ctype"
  avg.obj <- AverageExpression(obj, return.seurat=T, assays=assay)
  
  avg.dat <- GetAssayData(avg.obj, "data") %>% as.matrix
  
  avg.dat[is.na(avg.dat)] <- 0

  avg.scaled <- t(scale(t(avg.dat)))
  avg.scaled <- na.omit(avg.scaled)

  cellmarkers <- unique(unlist(lapply(colnames(avg.scaled),
                                      function(x) {strsplit(x, "_")[[1]][2]})))
  names(cellmarkers) <- cellmarkers

  return(list(scaleddat=avg.scaled, dat=avg.dat, cellmarkers=cellmarkers))
}

