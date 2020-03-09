#' Making network function
#' 
#' Given the dataset and genes, make a protein-protein interaction network
#' using BioGRID as a base network
#' 
#' @param data Gene expression data
#' @param filt.genes Genes that are nodes in the network
#' @param pathgenes List of KEGG pathway genes TO BE REMOVED
#' @param c Name of cell type to label the genes in the network
#' @return A network
#' @export
#' 
makeNetwork <- function(data, filt.genes, pathgenes, c, ppi=g.biogrid) {
  #rownames(data) <- toupper(rownames(data))
  
  # Selecting pathway genes from BioGRID
  pathgenes <- paste0(c, "_", pathgenes)
  pathgenes <- intersect(pathgenes, rownames(data))
  
  V(ppi)$name <- paste0(c, "_", V(ppi)$name)
  
  net <- induced_subgraph(ppi, intersect(setdiff(pathgenes, filt.genes), 
                                    V(ppi)$name))
  
  edgelist <- as_edgelist(net)
  
  if(nrow(edgelist) == 0) {
    print("No PPI edges")
    return(NULL)
  }
  
  # Calculating edge weights
  cor.edges <- cor(t(data[which(rownames(data) %in% V(net)$name),]), 
                   method="pearson")
  cor.weights <- rep(0, ecount(net))
  
  for(e in 1:ecount(net)) {
    cor.weights[e] <- abs(cor.edges[edgelist[e,1], edgelist[e,2]])
  }
  
  E(net)$weight <- cor.weights
  net
}

#' Clustering using label propagation
#' 
#' This is the initial clustering step in cabernet. The network is
#' seeded with receptors that had a high eigenvector centrality
#' measurement to produce clusters that reflect downstream signaling
#' activity seen in the dataset.
#' 
#' @param net LR network
#' @param clu Components in the network
#' @param clu.labeled Communities that contain at least two EC receptors
#' @param labelednodes List of EC receptors
#' @return Communities identified using label propagation
#' @export
#' 
clusterLabelProp <- function(net, clu, clu.labeled, labelednodes) {
  final.comm <- clu[clu %in% setdiff(unique(clu), clu.labeled)]
  
  for(cc in clu.labeled) {
    clu.net <- induced_subgraph(net, names(clu[clu==cc]))
    initclu <- seq(0, vcount(clu.net)-1)
    names(initclu) <- V(clu.net)$name
    initclu[-which(names(initclu) %in% labelednodes)] <- -1
    
    fixclu <- rep(FALSE, vcount(clu.net))
    fixclu[names(initclu) %in% labelednodes] <- TRUE
    clu.comm <- membership(cluster_label_prop(clu.net, initial=initclu))
    
    if(length(final.comm) != 0) {
      final.comm <- c(final.comm, clu.comm+max(final.comm))
    } else {
      final.comm <- clu.comm
    }
  }
  
  return(final.comm)
}

#' Clustering using Louvain
#' 
#' The large communities identified by label propagation are broken down 
#' into even smaller communities to ensure the degree of the networks
#' reflect the size of the dataset to maximize significance.
#' 
#' @param net LR network
#' @param commnums List of communities that had a degree larger than sample size
#' @param communities List of all communities in LR network
#' @return Clusters identified using louvain combined with all communities
#' @export
#' 
clusterLouvain <- function(net, commnums, communities) {
  final.comms <- communities[-which(communities %in% commnums)]
  for(cc in commnums) {
    comm.net <- induced_subgraph(net,names(communities[communities==cc]))
    louvain.comms <- membership(cluster_louvain(comm.net))
    final.comms <- c(final.comms, louvain.comms+max(final.comms))
  }
  return(final.comms)
}

#' Calculating the degree of each community in a network
#' 
#' Used to assess significance of the network relative to the
#' sample size of the dataset
#' 
#' @param net LR network
#' @param communites Communities of interest
#' @param dat Dataset of one cell type
#' @return List of maximum degree in each community
#' @export
#' 
calculateOversizedComms <- function(net, communities, dat) {
  deg.comm <- c()
  for(i in unique(communities)){
    deg.comm <- c(deg.comm, 
                  calculateNetVars(net, 
                                   names(communities[communities==i]))$deg)
  }
  names(deg.comm) <- unique(communities)
  max_num <- max(ncol(dat), 10)
  oversized.comms <- which(deg.comm > max_num)
  return(names(oversized.comms))
}

calculateNetVars<- function(net, genes) {
  comm.net <- induced_subgraph(net, genes)
  degree <- max(igraph::degree(comm.net))
  numedges <- vcount(comm.net)
  return(list(deg=degree, numedges=numedges))
}

calculateCommGlasso <- function(x, netlist) {
  net <- netlist$net
  lr <- netlist$mat
  
  n <- nrow(x)
  d <- ncol(x)
  #x <- scale(x) * sqrt((n - 1)/n)
  S <- cov(x)
  
  pmat.zero <- data.table(t(combn(seq(1,ncol(x)),2))) %>%
    mutate(A = colnames(x)[V1]) %>%
    mutate(B = colnames(x)[V2]) %>%
    mutate(AB = paste0(A, "_", B)) %>%
    mutate(L = ifelse(A %in% lr$L, 1, 0)) %>%
    mutate(R = ifelse(B %in% lr$L, 1, 0)) %>%
    filter(!(AB %in% c(lr$combo1, lr$combo2))) 
  
  pmat.zero <- as.matrix(pmat.zero %>% dplyr::select(V1, V2))
  
  if(nrow(pmat.zero) == 0) {
    pmat.zero <- NULL
  }
  
  bic <- c()
  lambda.max <- 0.9
  lambda.min <- 0.1
  lambdas = exp(seq(log(lambda.min), log(lambda.max), length = 20))
  for(l in lambdas) {
    res <- glasso(S, l, zero = pmat.zero, thr=1e-08)
    bic <- c(bic, EBIC(S, res$wi, n, 0, countDiagonal = F))
  }
  
  # Pick optimal lambda
  opt.lambda <- lambdas[which.min(bic)]
  
  e <- glasso(S, opt.lambda, zero=pmat.zero, thr=1e-08)
  
  #pred.edges <- as.matrix(wi2net(e$wi))
  pred.edges <- wi2net(e$wi)
  
  net.vars <- calculateNetVars(net, colnames(x))
  
  # Output
  ind <- which( upper.tri(pred.edges, diag=F) , arr.ind = TRUE )
  W <- matrix(0, nrow(ind), 6)
  colnames(W) <- c("node1", "node2", 
                   "weight", 
                   "cor", "deg", 
                   "numedges")
  W[,"node1"] <- colnames(x)[ind[,1]]
  W[,"node2"] <- colnames(x)[ind[,2]]
  W[,"weight"] <- getUpperTri(pred.edges, round = TRUE)
  W[,"deg"] <- net.vars$deg
  W[,"numedges"] <- net.vars$numedges
  W[,"cor"] <- getUpperTri(cor(x, method="pearson"), round = TRUE)
  return(data.frame(W, stringsAsFactors = F))
}

#' Log Likelihood
#' 
#' @param data Community gene expression matrix
#' @param theta Parameter for log-likelihood calculation
#' @return A network
#' @export
#' 
loglik_ave <- function(data, theta){
  return(-(log(det(theta)) - sum(diag(var(data) %*% theta))))
}

#' EBIC calculation
#' 
#' @param S Covariance matrix 
#' @param K Inverse covariance matrix
#' @param n Sample number
#' @return EBIC number
#' @export
#' 
EBIC <- function (S, K, n, gamma = 0.5, E, countDiagonal = FALSE) 
{
  L <- logGaus(S, K, n)
  if (missing(E)) {
    E <- sum(K[lower.tri(K, diag = countDiagonal)] != 0)
  }
  p <- nrow(K)
  -2 * L + E * log(n) + 4 * E * gamma * log(p)
}

#' Calculating log Gaussian
#' 
#' @param S Covariance matrix
#' @param K Inverse covariance matrix
#' @param n Sample size
#' @return Log Gaussian calculation
#' @export
#' 
logGaus <- function (S, K, n) 
{
  KS = K %*% S
  tr = function(A) sum(diag(A))
  return(n/2 * (log(det(K)) - tr(KS)))
}

#' Given graphical lasso matrix form output, extract edges by 
#' looking at rows vs. column elements
#' 
#' @param fit glasso matrix output
#' @return List form of predicted edges
#' @export
#' 
getUpperTri<- function(fit, round=TRUE) {
  fit.net <- sign(fit)
  fit.ind <- which(upper.tri(fit.net, diag=F) , arr.ind=T)
  if(round == TRUE) {
    return(round(fit[fit.ind],2))
  } else {
    return(fit[fit.ind])
  }
}
