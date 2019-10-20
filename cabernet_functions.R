#' Making network function
#' 
#' Given the dataset and genes, make a protein-protein interaction
#' using BioGRID as a base network
#' 
#' @param data Gene expression data
#' @param filt.genes Genes that are nodes in the network
#' @param pathgenes List of KEGG pathway genes TO BE REMOVED
#' @param c Name of cell type to label the genes in the network
#' @return A network
#' @export
#' 
makeNetwork <- function(data, filt.genes, pathgenes, c) {
  rownames(data) <- toupper(rownames(data))
  
  # Selecting pathway genes from BioGRID
  pathgenes <- paste0(c, "_", pathgenes)
  
  pathgenes <- intersect(pathgenes, rownames(data))
  
  temp.biogrid <- g.biogrid
  V(temp.biogrid)$name <- paste0(c, "_", V(temp.biogrid)$name)
  
  net <- induced_subgraph(temp.biogrid, 
                          intersect(setdiff(pathgenes, filt.genes), 
                                    V(temp.biogrid)$name))
  
  edgelist <- as_edgelist(net)
  
  if(nrow(edgelist) == 0) {
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
    #fixclu <- rep(FALSE, vcount(clu.net))
    #fixclu[names(initclu) %in% labelednodes] <- TRUE
    clu.comm <- membership(cluster_label_prop(clu.net, initial=initclu))
    final.comm <- c(final.comm, clu.comm+max(final.comm))
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
calculateDegree <- function(net, communities, dat) {
  deg.comm <- c()
  for(i in unique(communities)){
    comm.net <- induced_subgraph(net, names(communities[communities==i]))
    deg.comm <- c(deg.comm, max(igraph::degree(comm.net)))
  }
  names(deg.comm) <- unique(communities)
  max_num <- max(ncol(dat), 5)
  oversized.comms <- which(deg.comm > max_num)
  print("Degree")
  print(max(deg.comm))
  return(names(oversized.comms))
}

#' Run Graphical Lasso
#' 
#' Applying graphical lasso to each community
#' 
#' @param nodes genes in community
#' @param dat gene expression matrix
#' @param lr ligand receptor pairs
#' @param num community number
#' @return List of glasso results
#' @export 
#' 
runGlasso <- function(nodes, dat, lr, num) {
  lr.mat <- dat[which(rownames(dat) %in% nodes),]
  
  W <- calculateCommGlasso(t(lr.mat), lr)
  lr.W <- W %>% 
    mutate(commnum = num) %>%
    mutate(commsize = length(nodes)) %>%
    dplyr::select(node1, node2, weight, cor, commnum, commsize)
  
  return(lr.W)
}


#' Graphical lasso wrapper to EBICglasso 
#' 
#' 
#' @param x matrix of community gene expressions
#' @param lr list of LR pairs
#' @return A network
#' @export
#' 
calculateCommGlasso <- function(x, lr) {
  n <- nrow(x)
  d <- ncol(x)
  s <- cov(x)
  
  pmat.zero <- data.table(t(combn(seq(1,ncol(x)),2))) %>%
    mutate(A = colnames(x)[V1]) %>%
    mutate(B = colnames(x)[V2]) %>%
    mutate(AB = paste0(A, "_", B)) %>%
    mutate(L = ifelse(A %in% lr$L, 1, 0)) %>%
    mutate(R = ifelse(B %in% lr$L, 1, 0)) %>%
    filter(!(AB %in% c(lr$combo1, lr$combo2))) 
  
  pmat.zero <- pmat.zero %>% dplyr::select(V1, V2)
  
  if(nrow(pmat.zero) == 0) {
    e <- EBICglasso2(s, n=nrow(x), gamma=0, penalizeMatrix=NULL)
    if(nrow(s) == 2) {
      e <- cov2cor(s)
    }
    
  } else {
    e <- EBICglasso2(s, n=nrow(x), gamma=0, 
                     penalizeMatrix=as.matrix(pmat.zero))
  }
  
  ind <- which( upper.tri(e, diag=F) , arr.ind = TRUE )
  W <- matrix(0, nrow(ind),4)
  colnames(W) <- c("node1", "node2", "weight", "cor")
  W[,"node1"] <- colnames(x)[ind[,1]]
  W[,"node2"] <- colnames(x)[ind[,2]]
  W[,"weight"] <- getGlassoEdges(e)
  W[,"cor"] <- getCorEdges(cor(x, method="pearson"))
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

#' EBICglasso function from the EBICglasso package. 
#' Removed check for positive definite matrix
#' 
#' @param S covariance matrix
#' @param n sample number
#' @return EBICglasso output
#' @export
#' 
EBICglasso2 <- function (S, n, gamma = 0.5, penalize.diagonal = FALSE, 
                         nlambda = 100, lambda.min.ratio = 0.001, 
                         returnAllResults = FALSE, checkPD = TRUE, 
                         penalizeMatrix, countDiagonal = FALSE, 
                         refit = FALSE, threshold = FALSE, verbose = TRUE, ...) 
{
  EBICglassoCore2(S = S, n = n, gamma = gamma, penalize.diagonal = penalize.diagonal, 
                  nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, 
                  returnAllResults = returnAllResults, checkPD = checkPD, 
                  penalizeMatrix = penalizeMatrix, countDiagonal = countDiagonal, 
                  refit = refit, ebicMethod = "old", regularized = TRUE, 
                  threshold = threshold, verbose = verbose, ...)
}


EBICglassoCore2 <- function (S, n, gamma = 0.5, penalize.diagonal = FALSE, nlambda = 100, 
                             lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = TRUE, 
                             penalizeMatrix, countDiagonal = FALSE, refit = TRUE, ebicMethod = c("new", 
                                                                                                 "old"), regularized = TRUE, threshold = FALSE, verbose = TRUE, 
                             ...) 
{
  ebicMethod <- match.arg(ebicMethod)
  S <- cov2cor(S)
  lambda.max = max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  lambda.min = lambda.min.ratio * lambda.max
  
  lambda = exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  nlambda <- length(lambda)
  if (missing(penalizeMatrix)) {
    res <- glasso(S, lambda[i], trace=0,
                  penalize.diagonal = penalize.diagonal, ...)
  }
  else {
    glas_path <- list(w = array(0, c(ncol(S), ncol(S), length(lambda))), 
                      wi = array(0, c(ncol(S), ncol(S), length(lambda))), 
                      rholist = lambda)
    for (i in 1:nlambda) {
      res <- glasso(S, lambda[i], trace=0,
                    zero = penalizeMatrix,
                    penalize.diagonal = penalize.diagonal, ...)
      
      glas_path$w[, , i] <- res$w
      glas_path$wi[, , i] <- res$wi
    }
  }
  if (threshold) {
    for (i in 1:nlambda) {
      p <- ncol(glas_path$wi[, , i])
      thresh <- (log(p * (p - 1)/2))/sqrt(n)
      glas_path$wi[, , i] <- ifelse(abs(glas_path$wi[, 
                                                     , i]) < thresh, 0, glas_path$wi[, , i])
    }
  }
  if (ebicMethod == "old") {
    EBICs <- sapply(seq_along(lambda), function(i) {
      if (!regularized) {
        invSigma <- ggmFit(wi2net(glas_path$wi[, , i]), 
                           S, sampleSize = n, ebicTuning = gamma, refit = TRUE, 
                           verbose = FALSE)$invSigma
      }
      else {
        invSigma <- glas_path$wi[, , i]
      }
      EBIC(S, invSigma, n, gamma, countDiagonal = countDiagonal)
    })
  }
  else {
    EBICs <- sapply(seq_along(lambda), function(i) {
      fit <- ggmFit(wi2net(glas_path$wi[, , i]), S, n, 
                    ebicTuning = gamma, refit = !regularized, verbose = FALSE)
      fit$fitMeasures$ebic
    })
  }
  opt <- which.min(EBICs)
  
  net <- as.matrix(forceSymmetric(wi2net(glas_path$wi[, , 
                                                      opt])))
  colnames(net) <- rownames(net) <- colnames(S)
  
  if (refit) {
    if (verbose) 
      message("Refitting network without LASSO regularization")
    if (!all(net[upper.tri(net)] != 0)) {
      glassoRes <- suppressWarnings(glasso::glasso(S, 
                                                   0, zero = which(net == 0 & upper.tri(net), arr.ind = TRUE), 
                                                   trace = 0, penalize.diagonal = penalize.diagonal, 
                                                   ...))
    }
    else {
      glassoRes <- suppressWarnings(glasso::glasso(S, 
                                                   0, trace = 0, penalize.diagonal = penalize.diagonal, 
                                                   ...))
    }
    net <- as.matrix(forceSymmetric(wi2net(glassoRes$wi)))
    colnames(net) <- rownames(net) <- colnames(S)
    optwi <- glassoRes$wi
  }
  else {
    optwi <- glas_path$wi[, , opt]
  }
  
  if (returnAllResults) {
    return(list(results = glas_path, ebic = EBICs, optnet = net, 
                lambda = lambda, optwi = optwi))
  }
  else return(net)
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
getGlassoEdges <- function(fit) {
  fit.net <- sign(fit)
  fit.ind <- which(upper.tri(fit.net,diag=F) , arr.ind=T)
  return(round(fit[fit.ind],2))
}

#' Getting correlation edges given correlation matrix by looking
#' at row vs. column element
#' 
#' @param cor.mat matrix of correlation values
#' @return flattened list of correlation values
#' @export
#' 
getCorEdges <- function(cor.mat) {
  cor.ind <- which(upper.tri(cor.mat, diag=F), arr.ind=T)
  return(round(cor.mat[cor.ind],2))
}
