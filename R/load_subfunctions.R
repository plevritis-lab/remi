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
#'
makePPINetwork <- function(data, filt.genes, pathgenes, c, lrmat, ppi=g.biogrid) {
  # Selecting pathway genes from BioGRID
  pathgenes <- intersect(paste0(c, "_", pathgenes), rownames(data))
  pathgenes <- c(pathgenes, lrmat$R)
  igraph::V(ppi)$name <- paste0(c, "_", igraph::V(ppi)$name)

  overlap.genes <- intersect(setdiff(pathgenes, filt.genes),
                             igraph::V(ppi)$name)

  net <- igraph::induced_subgraph(ppi, vids=overlap.genes)
  edgelist <- igraph::as_edgelist(net)

  if(nrow(edgelist) < 5) {
    print(paste0(c, " has less than 5 PPI edges"))
    egraph <- igraph::make_empty_graph(n=0, directed=T)
    return(egraph)
  }

  # Calculating edge weights
  cor.edges <- cor(t(data[which(rownames(data) %in% igraph::V(net)$name),]),
                   method="pearson")
  cor.weights <- rep(0, igraph::ecount(net))

  for(e in 1:igraph::ecount(net)) {
    cor.weights[e] <- abs(cor.edges[edgelist[e,1], edgelist[e,2]])
  }

  igraph::E(net)$weight <- cor.weights
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
#'
findAdjacentCommEdges <- function(netlist, cgenes, comms) {
  all_edges <- igraph::E(netlist$net)[inc(cgenes)]
  all_edges_m <- data.table::data.table(igraph::get.edges(netlist$net, all_edges)) %>%
    dplyr::mutate(name1 = igraph::V(netlist$net)$name[V1]) %>%
    dplyr::mutate(name2 = igraph::V(netlist$net)$name[V2]) %>%
    dplyr::mutate(comm1 = comms$membership[name1]) %>%
    dplyr::mutate(comm2 = comms$membership[name2]) %>%
    dplyr::mutate(pairname = paste0(name1, "_", name2))

  return(all_edges_m)
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
#'
clusterLabelProp <- function(net, clu, clu.labeled, labelednodes) {
  final.comm <- clu[clu %in% setdiff(unique(clu), clu.labeled)]

  for(cc in clu.labeled) {
    clu.net <- igraph::induced_subgraph(net, names(clu[clu==cc]))
    initclu <- seq(0, igraph::vcount(clu.net)-1)
    names(initclu) <- igraph::V(clu.net)$name
    initclu[-which(names(initclu) %in% labelednodes)] <- -1

    fixclu <- rep(FALSE, igraph::vcount(clu.net))
    fixclu[names(initclu) %in% labelednodes] <- TRUE
    clu.comm <- igraph::membership(igraph::cluster_label_prop(clu.net,
                                                              fixed=fixclu,
                                                              initial=initclu))

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
#'
clusterLouvain <- function(net, commnums, communities) {
  suppressWarnings(
  if(unique(commnums) == unique(communities)) {
    cat("\nAll communities are larger than sample size\n")
    final.comms <- NA
  } else {
    final.comms <- communities[-which(communities %in% commnums)]
  })
  for(cc in commnums) {
    comm.net <- igraph::induced_subgraph(net,names(communities[communities==cc]))
    louvain.comms <- igraph::membership(igraph::cluster_louvain(comm.net))
    suppressWarnings(if(is.na(final.comms)) {
      final.comms <- louvain.comms
    } else {
      final.comms <- c(final.comms, louvain.comms+max(final.comms))
    })
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
#'
calculateOversizedComms <- function(net, communities, dat, maxNum) {
  deg.comm <- c()
  for(i in unique(communities)){
    deg.comm <- c(deg.comm,
                  calculateNetVars(net, names(communities[communities==i]))$numnodes)
  }
  names(deg.comm) <- unique(communities)
  if(is.null(maxNum)) { maxNum = ncol(dat) }
  max_num <- max(maxNum, 5)
  oversized.comms <- which(deg.comm > max_num)
  return(names(oversized.comms))
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
#'
calculateNetVars<- function(net, genes) {
  comm.net <- igraph::induced_subgraph(net, genes)
  degree <- max(igraph::degree(comm.net))
  numedges <- igraph::ecount(comm.net)
  numnodes <- igraph::vcount(comm.net)
  return(list(deg=degree, numedges=numedges, numnodes=numnodes))
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
#'
calculateCommGlasso <- function(S, x, netlist, lambda.max = 0.9,
                                lambda = NULL,
                                scale=F,
                                additional=NULL, verbose=T) {
  net <- netlist$net
  lr <- netlist$mat

  n <- nrow(x)
  d <- ncol(x)
  if(scale == T) {
    x <- scale(x)
    S <- cor(x)
  }

  allLRpairs <- c(lr$combo1, lr$combo2)
  filtallLRpairs <- setdiff(allLRpairs, additional)

  pmat.all <- data.table::data.table(t(combn(seq(1,ncol(x)),2))) %>%
    dplyr::mutate(A = colnames(x)[V1]) %>%
    dplyr::mutate(B = colnames(x)[V2]) %>%
    dplyr::mutate(AB = paste0(A, "_", B))

  pmat.zero <- pmat.all %>%
    dplyr::filter(!(AB %in% filtallLRpairs)) %>%
    dplyr::select(V1, V2) %>%
    as.matrix()

  if(nrow(pmat.zero) == 0) {
    pmat.zero <- NULL
  }

  if(is.null(lambda)) {
    bic <- c()
    lambda.max = max(max(S - diag(d)), -min(S - diag(d)))
    lambda.min = 0.1 * lambda.max
    lambdas = exp(seq(log(lambda.max), log(lambda.min), length = 20))

    #lambda.min <- 0.1
    #lambdas = exp(seq(log(lambda.min), log(lambda.max), length = 20))
    for(l in lambdas) {
      res <- glasso::glasso(S,
                            l,
                            zero = pmat.zero,
                            approx=F)

      bic <- c(bic, EBIC(S, res$wi, gamma=0, n))
    }

    # Pick optimal lambda
    opt.lambda <- lambdas[which.min(bic)]
  } else {
    opt.lambda <- as.numeric(lambda)
  }

  e <- glasso::glasso(S, opt.lambda,
                      zero=pmat.zero,
                      approx=F)

  pred.edges <- qgraph::wi2net(e$wi)
  #pred.edges <- e$wi

  net.vars <- calculateNetVars(net, colnames(x))

  # Output
  ind <- which( upper.tri(pred.edges, diag=F) , arr.ind = TRUE )
  W <- matrix(0, nrow(ind), 7)
  colnames(W) <- c("node1", "node2", "weight", "cor", "deg", "numgenes", "lambda")
  W[,"node1"] <- colnames(x)[ind[,1]]
  W[,"node2"] <- colnames(x)[ind[,2]]
  W[,"weight"] <- getUpperTri(pred.edges, round = TRUE)
  W[,"deg"] <- net.vars$deg
  W[,"numgenes"] <- net.vars$numnodes
  W[,"cor"] <- getUpperTri(cor(x, method="pearson"), round = TRUE)
  W[,"lambda"] <- opt.lambda

  return(list(W=data.frame(W, stringsAsFactors = F), S=S, lambda=opt.lambda))
}

#' Log Likelihood
#'
#' @param data Community gene expression matrix
#' @param theta Parameter for log-likelihood calculation
#' @return A network
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
#'
EBIC <- function (S, K, n, E, gamma = 1, countDiagonal = FALSE)
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

#' Generate null density given covariance matrix
#'
#' @param r density
#' @param S covariance matrix
#' @param i index for ligand
#' @param j index for receptor
#' @param n sample size
#' @return Null density
#'
null_density = function(r, S, i, j, n) {

  p = nrow(S)
  S_0 = 1 * S
  S_0[i,j] = 0
  S_0[j,i] = 0
  Theta_0 = solve(S_0)

  a = Theta_0[i,j]
  b = Theta_0[i,i]
  c = Theta_0[j,j]
  d = sqrt(S[i,i] * S[j,j])

  V = ((1 + r * d * a)^2 - r^2 * d^2 * b * c)
  V[V<0] = 0

  return(V^((n-p-2)/2))
}

#' Calculates p-value for a community
#'
#' @param R_ correlation matrix
#' @param S_ covariance matrix
#' @param D_ list of correlation perturbations and number of times graphical lasso predicted an edge to be present or not
#' @param i_ index for ligand
#' @param j_ index for receptor
#' @param n sample size
#' @param p number of features in community
#' @return p-value
#'
calculatePvalue <- function(R_, S_, D_, i_, j_, n, p) {
  if(p > n) {cat("More features in community than sample size. Can not calculate p-values\n")}

  clf_ = tree::tree(factor(Y) ~ T, data=D_)
  null_sample = rnorm(10000, sd=2 / sqrt(n-p))
  null_sample = dnorm(null_sample, sd=2/sqrt(n-p))

  weights_ = (as.numeric(predict(clf_, newdata=data.frame(T=null_sample),
                                 type='class')) - 1)

  weights_ = weights_ * (null_density(null_sample, S_, n=n, i=i_, j=j_)
                         / dnorm(null_sample, sd=2 / sqrt(n-p)))

  onesided_pvalue = (mean(weights_ * (null_sample > R_[i_, j_]))  + 1)/(mean(weights_)+ 1)

  twosided_pvalue = 2 * min(onesided_pvalue, 1 - onesided_pvalue)

  return(list(p=twosided_pvalue, w=weights_, w_n=(weights_ * (null_sample > R_[i_, j_]))))
}

