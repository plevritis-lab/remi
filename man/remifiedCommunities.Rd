% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_functions.R
\name{remifiedCommunities}
\alias{remifiedCommunities}
\title{remifiedCommunities}
\usage{
remifiedCommunities(net, dat.list, seed, verbose = T, cd = "Louvain", maxNum)
}
\arguments{
\item{net}{List of LR network and LR table}

\item{dat.list}{List of gene expression for each cell type}

\item{labelednodes}{List of high scoring eigenvector centrality receptors to seed
the label propagation clustering}
}
\value{
Communities in LR network
}
\description{
Community detection ft. label propagation and louvain clustering on
the LR network
}
