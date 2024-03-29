% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_functions.R
\name{setupData}
\alias{setupData}
\title{Clean the data}
\usage{
setupData(dat, cellmarkers, filter = T, var = 3)
}
\arguments{
\item{dat}{Data matrix of normalized log2 transformed TPM counts}

\item{cellmarkers}{list of cell markers as labeled in the columns of the dataset.}

\item{var}{Number of standard deviations user would like to cut the data by}
}
\value{
A list of gene expression per cell type
}
\description{
This function loads a file as a matrix,
removes genes with low expression,
and creates cell type-specific gene expression matrices.
}
