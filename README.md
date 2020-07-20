# REMI (Regularized Microenvironment Interactome)

## Quick Start
R Tutorial  
Single cell R tutorial   
Python Tutorial  

## Description
The REMI R package is designed to predict ligand receptor interactions within a microenvironment given RNA-sequencing data. It regularizes communities of ligand and receptors that have a high probability of being transcriptionally linked. 

REMI has been implemented in R (>3.6)

## Installation
```
install.packages("devtools")
library(devtools)
install_github("ayu1/remi")
```

## Usage

REMI takes in normalized and scaled RNA-sequencing data as an input, where the columns are samples and rows are genes. The column names are labeled as sample_celltype (i.e. S01_Bcell). The package has a built-in option to filter low-expressed genes, but it can also take it any pre-filtered scaled datasets (filter=F). 
```
interactome <- remi(data, filter=T)
```

## Contact
Email: ayu1@stanford.edu
