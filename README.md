# REMI (Regularized Microenvironment Interactome)

## Tutorials
[Bulk Flow-Sorted RNA-Seq in R](http://htmlpreview.github.io/?https://github.com/ayu1/remi/blob/master/vignettes/REMI_Tutorial.html)

- Single cell RNA-Seq data in R
- Python adaption for REMI

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

For usage tutorials, check the `vignettes` directory of the repo. 

## Contact
Email: ayu1@stanford.edu
