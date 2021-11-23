# REMI (Regularized Microenvironment Interactome)


## Description
The REMI package is written in R and is designed to predict ligand receptor interactions within a microenvironment given RNA-sequencing data. It regularizes communities of ligand and receptors that have a high probability of being transcriptionally linked. 

Manuscript: https://www.biorxiv.org/content/10.1101/2021.05.02.440071v1

## Installation

Install R (>= 3.6) 
```
install.packages("devtools")
library(devtools)
install_github("plevritis-lab/remi")
```

## Usage

REMI takes in normalized bulk or single-cell RNA-sequencing data as an input, where the columns are samples and rows are genes. The column names are labeled as sample_celltype (i.e. S01_Bcell). The package has a built-in option to filter low-expressed genes, but it can also take it any pre-filtered scaled datasets (filter=F). For single-cell RNA-sequencing data, REMI can be run directly from the Seurat object.

### Tutorials (also available in  `vignettes` directory of the repo):

[Bulk RNA-Seq in R](http://htmlpreview.github.io/?https://github.com/plevritis-lab/remi/blob/master/vignettes/REMI_Tutorial.html)

[Single-cell RNA-Seq data in R](http://htmlpreview.github.io/?https://github.com/ayu1/remi/blob/master/vignettes/singleCell_REMITutorial.html)


## Contact
Email: Alice Yu (ayu1@stanford.edu) or Sylvia Plevritis (sylvia.plevritis@stanford.edu)
