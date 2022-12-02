# REMI (Regularized Microenvironment Interactome)

![plot](https://github.com/plevritis-lab/remi/blob/master/extra/remi_figure.png?raw=true)

## Description
The REMI package is written in R and is designed to predict ecosystem-wide ligand receptor interactions within a microenvironment given RNA-sequencing data. REMI moves beyond pairwise interactions and accounts for the effect of 2+ cell types on system-level interactions. Specifically, it creates communities of multicellular genes and identifies which interactions have a high probability of being transcriptionally linked. 

Manuscript: https://www.science.org/doi/10.1126/sciadv.abi4757

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

To use a different ligand-receptor database (i.e. mouse), set the variable lr.adtabase as the new matrix/data frame. 
Column names must be PairName, Ligand, and Receptor where PairName separates the Ligand and Recpetor by "_".
Gene names in uploaded table must match gene names in input file. 

```
new.lr.table <- read_csv("newlrtable.csv")

new.lr.table

   Pair.Name    Ligand Receptor
   <chr>        <chr>  <chr>   
 1 A2M_LRP1     A2M    LRP1    
 2 AANAT_MTNR1A AANAT  MTNR1A  
 3 AANAT_MTNR1B AANAT  MTNR1B 

remi.res <- remi(obj, lr.database = new.lr.table)
```

## Contact
Email: Alice Yu (ayu1@alumni.stanford.edu) or Sylvia Plevritis (sylvia.plevritis@stanford.edu)
