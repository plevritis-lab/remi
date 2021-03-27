install.packages(
  c(
	'Seurat',
	'qgraph', 
	'data.table', 
	'igraph', 
	'glasso', 
	'tidyverse',
	'OneR', 
	'Matrix', 
	'circlize'
	'preprocessCore'
	'viridis'
	'networkD3'
  )
)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("msigdbr")