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
	'BiocManager'
  )
)

BiocManager::install(c('org.Hs.eg.db','msigdbr', 'clusterProfiler'))