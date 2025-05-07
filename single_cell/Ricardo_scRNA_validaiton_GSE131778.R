library(tidyverse)
library(dplyr)
library(readr)
library(AnnotationDbi)
library(textshape)
library(tibble)
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
library(pathview)
library(VennDiagram)
library(Seurat)
library(cowplot)
library(clusterProfiler.dplyr)
library(clusterProfiler)
require(ComplexHeatmap)
require(circlize)
library(viridis)
##Load DEG's table
######################################################################################
setwd("C:/Users/Ricardo/Desktop/Single_cell_R")

coronary_plaque <- readRDS("coronary_plaque_cell_annotated.rds")
#GENES_tree <- read_csv("GSE131778/Ricardo_scRNA_validation/genes_gtex_tree_on_depth_3.csv")
#GENES_tree_new <- read_csv("GSE131778/Ricardo_scRNA_validation/six_genes_from_gtex_tree.csv")
# set1 <- GENES_tree[c(1:20),]
# set2 <- GENES_tree[c(21:40),]
# set3 <- GENES_tree[c(41:52),]


genes <- read.csv('genes_07.csv')
genes <- genes$GENE_NAMES

cluster.averages <- AverageExpression(coronary_plaque, use.scale = FALSE) #average expression of raw data by tissue
EC_Expressed_tab <- as.data.frame(cluster.averages) #creates a dataframe from average data
EC_Expressed_tab <- select(EC_Expressed_tab, RNA.Endothelial.Cell.1, RNA.Endothelial.Cell.2) #selects only columns from endothelial tissue
EC_genes <- row.names(EC_Expressed_tab)
EC_Expressed_tab <- EC_Expressed_tab[EC_genes %in% genes,] #dataframe only with genes from correlation analysis
EC_Expressed_tab <- EC_Expressed_tab[EC_Expressed_tab$RNA.Endothelial.Cell.1 > 0.0125 & EC_Expressed_tab$RNA.Endothelial.Cell.2 > 0.0125,] #filtering by expression threshold 
EC_Expressed_tab <- EC_Expressed_tab[order(-EC_Expressed_tab$RNA.Endothelial.Cell.1),] #sorting by desc expression
endothelial_genes <- row.names(EC_Expressed_tab)

write.csv(endothelial_genes,"list_of_endothelial_genes.csv")

##cell vizualization
setwd("umap")
tiff('SDCBP_umap.tif', width = 7, height = 7, units = 'in', res=300)
FeaturePlot(coronary_plaque, features = "SDCBP" ,reduction ="umap",  cols = c("whitesmoke", "indianred"), 
            pt.size = 1, label = TRUE, label.size = 3)
#VlnPlot(coronary_plaque, features = 'SDCBP', pt.size=1)
dev.off()

tiff('PAPPA_umap.tif', width = 7, height = 7, units = 'in', res=300)
FeaturePlot(coronary_plaque, features = "PAPPA" ,reduction ="umap",  cols = c("whitesmoke", "indianred"), 
            pt.size = 1, label = TRUE, label.size = 3)
#VlnPlot(coronary_plaque, features = 'SDCBP', pt.size=1)
dev.off()

###annotating cell type

tiff("tUMAP_ALL_cell_type.tif",width = 7, height = 7, units = 'in', res=300)
DimPlot(coronary_plaque, reduction = "umap", label = TRUE, pt.size = 1)+ NoLegend()
dev.off()
