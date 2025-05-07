library(tidyverse)
library(dplyr)
library(readr)
library(AnnotationDbi)
library(textshape)
library(tibble)
library("org.Hs.eg.db")
library(fgsea)
library(clusterProfiler)
library(ggplot2)
library(pathview)
##Load DEG's table
STAT_vs_OSS_DEGs <- read_csv("tabelas/OSS_STAT_DEGs.csv")

#1. CREATE RANKED GENE LIST
#Add ID code
entrz_name_STAT_OSS <- AnnotationDbi::select(org.Hs.eg.db,
                                             key=STAT_vs_OSS_DEGs$SYMBOL,  
                                             columns="ENTREZID",
                                             keytype="SYMBOL")
STAT_vs_OSS_ID_entrz <- inner_join(STAT_vs_OSS_DEGs, entrz_name_STAT_OSS, by="SYMBOL")
#Create a vector with LogFC values
nw_t_STAT_OSS<-STAT_vs_OSS_ID_entrz
STAT_OSS_ID_entrz_logfc =nw_t_STAT_OSS$logFC
#Add Correspondent IDs info
names(STAT_OSS_ID_entrz_logfc)<-as.vector(nw_t_STAT_OSS$ENTREZID)
#put it in decreasing order
STAT_OSS_gene_list = sort(STAT_OSS_ID_entrz_logfc, decreasing = TRUE)
#Verify
STAT_OSS_gene_list[1:10]

##2.CLUSTERPROFILER ENRICHMENT

#2.1 KEGG
#define what organism you're working with (HUMAN, RAT, MOUSE, RABBIT...)
kegg_organism = "hsa"
#Enrichment KEGG
#Use the ranked gene list previously generated
STAT_OSS_kegg <- gseKEGG(geneList     = STAT_OSS_gene_list,
                         organism     = kegg_organism,
                         nPerm        = 1000,
                         minGSSize    = 3,
                         maxGSSize    = 800,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "none", ## You may add iff necessary
                         keyType       = "ncbi-geneid")
#Use the table give to generate and save ridgeplot (you can control how many categories you want)
tiff("Enrichment/STAT_OSS_KEGG_ridgeplot_Enriched_Pathways.tif",width = 12, height = 6, units = 'in', res=150)
ridgeplot(STAT_OSS_kegg, showCategory = 20, fill = "pvalue") + labs(x = "enrichment distribution")
dev.off()

#For the next plots you need to use a data frame
#rename entrz codes
STAT_OSS_kegg <- setReadable(STAT_OSS_kegg, 'org.Hs.eg.db', 'ENTREZID')
#create data frame file
STAT_OSS_kegg_tab <- as.data.frame(STAT_OSS_kegg)
write.csv(STAT_OSS_kegg_tab, "Enrichment/STAT_vs_OSS_kegg_tab.csv")

#generate Dotplot(use split if you want tto have "supressed" and "activated" pathways)
tiff("Enrichment/STAT_OSS_KEGG_dotplot_Enriched_Pathways.tif",width = 10, height = 6, units = 'in', res=150)
dotplot(STAT_OSS_kegg, showCategory = 10, title = "OSS vs STAT KEGG Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()

#generate the emmaplot
tiff("Enrichment/STAT_OSS_KEGG_emapplot_Enriched_Pathways.tif",width = 12, height = 6, units = 'in', res=150)
emapplot(STAT_OSS_kegg, showCategory = 42, color = "pvalue")
dev.off()

#2.2 GENE ONTOLOGY ENRICHEMENT
#use ranked gene list for gseGO
STAT_OSS_go <- gseGO(geneList     = STAT_OSS_gene_list,
                     ont = "BP",
                     OrgDb = org.Hs.eg.db,
                     nPerm        = 500,
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "none", #no p value adjusted
                     keyType       = "ENTREZID")

#Use the table give to generate and save ridgeplot (you can control how many categories you want)
tiff("Enrichment/STAT_OSS_GO_ridgeplot_Enriched_Pathways.tif",width = 12, height = 6, units = 'in', res=150)
ridgeplot(STAT_OSS_go, showCategory = 20, fill = "pvalue") + labs(x = "enrichment distribution")
dev.off()

#For the next plots you need to use a data frame
#rename entrz codes
STAT_OSS_go <- setReadable(STAT_OSS_go, OrgDb = org.Hs.eg.db) 
#create data frame file
STAT_OSS_GO_tab <- as.data.frame(STAT_OSS_go)
write.csv(STAT_OSS_GO_tab, "Enrichment/STAT_vs_OSS_GO_tab.csv")

#generate Dotplot(use split if you want tto have "supressed" and "activated" pathways)
tiff("Enrichment/STAT_OSS_GO_dotplot_Enriched_Pathways.tif",width = 10, height = 6, units = 'in', res=150)
dotplot(STAT_OSS_go, showCategory = 10, title = "OSS vs STAT GO Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()
#generate the emmaplot
tiff("Enrichment/STAT_OSS_GO_emapplot_Enriched_Pathways.tif",width = 12, height = 6, units = 'in', res=150)
emapplot(STAT_OSS_go, showCategory = 52, color = "pvalue")
dev.off()