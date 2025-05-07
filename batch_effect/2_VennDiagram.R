library(tidyverse)
library(dplyr)
library(readr)
library(gplots)
library(textshape)
library(tibble)
library(gameofthrones)
library(VennDiagram)  
#####################################################

# tables with all genes comparison
batch_correct_OSS_LSS_DEGs <-  read_csv("combat_OSS_LSS_DEGs.csv")
iguara_OSS_LSS_DEGs <- read_csv("iguara_OSS_LSS_DEGs.csv")
sarah_OSS_LSS_DEGs <- read_csv("sarah_OSS_LSS_DEGs.csv")


##venn diagram
List_overlaped_samples <- list( batch_correct_OSS_LSS_DEGs,iguara_OSS_LSS_DEGs, sarah_OSS_LSS_DEGs)
List_overlaped_samples = lapply(List_overlaped_samples, function(i) (i)$SYMBOL)


venn.diagram(List_overlaped_samples, filename = 'venn_sarah_iguara_combat.png', 
             category.names = c("ComBat (422 DEGs)", "Iguaracy (370 DEGs)",
                                "Sarah (762 DEGs)"),
             output = NULL,
             imagetype="png",
             height = 1100 ,
             width = 1100 ,
             resolution = 200,
             compression = "lzw",
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("#f7ca64", "#46bac2","#7e62a3"),
             cex = 1.2,
             fontfamily = "sans",
             cat.cex = 0.9,
             cat.default.pos = "outer",
             cat.pos = c(-10, 10, -190),
             cat.fontfamily = "sans",
             
             cat.col = c("#f7ca64", "#46bac2","#7e62a3"))

names_List_genes <- c("Batch (404 DEGs)", "Iguara (338 DEGs)",
                      "Sarah (649 DEGs)")

names(List_overlaped_samples) <- names_List_genes
intersect_names_List_genes <- venn(List_overlaped_samples)
intersect_names_List_to_subset <- attr(intersect_names_List_genes, "intersection") ##pegar todas as intersec?oes

teste <- as.data.frame(intersect_names_List_to_subset$`Batch (409 DEGs):Sarah (649 DEGs)`) ##selecionar qual intersec??o
colnames(teste)<- "SYMBOL"
teste1 <- inner_join(sarah_LSS_OSS_DEGs, teste, by = "SYMBOL")
write.csv(teste1, "batch_sarah_inter.csv")

