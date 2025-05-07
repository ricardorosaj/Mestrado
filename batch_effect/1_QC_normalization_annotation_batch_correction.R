library(oligo)
library(affy)
library(arrayQualityMetrics)
library(ggplot2)
library(Biobase)
library(pheatmap)
library(clariomdhumantranscriptcluster.db)
library(dplyr)
library(viridis)
library(limma)
library(cartography)
library(paletteer)
library(rcartocolor)
library(sva)
library(BatchQC)
library(rhdf5)
library(tools)
library(dendextend)
library(preprocessCore)
library(vioplot)
## Set working directory
setwd("samples")

#Import files names
list_out = list.celfiles(getwd(),full.names=TRUE, pattern = "CEL")

#Read cel files
files.cel<- read.celfiles(list_out)

#Add Pheno data informations

##Conditions
pData(files.cel)[,1]<- c("LSS", "LSS","LSS","LSS",
                        "LSS", "LSS","LSS","LSS",
                        "OSS","OSS","OSS",
                        "OSS","OSS","OSS","OSS")
# pData(files.cel)[,1]<- c("LSS1", "LSS1","LSS1",'LSS1',
#                          "LSS2", "LSS2","LSS2","LSS2")
# pData(files.cel)[,1]<- c("OSS1", "OSS1","OSS1",
#                          "OSS2", "OSS2","OSS2","OSS2")
# pData(files.cel)[,1]<- c("LSS", "LSS","LSS","LSS",
#                          "OSS","OSS","OSS","OSS")
##Sample names
pData(files.cel)[,2]<- c("LSS1","LSS2","LSS3","LSS4",
                         "LSS5","LSS6","LSS7","LSS8",
                         "OSS1","OSS2","OSS4",
                         "OSS5","OSS6","OSS7","OSS8")
# pData(files.cel)[,2]<- c("LSS1","LSS2","LSS3","LSS4",
#                          "LSS5","LSS6","LSS7","LSS8")
# pData(files.cel)[,2]<- c("OSS1","OSS2","OSS4",
#                          "OSS5","OSS6","OSS7","OSS8")
# pData(files.cel)[,2]<- c("LSS1","LSS2","LSS3","LSS4",
#                          "OSS1","OSS2","OSS3","OSS4")
##Batch of each sample
pData(files.cel)[,3]<- c("1","1","1","1",
                         "2","2","2","2",
                         "1","1","1",
                         "2","2","2","2")
# pData(files.cel)[,3]<- c("1","1","1",
#                          "2","2","2","2")
# pData(files.cel)[,3]<- c("1","1","1",
#                          "2","2","2","2")

## Change column names
colnames(pData(files.cel))[1]="Conditions"
colnames(pData(files.cel))[2]="Samples_name"
colnames(pData(files.cel))[3]="BatchID"

pheno <- pData(files.cel)
pheno

# #Quality control check
# arrayQualityMetrics(expressionset = files.cel,
#                    force = TRUE, do.logtransform = TRUE,
#                    intgroup = c("Samples_name"))

#Normalizations

#Log2 transformations
exp_raw <- log2(Biobase::exprs(files.cel))

#Coordenates for PCA prior normalization (Group cluster)
PCA_raw <- prcomp(t(exp_raw), scale. = T, center=T)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG_non_normalized <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                                    Conditions = pData(files.cel)$Conditions,
                                    Samples_name = pData(files.cel)$Samples_name,
                                    BatchID = pData(files.cel)$BatchID)


ggplot(dataGG_non_normalized, aes(PC1, PC2)) +
  geom_point(aes(colour = Conditions), size=6) +
  ggtitle("PCA plot of the log-transformed raw expression data")+
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) 

#Boxplot of expression deviation
oligo::boxplot(files.cel)

#RMA normalization (Background correction, Quantile normalization)
norm.file.cel <- oligo::rma(files.cel, normalize = T)

oligo::boxplot(norm.file.cel)

#Quality control check
# arrayQualityMetrics(expressionset = norm.file.cel,
#                    force = TRUE, do.logtransform = TRUE,
#                    intgroup = c("Samples_name"))

## Filter lowly expressed probes
norm.expr.median<- rowMedians(Biobase::exprs(norm.file.cel))

#tiff("imagens/intensity_hist.tif",width = 13, height =13, units = 'in', res=250)
hist_res<- hist(norm.expr.median, 100, freq = FALSE, 
                main = "Histogram of the median intensities - Both Batches", 
                border = "antiquewhite4",
                xlab = "Median intensities")
manual_threshold <- 3.0
abline(v = manual_threshold, col = "red", lwd = 2)
#dev.off()


#Define the number of replicated in treatment to define the cutoff
no_of_samples <-table(paste0(pData(norm.file.cel)$Conditions))

pData(norm.file.cel)

no_of_samples 

samples_cutoff <- min(no_of_samples)

#Function to define the number of probes that do not passes minimum median fluorescence intensity in the data
idx_man_threshold <- apply(Biobase::exprs(norm.file.cel), 1,
                           function(x){
                             sum(x > manual_threshold) >= samples_cutoff})
table(idx_man_threshold)

#sUBSET THE DATA WITH ONLY PROBES THAT PASSES THE FILTERING STEP
array.norm.filtered<- subset(norm.file.cel, idx_man_threshold)

oligo::boxplot(array.norm.filtered)

####
#PCA after normalization
#Log2 transformations
exp_norm.filtered <- log2(Biobase::exprs(array.norm.filtered))

EXP.norm.nolog <- (Biobase::exprs(array.norm.filtered))

#Coordinates for PCA after normalization and probe filtering (Group cluster)
PCA_norm.filtered <- prcomp(t(exp_norm.filtered), scale. = T, center=T)

percentVar.filtered <- round(100*PCA_norm.filtered$sdev^2/sum(PCA_norm.filtered$sdev^2),1)

sd_ratio.filtered <- sqrt(percentVar.filtered[2] / percentVar.filtered[1])

dataGG_non_normalized.filtered <- data.frame(PC1 = PCA_norm.filtered$x[,1], PC2 = PCA_norm.filtered$x[,2],
                                    Conditions = pData(array.norm.filtered)$Conditions,
                                    Samples_name = pData(array.norm.filtered)$Samples_name,
                                    BatchID = pData(array.norm.filtered)$BatchID)
#choose color pallet 
my.colors1=carto_pal(7, "TealGrn")
my.colors2=carto_pal(7, "BluGrn")

#plot the pca
#tiff("imagens/PCA_conditions.tif",width = 15, height =15, units = 'in', res=250)
ggplot(dataGG_non_normalized.filtered, aes(PC1, PC2)) +
  geom_point(aes(colour = Conditions), size=7) +
  scale_color_manual(values = c("#96d2a4","#257d98","#4cc8a3"))+
  ggtitle("Principal components analysis (PCA) of filtered expression data")+
  xlab(paste0("PC1, VarExp: ", percentVar.filtered[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar.filtered[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio.filtered)+
  theme_bw()
#dev.off()
#Create distance matrix to visualize data clustering
sub_dists <- as.matrix(dist(t(exp_norm.filtered), method = "manhattan"))

rownames(sub_dists) <- pData(files.cel)$Samples_name
colnames(sub_dists)<- pData(files.cel)$Samples_name

#Change diagonal data in matrix to zero
diag(sub_dists) <- 0

#Heatmap of expressions distance after normalization
#tiff("imagens/heatmap.tif",width = 18, height =18, units = 'in', res=250)
pheatmap(sub_dists,
         color = mako(10),
         angle_col=45,
         legend_breaks = c(min(sub_dists, na.rm = TRUE), 
                          max(sub_dists, na.rm = TRUE)), 
         legend_labels = (c("small distance", "large distance")),
         treeheight_row = 30,
         treeheight_col = 30,
         cellwidth= 50,
         cellheight= 50,
         border_color="white",
         main = "Gene Expression Profile Distance",
         number_color= "blue4",
         fontsize=16)
#dev.off()


#Annotate probe to gene names and symbols
annot.clariom.ht <- AnnotationDbi::select(clariomdhumantranscriptcluster.db,
                                       keys = (featureNames(array.norm.filtered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

#Remove NA gene names that downloaded from the database
annotation_clariom <- subset(annot.clariom.ht, !is.na(SYMBOL))
annotation_clariom<-annotation_clariom[!duplicated(annotation_clariom$SYMBOL),]

#Steps form Removing multiple mapping from probes

#1- Group data using the probeID
anno_grouped <- annotation_clariom %>%
  dplyr::group_by(PROBEID) 

#2- Check which probes have multiple gene mapping
anno_summarized <- anno_grouped %>%
  dplyr::summarize(no_of_matches = n_distinct(SYMBOL)) 

#3- Remove those who have more than 1 mapping per gene
anno_filtered <- anno_summarized %>%
  filter(no_of_matches > 1) 

#4- Define which IDs are going to be excluded
ids_to_exlude <- (featureNames(array.norm.filtered) %in% anno_filtered$PROBEID)

#5- Generate the final dataset without multiple mapping probes
array.norm.filtered.nomapping <- array.norm.filtered %>%
  subset(!ids_to_exlude)

#Add the row names as a featureData for the expressionset objetc
fData(array.norm.filtered.nomapping)$PROBEID <- rownames(fData(array.norm.filtered.nomapping))

#Add the columns GENE and SYMBOL to the featureData for only the probes that have been filtered
fData(array.norm.filtered.nomapping)<- dplyr::left_join(fData(array.norm.filtered.nomapping), 
                                                        annotation_clariom, "PROBEID")

#restore the probeIDs names to the dataset
rownames(fData(array.norm.filtered.nomapping)) <- fData(array.norm.filtered.nomapping)$PROBEID 

array.norm.filtered.nomapping<- array.norm.filtered.nomapping %>%
  subset(!is.na(fData(array.norm.filtered.nomapping)$SYMBOL))

##########################
#BATCH CORRECTION#

expr.matrix <- exprs(array.norm.filtered.nomapping)
#write.exprs(array.norm.filtered.nomapping, file='exprs_no_batch.txt', sep='\t')

batch <- pheno$BatchID
cond <- pheno$Conditions
# batchQC(expr.matrix, batchid, condition = cond)

#Dendrogram plot before batch correction
cc = cor(exp)
dend <- as.dendrogram(hclust(as.dist(1-cc)))
useries = unique(batchid)
series_match = batchid
colos <- colorspace::rainbow_hcl(length(useries), c = 160, l  = 50)
names(colos) = useries
series_color <- colos[series_match]
clu = cutree(dend, h=0.25)
labels_colors(dend) <- series_color[order.dendrogram(dend)]
dend <- color_branches(dend, h = 0.25)
par(mar = c(4,1,1,12))
plot(dend, horiz = TRUE, xlab='Distance between samples (1-cc)', main='Before Batch Correction')
colored_bars(cbind(clu, series_color), dend, rowLabels = c("Cluster", "Series"), horiz = TRUE)
legend("topleft", legend = useries, fill = colos, bg="white", cex=0.6)

# Executing SVA correction
modmatrix = model.matrix(~as.factor(Conditions), data=pheno)
sva.object <- batchQC_sva(expr.matrix, mod=modmatrix)
correctedExpression <- batchQC_svregress_adjusted(expr.matrix,modmatrix,sva.object)
# write.table(correctedExpression, file='exprs_sva.txt', sep='\t')

#Dendrogram plot after batch correction
cc = cor(correctedExpression)
dend <- as.dendrogram(hclust(as.dist(1-cc)))
useries = unique(batchid)
series_match = useries[match(batchid, useries)]
colos <- colorspace::rainbow_hcl(length(useries), c = 160, l  = 50)
names(colos) = useries
series_color <- colos[series_match]
clu = cutree(dend, h=0.25)
labels_colors(dend) <- series_color[order.dendrogram(dend)]
dend <- color_branches(dend, h = 0.25)
par(mar = c(4,1,1,12))
plot(dend, horiz = TRUE, xlab='Distance between samples (1-cc)', xlim=c(0.1,0), main='After Batch Correction - SVA')
colored_bars(cbind(clu, series_color), dend, rowLabels = c("Cluster", "Series"), horiz = TRUE)
legend("topleft", legend = useries, fill = colos, bg="white", cex=0.6)

#Correlation of samples 
cc = cor(exp)
dend <- as.dendrogram(hclust(as.dist(1-cc)))
clu = cutree(dend, h=0.25)
largest_cluster = names(rev(sort(table(clu))))[1]
ww = which(clu == largest_cluster)
plot(density(cor(exp[,ww])), lwd=3, main="correlation of leftover samples", ylim=c(0,30))
lines(density(cor(correctedExpression)), lwd=3, main="correlation of leftover samples", col="red")
legend("topleft", legend=c("uncorrected","corrected - SVA"), lty=1, lwd=3, col=c("black","red"))

##########################

########Find DEGs
#1- Define the groups to compare
groups = pData(array.norm.filtered.nomapping)$Conditions
f = factor(groups,levels=c("LSS","OSS"))
# f = factor(groups,levels=c("OSS1","OSS2"))
# f = factor(groups,levels=c("LSS1","LSS2"))

#2- Create Design Matrix
design = model.matrix(~ 0 + f)
colnames(design) = c("LSS","OSS")
# colnames(design) = c("OSS1","OSS2")
# colnames(design) = c("LSS1","LSS2")

#3- Fit the contrasts to Lineal model
# data.fit = lmFit(array.norm.filtered.nomapping,design)
data.fit = lmFit(correctedExpression,design)
data.fit$coefficients[1:10,]

#4- Define the contrasts (what to compare)
contrast.matrix = makeContrasts(OSS-LSS,
                                LSS-OSS,
                                levels=design)

# contrast.matrix = makeContrasts(OSS1-OSS2,
#                                 OSS2-OSS1,
#                                 levels=design)

# contrast.matrix = makeContrasts(LSS1-LSS2,
#                                 LSS2-LSS1,
#                                 levels=design)

#5- Fit the contrasts
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

#save(data, file = "Data_matrix_exp_Control_Samples1_combat.RData")

#6-DEGs summary
#DEresults = decideTests(data.fit.eb,coef=1,method="separate",adj.p.Val=0.05,lfc=1.3,adjust.method="none")
#summary(DEresults)

#save all comparison
OSS_vs_LSS <- topTreat(data.fit.eb, coef=1, n=Inf)
#create columns for venn diagram (batch correction for some reason removes then)
pos_of_probeid <- match(rownames(OSS_vs_LSS), featureData(array.norm.filtered.nomapping)$PROBEID)
OSS_vs_LSS$PROBEID <- featureData(array.norm.filtered.nomapping)$PROBEID[pos_of_probeid]
OSS_vs_LSS$SYMBOL <- featureData(array.norm.filtered.nomapping)$SYMBOL[pos_of_probeid]
OSS_vs_LSS$GENENAME <- featureData(array.norm.filtered.nomapping)$GENENAME[pos_of_probeid]
OSS_vs_LSS <- OSS_vs_LSS[,c(7,8,9,1,2,3,4,5,6)]
write.csv(OSS_vs_LSS, "_tabelas/combat_OSS_vs_LSS_allgenes.csv")

LSS_vs_OSS <- topTreat(data.fit.eb, coef=2, n=Inf)
# create columns for venn diagram (batch correction for some reason removes then)
pos_of_probeid <- match(rownames(LSS_vs_OSS), featureData(array.norm.filtered.nomapping)$PROBEID)
LSS_vs_OSS$PROBEID <- featureData(array.norm.filtered.nomapping)$PROBEID[pos_of_probeid]
LSS_vs_OSS$SYMBOL <- featureData(array.norm.filtered.nomapping)$SYMBOL[pos_of_probeid]
LSS_vs_OSS$GENENAME <- featureData(array.norm.filtered.nomapping)$GENENAME[pos_of_probeid]
LSS_vs_OSS <- LSS_vs_OSS[,c(7,8,9,1,2,3,4,5,6)]
write.csv(LSS_vs_OSS, "_tabelas/combat_LSS_vs_OSS_allgenes.csv")

# OSS1_vs_OSS2 <- topTreat(data.fit.eb, coef=1, n=Inf)
# write.csv(OSS1_vs_OSS2, "_tabelas/OSS1_vs_OSS2_allgenes.csv")
# OSS2_vs_OSS1 <- topTreat(data.fit.eb, coef=1, n=Inf)
# write.csv(OSS2_vs_OSS1, "_tabelas/OSS2_vs_OSS1_allgenes.csv")

# LSS1_vs_LSS2 <- topTreat(data.fit.eb, coef=1, n=Inf)
# write.csv(LSS1_vs_LSS2, "_tabelas/LSS1_vs_LSS2_allgenes.csv")
# LSS2_vs_LSS1 <- topTreat(data.fit.eb, coef=1, n=Inf)
# write.csv(LSS2_vs_LSS1, "_tabelas/LSS2_vs_LSS1_allgenes.csv")

#save threshold tables
 OSS_vs_LSS %>%
   filter(abs(logFC)>= 1.3, adj.P.Val<0.05) %>%
   write.csv("_tabelas/combat_OSS_LSS_DEGs.csv")

 LSS_vs_OSS %>%
   filter(abs(logFC)>= 1.3,adj.P.Val<0.05) %>%
   write.csv("_tabelas/combat_LSS_OSS_DEGs.csv")

 LSS_vs_OSS %>%
   filter(logFC>= 1.3,adj.P.Val<0.05) %>%
   write.csv("_tabelas/combat_LSS_OSS_DEGsUP.csv")

# OSS1_vs_OSS2 %>%
#   filter(abs(logFC)>= 1.3, adj.P.Val<0.05) %>%
#   write.csv("_tabelas/OSS1_OSS2_DEGs.csv")
# 
# OSS2_vs_OSS1 %>%
#   filter(abs(logFC)>= 1.3, adj.P.Val<0.05) %>%
#   write.csv("_tabelas/OSS2_OSS1_DEGs.csv")

# LSS1_vs_LSS2 %>%
#   filter(abs(logFC)>= 1.3, adj.P.Val<0.05) %>%
#   write.csv("_tabelas/LSS1_LSS2_DEGs.csv")
# 
# LSS2_vs_LSS1 %>%
#   filter(abs(logFC)>= 1.3, adj.P.Val<0.05) %>%
#   write.csv("_tabelas/LSS2_LSS1_DEGs.csv")

setwd('..')
