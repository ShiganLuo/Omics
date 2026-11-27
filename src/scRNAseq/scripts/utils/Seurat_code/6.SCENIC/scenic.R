####该代码基于pbmc结果用R版本SCENIC进行分析，参考https://cloud.tencent.com/developer/article/1692240



library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)




setwd("~/精神疾病/data/388/Seurat/data/SCENIC")
dir.create("ALL")
dir.create("ALL/int")
scRNA <- readRDS("~/精神疾病/data/388/Seurat/data/SCTed_new.rds")

scRNA@meta.data$cell_type<-Idents(scRNA)
cellInfo <- data.frame(scRNA@meta.data)
cellInfo<- cellInfo[, c('group','cell_type')]#'seurat_clusters','cell_type'
colnames(cellInfo)=c('group','celltype')
saveRDS(cellInfo, "~/精神疾病/data/388/Seurat/data/SCENIC/ALL/cellInfo.rds")








scRNA1 <- subset(x = scRNA, group == "CON")  
subcell1 <- sample(colnames(scRNA1),400)
scRNA2 <- subset(x = scRNA, group == "ASD")  
subcell2 <- sample(colnames(scRNA2),400)
scRNA3 <- subset(x = scRNA, group == "BD")  
subcell3 <- sample(colnames(scRNA3),400)
scRNA4 <- subset(x = scRNA, group == "MDD")  
subcell4 <- sample(colnames(scRNA4),400)
scRNA5 <- subset(x = scRNA, group == "PTSD")  
subcell5 <- sample(colnames(scRNA5),400)
scRNA6 <- subset(x = scRNA, group == "SZ")  
subcell6 <- sample(colnames(scRNA6),400)

scRNAsub <- scRNA[,c(subcell1,subcell2,subcell3,subcell4,subcell5,subcell6)]
saveRDS(scRNAsub, "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/scRNAsub.rds")
exprMat <- as.matrix(scRNAsub@assays[["RNA"]]@counts)

table(scRNAsub@meta.data$group)





mydbDIR <- "~/精神疾病/data/388/Seurat/data/SCENIC/ALL"
mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")#注意下载old下面支持R
names(mydbs) <- c("500bp", "10kb")
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "123")
saveRDS(scenicOptions, "~/精神疾病/data/388/Seurat/data/SCENIC/ALL/scenicOptions.rds")







genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)



runGenie3(exprMat_filtered_log, scenicOptions)





runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)#真正跑用这个













library(doParallel)
exprMat_all <- as.matrix(scRNA[["RNA"]]@counts)
exprMat_all <- log2(exprMat_all+1)
scenicOptions@settings$nCores <- 1
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)



AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA, AUCmatrix)
saveRDS(scRNAauc,'scRNAauc.rds')

BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA, BINmatrix)
saveRDS(scRNAbin, 'scRNAbin.rds')

scRNAauc <- RunPCA(scRNAauc, features = VariableFeatures(object = scRNAauc), npcs = 30)
scRNAauc <- RunUMAP(scRNAauc, dims = 1:30)

p1 = FeaturePlot(scRNAauc, features='POU2F1_479', label=T, reduction ='umap')
p2 = FeaturePlot(scRNAbin, features='POU2F1_479', label=T, reduction = 'umap')
p3 = DimPlot(scRNA, reduction = 'umap', group.by = "cell_type", label=T)
p1|p2|p3

library(pheatmap)
cellInfo <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/ALL/cellInfo.rds")
celltype = subset(cellInfo,select = 'group')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

my.regulons <- c('SREBF1_extended 2334g','ZNF274 extended 638g','BRF1 extended_1583g')##需要根据自己结果调整选择
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]


pheatmap(myAUCmatrix, 
         show_colnames=F, 
         annotation_col=celltype, 
         width = 6, height = 5)
pheatmap(myBINmatrix, 
         show_colnames=F, 
         annotation_col=celltype,
         color = colorRampPalette(colors = c("white","black"))(100),
         width = 6, height = 5)
