####该代码基于pbmc结果用R版本SCENIC进行分析，参考https://cloud.tencent.com/developer/article/1692240

#SCENIC分析流程，官方介绍的主要分析有四步：
#GENIE3/GRNBoost：基于共表达情况鉴定每个TF的潜在靶点；
#RcisTarget：基于DNA-motif 分析选择潜在的直接结合靶点；
#AUCell：分析每个细胞的regulons活性；
#细胞聚类：基于regulons的活性鉴定稳定的细胞状态并对结果进行探索

#SCENIC在R中实现基于三个R包：
# GENIE3：推断基因共表达网络
# RcisTarget：用于分析转录因子结合motif
# AUCell：用于鉴定scRNA-seq数据中具有活性基因集（基因网络）的细胞

# rm(list=ls())  #删除当前R环境中的所有对象（变量、函数等）
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)




##分析准备
setwd("~/精神疾病/data/388/Seurat/data/SCENIC")
dir.create("ALL")
dir.create("ALL/int")
scRNA <- readRDS("~/精神疾病/data/388/Seurat/data/SCTed_new.rds")
# scRNA11 <- readRDS("~/精神疾病/data/388/Seurat/data/SCTed_new.rds")
# scRNA0 <- readRDS("/home/chenzh/brain/SCTed.rds")

##准备细胞meta信息
scRNA@meta.data$cell_type<-Idents(scRNA)
cellInfo <- data.frame(scRNA@meta.data)
cellInfo<- cellInfo[, c('group','cell_type')]#'seurat_clusters','cell_type'
colnames(cellInfo)=c('group','celltype')
saveRDS(cellInfo, "~/精神疾病/data/388/Seurat/data/SCENIC/ALL/cellInfo.rds")







##-----------准备表达矩阵------------##

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
# CON  ASD   BD  MDD PTSD   SZ 
# 400  400  400  400  400  400

#scRNA_B<- subset(scRNA, subset = cell_type=="B")
#exprMat <- as.matrix(scRNA_B[["RNA"]]$count[1:500,1:100])#取少部分基因和细胞，保证能运行过去，但结果不准确
#exprMat <- as.matrix(scRNA[["RNA"]]$count[1:500,1:100])#取少部分基因和细胞，保证能运行过去，但结果不准确




##设置分析环境
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
#nCores=8，代表开启8个线程计算，参数根据自己电脑CPU情况设置；
#mydbDIR <- "./SCENIC"，设置的是存放数据库的目录，需要填写自己存放数据库的目录；dbs = mydbs，在变量设置中，把hg38数据库文件赋值给了mydbs，用hg19数据库或小鼠的数据库需要相应调整。
#人org="hgnc", or 小鼠org="mgi", or 果蝇"dmel"







##转录调控网络推断
#install.packages("doRNG")
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
#根据表达数据推断潜在的转录因子靶标，使用 GENIE3 或 GRNBoost，GENIE3 非常耗时且计算量大（在 3-5k 单元的数据集上需要几个小时或几天的时间）
#GRNboost可在很短的时间内提供与 GENIE3 类似的结果，这使用的R，选择GENIC3
##nParts参数，是把表达矩阵分成n份分开计算，如nParts = 10



#这一步消耗的计算资源非常大，个人电脑需要几个小时的运行时间
runGenie3(exprMat_filtered_log, scenicOptions)

#以上代码运行后，int目录下有不少中间结果产生，简要解释一下：
#1.2_corrMat.Rds：基因之间的相关性矩阵
#1.3_GENIE3_weightMatrix_part_1.Rds等：GENIE3的中间结果
#1.4_GENIE3_linkList.Rds：GENIE3最终结果，是把“1.3_”开头的文件合并在一起




#推断共表达模块
runSCENIC_1_coexNetwork2modules(scenicOptions)
##推断转录调控网络（regulon），此步结果在output文件夹，之前的都在int文件夹。内存不够的电脑在这一步运行时会报错
runSCENIC_2_createRegulons(scenicOptions)#真正跑用这个
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量
#runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))#举例加参数

#SCENIC结果中regulon的名称有两种，一种是TF名称+extended+靶基因数目，
#另一种是TF名称+靶基因数目。第一种是转录因子与所有靶基因组成的基因调控网络，第二种是转录因子与高可信靶基因（即highConfAnnot=TRUE的基因）组成的基因调控网络。










##------------regulon活性评分与可视化------------##
#regulons计算AUC值并进行下游分析


library(doParallel)
exprMat_all <- as.matrix(scRNA[["RNA"]]@counts)
exprMat_all <- log2(exprMat_all+1)
scenicOptions@settings$nCores <- 1 ##重新设置核数，否则会报doMC 包的问题
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

#二进制转换及衍生分析（将评分转换为表示细胞是否在某个调控模块上活跃的0或1） 
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)

#SCENIC结果可视化
#runSCENIC_3_scoreCells()和runSCENIC_4_aucell_binarize()
#运行之后会生成一些可视化的heatmap图与tSNE图，但是他们既不容易与seurat分析的结果联系起来，又不容易调整图形参数和分析内容。
#我们可以调用SCENIC的分析结果，使用seurat和pheatmap进行可视化。
#Seurat可视化SCENIC结果：把SCENIC结果中最重要的regulonAUC矩阵导入Seurat，这样得到的可视化结果更容易与我们之前的分析联系起来


##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA, AUCmatrix)
# scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA, BINmatrix)
# scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

##利用Seurat可视化AUC
#FeaturePlot
scRNAauc <- RunPCA(scRNAauc, features = VariableFeatures(object = scRNAauc), npcs = 30)
scRNAauc <- RunUMAP(scRNAauc, dims = 1:30)

p1 = FeaturePlot(scRNAauc, features='POU2F1_479', label=T, reduction ='umap')
p2 = FeaturePlot(scRNAbin, features='POU2F1_479', label=T, reduction = 'umap')
p3 = DimPlot(scRNA, reduction = 'umap', group.by = "cell_type", label=T)
p1|p2|p3

#pheatmap可视化SCENIC结果 
library(pheatmap)
cellInfo <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/ALL/cellInfo.rds")
celltype = subset(cellInfo,select = 'group')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

#挑选部分感兴趣的regulons
my.regulons <- c('SREBF1_extended 2334g','ZNF274 extended 638g','BRF1 extended_1583g')##需要根据自己结果调整选择
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]


#使用regulon原始AUC值绘制热图
pheatmap(myAUCmatrix, 
         show_colnames=F, 
         annotation_col=celltype, 
         width = 6, height = 5)
#使用regulon二进制AUC值绘制热图
pheatmap(myBINmatrix, 
         show_colnames=F, 
         annotation_col=celltype,
         color = colorRampPalette(colors = c("white","black"))(100),
         width = 6, height = 5)
