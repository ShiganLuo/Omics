################################# 加载包 ################################
# rm(list = ls())
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) 
library(data.table) 
library(Seurat)
getwd()
setwd("~/精神疾病/data/388/Seurat/data/WGCNA")
### 启用WGCNA多核计算
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 
getwd()



################################# 0.输入数据准备 ################################
seurat.obj <- readRDS("/home/maxr/精神疾病/data/388/Seurat/data/SCTed_new.rds")


# 查看表达矩阵中的细胞数量
n_cells_counts <- ncol(seurat.obj@assays$RNA@counts)

# 查看元数据中的细胞数量
n_cells_meta <- nrow(seurat.obj@meta.data)







# 查看当前分组
Idents(seurat.obj) <- "group"
table(Idents(seurat.obj))
#  CON  ASD   BD  MDD PTSD   SZ 
# 2719 2378 2868 1551 2513 2838


#按照样本取平均，获得表达矩阵
# 对 Seurat 对象进行对数归一化处理
seurat.obj = NormalizeData(object = seurat.obj, assay = "RNA")
expr <- AverageExpression(seurat.obj, group.by = "group", assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
dim(expr)
# [1] 24664     6

# 更改矩阵行名为样本＋疾病，为表型数据做准备
seurat.obj@meta.data$sample = rownames(seurat.obj@meta.data)

expr2 <- as.data.frame(AverageExpression(seurat.obj, group.by = c("group"), assays = "RNA")[["RNA"]])
name2 <- as.data.frame(colnames(expr2))
name2 <- str_split_fixed(name2$`colnames(expr2)`,"_",2)
name2 <- as.data.frame(name2)
table(name2$V1)
# ASD   BD  CON  MDD PTSD   SZ 
# 1    1    1    1    1    1 
# ASD   BD  CON  MDD PTSD   SZ 
# 2378 2868 2719 1551 2513 2838 


expr <- t(expr)
expr <- as.data.frame(expr)
rownames(expr) <- paste0(rownames(expr),"#",name2$V1)
expr[1:4,1:4]
#          AL627309.1 AL627309.5  LINC01409  LINC01128
# ASD#ASD 0.01300634 0.01156119 0.03528014 0.02901486
# BD#BD   0.01078420 0.03175346 0.09250293 0.06694964
# CON#CON 0.01182909 0.02654205 0.10965935 0.06104129
# MDD#MDD 0.01107854 0.04542202 0.08104490 0.07360202


### 筛选数据前7000的基因
expr <- t(expr)
keep_data <- expr[order(apply(expr, 1, mad),decreasing = T)[1:7000],]
dim(keep_data)
keep_data[1:4,1:4]
keep_data <- t(keep_data)
keep_data <- as.data.frame(keep_data)
dim(keep_data)
keep_data[1:4,1:4]

#           PLXDC2        CALM1       KCNIP4        LSAMP
# ASD#ASD 1.091874e+14 2.950251e+01 1.163843e+04 2.329951e+31
# BD#BD   7.594853e+43 3.662962e+00 1.787717e+61 6.064693e+14
# CON#CON 1.699494e+36 3.024278e+33 4.092916e+32 2.580764e+17
# MDD#MDD 1.046094e+22 1.441165e+34 1.630571e-02 6.868630e+29

### 创建datTraits，包含分组、表型等信息
df <- data.frame(rownames(keep_data))
a <- as.data.frame(str_split_fixed(df$rownames.keep_data.,"#",2))
df <- cbind(df,a)
colnames(df) <- c("AverageExpression_group","Sample","group")

datTraits <- data.frame(row.names = rownames(keep_data),
                        group = df$group)
table(datTraits$group)
# ASD   BD  CON  MDD PTSD   SZ 
#  1    1    1    1    1    1 

### 给分组加上编号
grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
head(datTraits)
#            group groupNo
# ASD#ASD     ASD       1
# BD#BD        BD       2
# CON#CON     CON       3
# MDD#MDD     MDD       4
# PTSD#PTSD  PTSD       5
# SZ#SZ        SZ       6
table(datTraits$group)
table(datTraits$groupNo)

#使用输入文件
datExpr0 <- as.data.frame(keep_data)




############################## 1.判断数据质量 ################################
### 判断数据质量--缺失值
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){ #如果存在异常样本或基因
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) # 异常的基因
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) # 异常的样本
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:从数据中删除违规基因和样本
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

### 绘制样品的系统聚类树
if(T){
  #针对样本做聚类树
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  pdf("step1_Sample dendrogram.pdf",width = 8,height = 6)
  p <- plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
            cex.axis = 1, cex.main = 1,cex.lab=1)
  print(p)
  dev.off()
  # ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)),
                                  colors = rainbow(length(table(datTraits$group))),
                                  signed = FALSE)
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6)
  p2 <- plotDendroAndColors(sampleTree, sample_colors,
                            groupLabels = "trait",
                            cex.dendroLabels = 0.8,
                            marAll = c(1, 4, 3, 1),
                            cex.rowText = 0.01,
                            main = "Sample dendrogram and trait"
  )
  print(p2)
  dev.off()
  # # Plot a line to show the cut
  # abline(h = 23500, col = "red") #根据实际情况而定
  # dev.off()
}


