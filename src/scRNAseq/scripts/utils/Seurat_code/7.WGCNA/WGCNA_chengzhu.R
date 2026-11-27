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
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 
getwd()



seurat.obj <- readRDS("/home/maxr/精神疾病/data/388/Seurat/data/SCTed_new.rds")


n_cells_counts <- ncol(seurat.obj@assays$RNA@counts)

n_cells_meta <- nrow(seurat.obj@meta.data)







Idents(seurat.obj) <- "group"
table(Idents(seurat.obj))


seurat.obj = NormalizeData(object = seurat.obj, assay = "RNA")
expr <- AverageExpression(seurat.obj, group.by = "group", assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
dim(expr)

seurat.obj@meta.data$sample = rownames(seurat.obj@meta.data)

expr2 <- as.data.frame(AverageExpression(seurat.obj, group.by = c("group"), assays = "RNA")[["RNA"]])
name2 <- as.data.frame(colnames(expr2))
name2 <- str_split_fixed(name2$`colnames(expr2)`,"_",2)
name2 <- as.data.frame(name2)
table(name2$V1)


expr <- t(expr)
expr <- as.data.frame(expr)
rownames(expr) <- paste0(rownames(expr),"#",name2$V1)
expr[1:4,1:4]


expr <- t(expr)
keep_data <- expr[order(apply(expr, 1, mad),decreasing = T)[1:7000],]
dim(keep_data)
keep_data[1:4,1:4]
keep_data <- t(keep_data)
keep_data <- as.data.frame(keep_data)
dim(keep_data)
keep_data[1:4,1:4]


df <- data.frame(rownames(keep_data))
a <- as.data.frame(str_split_fixed(df$rownames.keep_data.,"#",2))
df <- cbind(df,a)
colnames(df) <- c("AverageExpression_group","Sample","group")

datTraits <- data.frame(row.names = rownames(keep_data),
                        group = df$group)
table(datTraits$group)

grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
head(datTraits)
table(datTraits$group)
table(datTraits$groupNo)

datExpr0 <- as.data.frame(keep_data)




gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

if(T){
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  pdf("step1_Sample dendrogram.pdf",width = 8,height = 6)
  p <- plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
            cex.axis = 1, cex.main = 1,cex.lab=1)
  print(p)
  dev.off()
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)),
                                  colors = rainbow(length(table(datTraits$group))),
                                  signed = FALSE)
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
}


