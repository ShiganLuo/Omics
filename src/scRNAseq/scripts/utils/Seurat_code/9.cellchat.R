#######-----------配置--------------#######
# install.packages("stringr")
library(stringr)
library(dplyr)
setwd("~/精神疾病/data/388/Seurat/data/Cellchat/")


rm(list = ls())  

data1 <- read.table("~/精神疾病/data/388/数据集/SZBDMulti-Seq/SZ1-annotated_matrix.txt.gz",sep="\t", row.names = 1,header = T)
data2 <- read.table("~/精神疾病/data/388/数据集/SZBDMulti-Seq/SZ2-annotated_matrix.txt.gz",sep="\t", row.names = 1,header = T)
data3 <- read.table("~/精神疾病/data/388/数据集/SZBDMulti-Seq/SZ5-annotated_matrix.txt.gz",sep="\t", row.names = 1,header = T)


data.combined <- bind_cols(data1, data2, data3)  
set.seed(42)
random.columns <- sample(ncol(data.combined), 10000)

data.combined <- data.combined[, random.columns]

dim(data.combined)

path <- "~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat.txt"  
write.table(data.combined, file=path, sep="\t", quote=FALSE, row.names=T)

Seurat_CON <- CreateSeuratObject(counts = data.combined, project = "SZ", min.cells = 3, min.features = 200)
saveRDS(Seurat_CON, "~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat.rds")


library(CellChat)
library(patchwork)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat)
library(tidyverse)
source("/home/maxr/exercise/code/custom_function.R")  #师兄写的包，功能及参数，直接source
rm(list = ls())  
outdir <- "~/精神疾病/data/388/Seurat/data/Cellchat/"

seurat_obj <-
  readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat.rds")
seurat_obj <- NormalizeData(seurat_obj)
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")#

seurat_obj@meta.data$cell_type <- rownames(seurat_obj@meta.data)
head(seurat_obj@meta.data)


seurat_obj@meta.data$cell_type = "others"

seurat_obj@meta.data[grep("OPC",names(seurat_obj@active.ident)),]$cell_type = "OPC"
seurat_obj@meta.data[grep("Oligo",names(seurat_obj@active.ident)),]$cell_type = "Oligo"
seurat_obj@meta.data[grep("Astro",names(seurat_obj@active.ident)),]$cell_type = "Astro"
seurat_obj@meta.data[grep("Micro",names(seurat_obj@active.ident)),]$cell_type = "Micro"
seurat_obj@meta.data[grep("Endo",names(seurat_obj@active.ident)),]$cell_type = "Endo"
seurat_obj@meta.data[grep("IT|CT|ET|NP|L6b",names(seurat_obj@active.ident)),]$cell_type = "Exc"
seurat_obj@meta.data[grep("Sst|Pvalb|Chan|Pax6|Lamp|Sncg|Vip",names(seurat_obj@active.ident)),]$cell_type = "Inh"

View(seurat_obj@meta.data)

saveRDS(seurat_obj,"~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat_Normalized.rds")


rm(list = ls()) 
outdir <- "~/精神疾病/data/388/Seurat/data/Cellchat/"

seurat_obj <-
  readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat_Normalized.rds")

data.input <- seurat_obj@assays$RNA@data
meta <- as.data.frame(Idents(seurat_obj))
meta$cell_type = seurat_obj@meta.data$cell_type


cellchat <-
  createCellChat(object = data.input,
                 meta = meta,
                 group.by = "cell_type")
cellchat <- addMeta(cellchat, meta = meta)

cellchat <-
  setIdent(cellchat, ident.use = "cell_type") # 将“cell_type”设置为默认单元格标识
levels(cellchat@idents)
groupSize <-
  as.numeric(table(cellchat@idents))
print(groupSize)


CellChatDB <-
  CellChatDB.human

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <-
  subsetData(cellchat)
future::plan("multisession", workers = 10) # 并行
cellchat <-
  identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf(file = paste0(outdir, "SZ_netVisual_circle_by_count.pdf"))
par(mar = c(6, 6, 6, 6))
a = netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  
  
  title.name = "Number of interactions"
)
print(a)
dev.off()


pdf(file = paste0(outdir, "SZ_netVisual_circle_by_weight.pdf"))
par(mar = c(5, 5, 5, 5))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
  
)
dev.off()

pdf(file = paste0(outdir, "SZ_netVisual_circle_by_cell_type.pdf"),width = 10,height = 10)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <-
    matrix(
      0,
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
  mat2[i,] <- mat[i,]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = T,
    arrow.width = 0.2,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()

saveRDS(cellchat, file = paste0(outdir, "SZ_cellchat_finished.rds"))

library(CellChat)
library(patchwork)
library(cowplot)
rm(list = ls())
outdir <- "~/精神疾病/data/388/Seurat/data/Cellchat/ALL/"

cellchat.control  <- readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/CON_cellchat_finished.rds")
cellchat.case  <- readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/SZ/SZ_cellchat_finished.rds")


object.list <- list(CON = cellchat.control,
                    SZ = cellchat.case)#####改疾病组

cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),
  cell.prefix = TRUE
)




gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2, 3, 4, 5, 6))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2, 3, 4, 5, 6), measure = "weight")

gg_combined <- gg1 + gg2
ggsave("1-配受体比较_ALL.pdf", plot = gg_combined, width = 6, height = 3)




pdf("2-相互作用强度_ALL.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()




gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
p <- gg1 + gg2
pdf("3-热图_SZ.pdf",p,width = 8,height = 4)
print(p)
dev.off()



weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("4-网络图_SZ.pdf",height = 5,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale = T,
                   label.edge= F,
                   edge.weight.max = weight.max[2],
                   edge.width.max = 5, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()



load("ALL_cellchat_object.list.RData")
load("ALL_cellchat_merged_.RData")

gg1 <- rankNet(cellchat,
               mode = "comparison", stacked = T,
               comparison = c(1,2,3,4,5,6),#展示所有组加这一行
               do.stat = TRUE) +coord_flip()
gg2 <- rankNet(cellchat,
               mode = "comparison", stacked = F,
               comparison = c(1,2,3,4,5,6),#展示所有组加这一行
               do.stat = TRUE) +coord_flip()

gg1+gg2
ggsave("5-堆叠柱形图_ALL.pdf", plot = gg1+gg2, width = 12,height = 10)



table(object.list$SZ@meta)

netVisual_bubble(cellchat,
                 sources.use = c(5),
                 targets.use = c(1:8),
                 comparison = c(1:2), angle.x = 45)
ggsave("6-气泡图1_SZ.pdf",height = 8,width = 6)

netVisual_bubble(cellchat,
                 sources.use = c(1:8),
                 targets.use = c(5),
                 comparison = c(1:2), angle.x = 45)
ggsave("6-气泡图2_SZ.pdf",height = 8,width = 6)


p1 <- netVisual_bubble(cellchat,
                       sources.use = 5,
                       targets.use = c(1:8),
                       comparison = c(1, 2),
                       max.dataset = 2,
                       title.name = "Increased signaling In SZ",
                       angle.x = 45,
                       remove.isolate = T)
p2 <- netVisual_bubble(cellchat,
                       sources.use = 5,
                       targets.use = c(1:8),
                       comparison = c(1, 2),
                       max.dataset = 1,
                       title.name = "Decreased signaling In SZ",
                       angle.x = 45,
                       remove.isolate = T)
p1 + p2

ggsave("配受体上下调_SZ.pdf",width = 10,height = 10)



setwd("~/精神疾病/data/388/Seurat/data/Cellchat/SZ/")
load("SZVSCON_cellchat_object.list.RData")
load("SZ_cellchat_merged_.RData")
outdir <- "~/精神疾病/data/388/Seurat/data/Cellchat/SZ/"




pathways.show <- c("APP")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("网络_APP_SZ.pdf",height = 5,width = 8)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "circle",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()





pathways.show <- c("APP")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("和弦_APP_SZ.pdf",height = 15,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()

pathways.show <- c("APP") 
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",
                                                  names(object.list)[i]))
}
p <- ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pdf("热图_APP_SZ.pdf",height = 5,width = 10)
print(p)
dev.off()

levels(object.list[[1]]@idents)
pathways.show <- c("PECAM1")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("PECAM14_ASD.pdf",height = 20,width = 15)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(2:3), targets.use = c(2:3), 
                       lab.cex = 0.5, title.name = paste0("Signaling from Macrophage/to Endothelial cell - ", 
                                                          names(object.list)[i]))
}
dev.off()



cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CON", "ASD","BD","MDD","PTSD","SZ")) # set factor level
pdf("组间差异APP.pdf",height = 4,width = 10)
plotGeneExpression(cellchat, signaling = "APP", split.by = "datasets", 
                   colors.ggplot = T)+
  scale_fill_manual(values = c('#00AFBB',"#e64b35ff"))
dev.off()

save(object.list, file = "~/精神疾病/data/388/Seurat/data/Cellchat/SZVSCON_cellchat_object.list.RData")
save(cellchat, file = "~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat_merged_.RData")
