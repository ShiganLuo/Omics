#######-----------配置--------------#######
# install.packages("stringr")
library(stringr)
library(dplyr)
setwd("~/精神疾病/data/388/Seurat/data/Cellchat/")


###----------数据读取------------
rm(list = ls())  

data1 <- read.table("~/精神疾病/data/388/数据集/SZBDMulti-Seq/SZ1-annotated_matrix.txt.gz",sep="\t", row.names = 1,header = T)
data2 <- read.table("~/精神疾病/data/388/数据集/SZBDMulti-Seq/SZ2-annotated_matrix.txt.gz",sep="\t", row.names = 1,header = T)
data3 <- read.table("~/精神疾病/data/388/数据集/SZBDMulti-Seq/SZ5-annotated_matrix.txt.gz",sep="\t", row.names = 1,header = T)


# 使用bind_rows()合并数据框  
data.combined <- bind_cols(data1, data2, data3)  
# 抽样
set.seed(42)  # 设置随机种子以保证可重复性
random.columns <- sample(ncol(data.combined), 10000)

# 创建新的数据框
data.combined <- data.combined[, random.columns]

# 查看结果
dim(data.combined)
# [1] 33877 10000

path <- "~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat.txt"  
write.table(data.combined, file=path, sep="\t", quote=FALSE, row.names=T)

Seurat_CON <- CreateSeuratObject(counts = data.combined, project = "SZ", min.cells = 3, min.features = 200)
saveRDS(Seurat_CON, "~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat.rds")


###加载以下安装包
library(CellChat)
library(patchwork)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat)
library(tidyverse)
source("/home/maxr/exercise/code/custom_function.R")  #师兄写的包，功能及参数，直接source
rm(list = ls())  
outdir <- "~/精神疾病/data/388/Seurat/data/Cellchat/"

#导入seurat对象
seurat_obj <-
  readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat.rds")
###归一化
seurat_obj <- NormalizeData(seurat_obj)
##基因表达
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")#

##添加行名作为新的一列
seurat_obj@meta.data$cell_type <- rownames(seurat_obj@meta.data)
head(seurat_obj@meta.data)


##分细胞类型
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


##----------------通讯分析--------------##
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
levels(cellchat@idents) # 显示单元格标签的因子级别
# "Astro"  "Endo"   "Exc"    "Inh"    "Micro"  "Oligo"  "OPC"    "others"
groupSize <-
  as.numeric(table(cellchat@idents)) # 每个细胞组中的细胞数
print(groupSize)
# [1]  715   35 4906 2520  389  779  590   66


CellChatDB <-
  CellChatDB.human # 如果是鼠数据，请使用CellChatDB.mouse
####展示数据库情况
# showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # 只需使用默认的CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <-
  subsetData(cellchat) # 将信号基因的表达数据子集化以节省计算成本
future::plan("multisession", workers = 10) # 并行
cellchat <-
  identifyOverExpressedGenes(cellchat) # take a short time
##识别细胞类型交互marker基因
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) # take a short time

#### 计算通信概率并推断细胞间通信网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# 如果某些细胞组中只有少量细胞，则过滤掉细胞间通信
cellchat <- filterCommunication(cellchat, min.cells = 10)
#### 提取推断的细胞间通信网络作为数据帧
#pass
#### 从信号通路水平推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)
#### 计算聚合的细胞间通信网络
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf(file = paste0(outdir, "SZ_netVisual_circle_by_count.pdf"))
par(mar = c(6, 6, 6, 6))  # 增加边距
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
par(mar = c(5, 5, 5, 5))  # 增加边距
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
par(mfrow = c(3, 4), xpd = TRUE)   ###参数影响线的粗细
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

###-----------------multicellchat组间-----------------
library(CellChat)
library(patchwork)
library(cowplot)
rm(list = ls())
# load("    .Rdata")
outdir <- "~/精神疾病/data/388/Seurat/data/Cellchat/ALL/"

cellchat.control  <- readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/CON_cellchat_finished.rds")
cellchat.case  <- readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/SZ/SZ_cellchat_finished.rds")
# cellchat.case2  <- readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/BD/BD_cellchat_finished.rds")
# cellchat.case3  <- readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/MDD/MDD_cellchat_finished.rds")
# cellchat.case4  <- readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/PTSD/PTSD_cellchat_finished.rds")
# cellchat.case5  <- readRDS("~/精神疾病/data/388/Seurat/data/Cellchat/SZ/SZ_cellchat_finished.rds")


# #合并cellchat对象: #对照组在前，比较组在后
object.list <- list(CON = cellchat.control,
                    # ASD = cellchat.case1,
                    # BD = cellchat.case2,
                    # MDD = cellchat.case3,
                    # PTSD = cellchat.case4,
                    SZ = cellchat.case)#####改疾病组

cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),
  cell.prefix = TRUE  # 添加前缀以避免重复
)



####-----图一，呈现两组之间配受体数量的差异----
# 总体比较：通讯数目与通讯强度差异 呈现两组之间配受体数量的差异
# custom_colors <- c("#F28482", "#F7EDE2", "#F6BD60", "#F5CAC3", "#84A59D","#A8D1E7")  # 根据需要定义更多颜色
# gg_combined <- gg_combined +
#   scale_fill_manual(values = custom_colors) +  # 如果使用填充颜色
#   scale_color_manual(values = custom_colors)    # 如果使用边框颜色

gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2, 3, 4, 5, 6))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2, 3, 4, 5, 6), measure = "weight")

gg_combined <- gg1 + gg2
ggsave("1-配受体比较_ALL.pdf", plot = gg_combined, width = 6, height = 3)




####-----图二，网络图Circle Plot绘制两组的，第二个组相较于第一个组细胞通讯的变化
# （low-high,相减的过程），粗细-相互作用强度，红色-在low上调的信号通路，蓝色-下调；----
pdf("2-相互作用强度_ALL.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
# netVisual_diffInteraction(cellchat, weight.scale = T,label.edge = T)
# netVisual_diffInteraction(cellchat, weight.scale = T, label.edge = T,measure = "weight")




####-----图三，热图Heatmap----  
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
p <- gg1 + gg2
pdf("3-热图_SZ.pdf",p,width = 8,height = 4)
print(p)
dev.off()



####-----图四，网络图Circle Plot----
#单独展示两个组之间的通讯强度比较
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



####-----图五，堆叠柱形图2个----
#细胞通讯在两组之间的通讯强度 黑色的是没有显著差异的不看
##基于信息流或互作数对信号通路进行排序
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



####-----图六，总配受体对概率差异气泡图----
table(object.list$SZ@meta)
#                   cell_type
# Idents(seurat_obj) Astro Endo  Exc  Inh Micro Oligo  OPC  others
#               CON   715   35  4906 2520  389   779   590     66
#               ASD   967  516  5036 1200   39  1784   292    166
#                BD  1404  196  2521 1209  460  3651   390    169
#               MDD  1171   34  4184 1431  385  2247   470     78
#              PTSD   694   58  4211 1943  538  2008   480     68
#                SZ   416   15  5863 1494  107  1753   313     39

netVisual_bubble(cellchat,
                 sources.use = c(5),  #2 ##箭头左边的细胞类型数目，第几个细胞类型
                 targets.use = c(1:8),  #3:8 ##箭头→右边的细胞类型数目
                 # signaling = c("PTPRM","SEMA3"),
                 comparison = c(1:2), angle.x = 45)   #comparison = c(3)可修改样本对象，对应day哪天;comparison = c(1:6)全部的day；comparison = c(1, 2)就是group里面有几组 2组 Before Late,
ggsave("6-气泡图1_SZ.pdf",height = 8,width = 6)

netVisual_bubble(cellchat,
                 sources.use = c(1:8),  #2 ##箭头左边的细胞类型数目，第几个细胞类型
                 targets.use = c(5),  #3:8 ##箭头→右边的细胞类型数目
                 # signaling = c("PTPRM","SEMA3"),
                 comparison = c(1:2), angle.x = 45)   #comparison = c(3)可修改样本对象，对应day哪天;comparison = c(1:6)全部的day；comparison = c(1, 2)就是group里面有几组 2组 Before Late,
ggsave("6-气泡图2_SZ.pdf",height = 8,width = 6)
# 


###------------区分上下调配受体对----------
p1 <- netVisual_bubble(cellchat,
                       sources.use = 5,
                       targets.use = c(1:8),
                       comparison = c(1, 2),
                       max.dataset = 2,
                       title.name = "Increased signaling In SZ",
                       angle.x = 45,
                       remove.isolate = T) #Increased为比较组通讯概率更强的配受体对信息
p2 <- netVisual_bubble(cellchat,
                       sources.use = 5,
                       targets.use = c(1:8),
                       comparison = c(1, 2),
                       max.dataset = 1,
                       title.name = "Decreased signaling In SZ",
                       angle.x = 45,
                       remove.isolate = T) #Decreased为对照组通讯概率更强的配受体对信息
p1 + p2

ggsave("配受体上下调_SZ.pdf",width = 10,height = 10)



####-----图七，网络图展示信号通路----
setwd("~/精神疾病/data/388/Seurat/data/Cellchat/SZ/")
load("SZVSCON_cellchat_object.list.RData")
load("SZ_cellchat_merged_.RData")
outdir <- "~/精神疾病/data/388/Seurat/data/Cellchat/SZ/"

#关注通路：APP（CD74）、CD22（PTPRC）、PSAP（GPR37L1）、PTPRM、SPP1（ITGAV、）



pathways.show <- c("APP")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("网络_APP_SZ.pdf",height = 5,width = 8)
par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "circle",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()

# # 一个组中有
# pathways.show <- c("SPP1")
# weight.max <- getMaxWeight(object.list[1],slot.name = c('netP') ,attribute =pathways.show)
# pdf("SPP1_ASD.pdf",height = 5,width = 8)
# par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
# for (i in 1: length(object.list)) {
#   netVisual_aggregate(object.list[[1]],signaling = pathways.show,layout = "circle",
#                       edge.weight.max = weight.max[1],edge.width.max = 10,
#                       signaling.name = paste(pathways.show,names(object.list)[1]))
# }
# dev.off()




####-----图八，和弦图Chord展示信号通路----
# par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
pathways.show <- c("APP")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("和弦_APP_SZ.pdf",height = 15,width = 10)
par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()

# 
pathways.show <- c("APP") 
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
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

# --------其他可视化------
levels(object.list[[1]]@idents)
pathways.show <- c("PECAM1")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("PECAM14_ASD.pdf",height = 20,width = 15)
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(2:3), targets.use = c(2:3), 
                       lab.cex = 0.5, title.name = paste0("Signaling from Macrophage/to Endothelial cell - ", 
                                                          names(object.list)[i]))
}
dev.off()



#-----图九：基因表达组间差异的小提琴图
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CON", "ASD","BD","MDD","PTSD","SZ")) # set factor level
pdf("组间差异APP.pdf",height = 4,width = 10)
plotGeneExpression(cellchat, signaling = "APP", split.by = "datasets", 
                   colors.ggplot = T)+
  scale_fill_manual(values = c('#00AFBB',"#e64b35ff"))
dev.off()

# -------保存运行数据--------
save(object.list, file = "~/精神疾病/data/388/Seurat/data/Cellchat/SZVSCON_cellchat_object.list.RData")
save(cellchat, file = "~/精神疾病/data/388/Seurat/data/Cellchat/SZ_cellchat_merged_.RData")
