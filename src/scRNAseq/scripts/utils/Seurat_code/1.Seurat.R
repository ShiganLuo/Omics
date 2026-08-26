library(Seurat)
##———— ——————————————————Seurat——————————————————————————###

setwd("/home/lsg/Data/glioblastoma") 

# ##————————先Seurat后merge合并————————————###

# merged_data_CON <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/CON_Micro_merged_data1.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# merged_data_ASD1 <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/ASD_Micro_merged_data1.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# merged_data_ASD2 <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/ASD_Micro_merged_data2.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# merged_data_ASD3 <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/ASD_Micro_merged_data3.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# merged_data_BD1 <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/BD_Micro_merged_data1.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# merged_data_BD2 <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/BD_Micro_merged_data2.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# # merged_data_MDD <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/MDD_Micro_merged_data1.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# merged_data_PTSD <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/PTSD_Micro_merged_data1.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# merged_data_SZ1 <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/SZ_Micro_merged_data1.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# merged_data_SZ2 <- read.table("/home/lsg/Data/glioblastoma/data/merged_data/SZ_Micro_merged_data2.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")



# # 假设pbmc.data已经是一个存在的DataFrame，并且其列名是barcodes  
# # 使用paste0来添加前缀，这里不需要分隔符
# new_colnames <- paste0("CON-", colnames(merged_data_CON))  
# colnames(merged_data_CON) <- new_colnames  

# new_colnames <- paste0("ASD-", colnames(merged_data_ASD1))  
# colnames(merged_data_ASD1) <- new_colnames  # 更新DataFrame的列名 
# new_colnames <- paste0("ASD-", colnames(merged_data_ASD2))  
# colnames(merged_data_ASD2) <- new_colnames 
# new_colnames <- paste0("ASD-", colnames(merged_data_ASD3))  
# colnames(merged_data_ASD3) <- new_colnames 


# new_colnames <- paste0("BD-", colnames(merged_data_BD1))  
# colnames(merged_data_BD1) <- new_colnames  
# new_colnames <- paste0("BD-", colnames(merged_data_BD2))  
# colnames(merged_data_BD2) <- new_colnames  

# # new_colnames <- paste0("MDD_", colnames(merged_data_MDD))  
# # colnames(merged_data_MDD) <- new_colnames  

# new_colnames <- paste0("PTSD-", colnames(merged_data_PTSD))  
# colnames(merged_data_PTSD) <- new_colnames  

# new_colnames <- paste0("SZ-", colnames(merged_data_SZ1))  
# colnames(merged_data_SZ1) <- new_colnames 
# new_colnames <- paste0("SZ-", colnames(merged_data_SZ2))  
# colnames(merged_data_SZ2) <- new_colnames 


# CON <- CreateSeuratObject(counts = merged_data_CON, project = "CON", min.cells = 3, min.features = 200)
# ASD1 <- CreateSeuratObject(counts = merged_data_ASD1, project = "ASD", min.cells = 3, min.features = 200)
# ASD2 <- CreateSeuratObject(counts = merged_data_ASD2, project = "ASD", min.cells = 3, min.features = 200)
# ASD3 <- CreateSeuratObject(counts = merged_data_ASD3, project = "ASD", min.cells = 3, min.features = 200)
# BD1 <- CreateSeuratObject(counts = merged_data_BD1, project = "BD", min.cells = 3, min.features = 200)
# BD2 <- CreateSeuratObject(counts = merged_data_BD2, project = "BD", min.cells = 3, min.features = 200)
# # MDD <- CreateSeuratObject(counts = merged_data_MDD, project = "MDD", min.cells = 3, min.features = 200)
# PTSD <- CreateSeuratObject(counts = merged_data_PTSD, project = "PTSD", min.cells = 3, min.features = 200)
# SZ1 <- CreateSeuratObject(counts = merged_data_SZ1, project = "SZ", min.cells = 3, min.features = 200)
# SZ2 <- CreateSeuratObject(counts = merged_data_SZ2, project = "SZ", min.cells = 3, min.features = 200)

# # merged<-merge(CON,c(ASD1,ASD2,ASD3,BD1,BD2,MDD,PTSD,SZ1,SZ2))
# merged<-merge(CON,c(ASD1,ASD2,ASD3,BD1,BD2,PTSD,SZ1,SZ2))

# # 重新设置'group'的水平顺序 
# merged$group <- merged$orig.ident
# # new_levels <- c("CON", "ASD", "BD", "MDD","PTSD","SZ") 
# new_levels <- c("CON", "ASD", "BD", "PTSD","SZ") 
# merged$group <- factor(merged$group, levels = new_levels)  

# table(merged$group)
# # CON  ASD   BD  MDD PTSD   SZ 
# # 10633 16715  5778  1551  3575  6477 

# saveRDS(merged,
#         file = "/home/lsg/Data/glioblastoma/data/Seurat/merged_only_无MDD.rds")


merged <- readRDS("/home/lsg/Data/glioblastoma/data/Seurat/merged_only_无MDD.rds")
str(merged)
##提取线粒体基因
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
#merged[["percent.rp"]] = PercentageFeatureSet(merged, pattern = "^RP[SL][[:digit:]]") 提取核糖体蛋白基因

##
all.genes <- rownames(merged)
print(all.genes)
merged <- ScaleData(merged, features = all.genes)
# 假设 seurat_object 是你的 SeuratObject 对象  
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
# merged <- RunUMAP(merged, features = VariableFeatures(object = merged))#, dims = 1:10, dims = 0

# 检查数据中是否存在非有限值
print(sum(!is.finite(merged@assays$RNA@data)))

# UMAP
merged <- RunUMAP(merged, dims = 1:10)
pdf("umap_无MDD.pdf", width = 7, height = 6)
DimPlot(merged, reduction = "umap", label = T)
dev.off()






##小提琴图##
merged <- readRDS("/home/lsg/Data/glioblastoma/data/Seurat/merged_only_无MDD.rds")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "group")

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

########根据小提琴图过滤######
# merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

# 然后根据 group 和 nFeature_RNA 的条件进行过滤
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

merged <- subset(merged,
                 (group == "CON" & nFeature_RNA < 3000) |
                   (group == "ASD" & nFeature_RNA < 5000) |
                   (group == "BD" & nFeature_RNA < 2500) |
                   # (group == "MDD" & nFeature_RNA < 2000) |
                   (group == "PTSD" & nFeature_RNA < 1900) |
                   (group == "SZ" & nFeature_RNA < 2500) &
                   nFeature_RNA > 200
)

pdf("vln_plot_无MDD_过滤.pdf", width = 10, height = 5)
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "group")
plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

saveRDS(merged,file = "~/脑课题/result/Seurat/merged_after_Vln_无MDD.rds")

#####————————————————————————————————#####标准化#####————————————————————————————————————#####

##---------------- 方法一 ----------------##
merged <- readRDS("~/脑课题/result/Seurat/merged_after_Vln.rds")

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(merged, file = "~/脑课题/result/Seurat/merged_Normalized.rds")

##找高变基因
# merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
# top10 <- head(VariableFeatures(merged), 10)
# top10


# all.genes <- rownames(merged)
# merged <- ScaleData(merged, features = all.genes)
# 
# merged <- RunPCA(merged, features = VariableFeatures(object = merged))
# print(merged[["pca"]], dims = 1:5, nfeatures = 5)
# 
# VizDimLoadings(merged, dims = 1:2, reduction = "pca")
# DimPlot(merged, reduction = "pca")
# 
# 
# DimHeatmap(merged, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(merged, dims = 1:15, cells = 500, balanced = TRUE)



##---------------- 方法二 ----------------##
#### 样本整合 ----
library(glmGamPoi)
options(future.globals.maxSize = 1 * 1024^2 * 1024)  # 1 GB


merged <- readRDS("~/脑课题/result/Seurat/merged_after_Vln.rds")


seurat.list <- SplitObject(merged, split.by = "orig.ident")
seurat.list <- lapply(
  X = seurat.list,
  FUN = function(x) {
    x <- SCTransform(
      x,
      method = "glmGamPoi",
      vars.to.regress = c("percent.mt"),
      verbose = T
    )
  }
)

seurat.list <- lapply(seurat.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)  # 对于 SCT，常用的选择方法是 "vst"
  return(x)
})


features <-
  SelectIntegrationFeatures(object.list = seurat.list,
                            nfeatures = 3000)
seurat.list <-
  PrepSCTIntegration(object.list = seurat.list,
                     anchor.features = features)

immune.anchors <- FindIntegrationAnchors(
  object.list = seurat.list,
  normalization.method = "SCT",
  anchor.features = features
)

merged  <- IntegrateData(anchorset = immune.anchors,
                         normalization.method = "SCT")
saveRDS(merged,
        file = "~/脑课题/Seurat/Scted_new.rds")


####---------------- 降维聚类分群 --------------------####
# PCA降维
merged <- RunPCA(merged)
# 肘部图确定最佳dims
ElbowPlot(merged, ndims = 50)
# t-SNE
merged <- RunTSNE(merged, dims = 1:10)
# UMAP
merged <- RunUMAP(merged, dims = 1:10)

# 聚类分析
# merged <- FindNeighbors(merged, dims = 1:15)
# merged <- FindClusters(merged,resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
# clustree(merged@meta.data,prefix = "integrated_snn_res.")
# merged <- FindClusters(merged,
#                        resolution = 0.1) # 分辨率值越大，cluster 越多0~1之间
DimPlot(merged,
        reduction = "umap",
        label = T)

DimPlot(merged,
        reduction = "tsne", # tsne, umap, pca
        label = T)
DimPlot(merged,reduction = "umap",label = TRUE,pt.size = 1.0,group.by = "group")
DimPlot(merged,reduction = "umap",label = TRUE,split.by="group",pt.size = 1,ncol = 3)

# # 差异分析（寻找差异表达的maker）
# markers <- FindAllMarkers(
#   merged,
#   only.pos = TRUE,#只选取正数（正数上调，负数下调）
#   min.pct = 0.25,
#   logfc.threshold = 0.25
# )
# write.csv(markers,"markers.csv")
# top.markers <- markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# 
# write.table(top.markers[top.markers$cluster == "2",]$gene,
#             quote = F, row.names = F)
# 
# DotPlot(merged, features = top.markers$gene,
#         col.min = 0)+theme(axis.text.x = element_text(angle = 45,vjust = 0))
# 
# DoHeatmap(merged, features = top.markers$gene)




# #### 细胞注释 ----
# features <- c("CD3D",# T cells
#               "GNLY",# NK cells
#               "LYZ",# Myeloid cells
#               "CD24",# B cells
#               "COL1A1",# Fibroblasts
#               "ACTA2",# SMC
#               "PDGFRB",  #周细胞
#               "CLDN5") # Endothelial cells
# FeaturePlot(
#   merged,
#   features = features,
#   reduction = "umap",
#   order = T,
#   ncol = 4,
#   min.cutoff = 0
# )
# 
# DotPlot(merged, features = features, col.min = 0)
# 
# 
# # 分配细胞名字
# new.cluster.ids <- c(
#   "Smooth muscle cells",
#   # 00
#   "Endothelial Cells",
#   # 01
#   "Pericytes",
#   # 02
#   "Smooth muscle cells",
#   # 03
#   "Myeloid cells",
#   # 04
#   "Fibroblasts",
#   # 05
#   "T cells"
#   # 06
# )
# names(new.cluster.ids) <- levels(merged)
# merged <- RenameIdents(merged,
#                        new.cluster.ids)
# 
# merged$cell_type <- Idents(merged)
# 
# DimPlot(merged,
#         label = T,
#         repel = T,
#         reduction = "umap") +
#   NoLegend()










#### 保存结果 ----
saveRDS(merged,
        file = "~/精神疾病/data/388/Seurat/data/merged_SCTed.rds")


