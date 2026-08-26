##改代码通过harmony消除异质性
#install.packages("harmony")#未安装需运行
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(harmony)
#### 亚群重聚类 ----
#endo <- subset(seurat.obj, idents ='endothelial cells', invert=F)
#saveRDS(endo,
#        file = "./results/seurat/merge/0922/endo.rds")

seurat_obj <- readRDS("~/脑课题/result/Seurat/merged_after_Vln_无MDD.rds")

# 分别找高变基因 ----
sort(table(seurat_obj$group), decreasing = T)
table(seurat_obj$group)
# CON   ASD    BD  PTSD    SZ 
# 10536 16326  5524  3523  5983  
seurat_obj_list <- SplitObject(seurat_obj, split.by = 'group')


seurat_obj_list <- lapply(seurat_obj_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 3000)
  return(x)
})
features <- SelectIntegrationFeatures(seurat_obj_list,
                                      nfeatures = 3000)
head(features, n = 20)


# --------------从 Seurat 对象列表中排除一些基因---------------
exclude_genes <- read.table('~/脑课题/result/Seurat/exclude_genes.txt') |> pull()
# #HVGs <- setdiff(features, exclude_genes)
HVGs <- features[!(features %in% exclude_genes)][1:1000]
seurat_obj_list <- lapply(seurat_obj_list, function(x) {
  VariableFeatures(x) <- HVGs
  x <- ScaleData(x, features = HVGs,
                 vars.to.regress = c("nCount_RNA", "percent.mt")
  )
  return(x)
})

seurat_obj <- merge(seurat_obj_list[[1]],
                    y = seurat_obj_list[-1],
                    merge.data = TRUE)

scale_data_list <- lapply(seurat_obj_list, function(x) {
  z <- x@assays$RNA@scale.data
  print(dim(z))
  return(z)
})
# [1]  1000 10536
# [1]  1000 16326
# [1] 1000 5524
# [1] 1000 3523
# [1] 1000 5983
scale_data <- Reduce(cbind, scale_data_list)
seurat_obj@assays$RNA@scale.data <- scale_data
VariableFeatures(seurat_obj) <- HVGs


# seurat_obj <- readRDS("./data/merge/sub_harmony.rds")

##PCA
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunHarmony(seurat_obj, "group")#可以指定按组还是按样本

# 肘部图确定最佳dims
pdf("elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(seurat_obj, ndims = 50)
dev.off()

##运行之后在reductions下多出harmony值，后续的聚类及可视化调用harmony
seurat_obj <- FindNeighbors(seurat_obj,reduction = "harmony",dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) 
library(clustree)
pdf("clustree.pdf", width = 8, height = 10)
clustree(seurat_obj@meta.data,prefix = "RNA_snn_res.")
dev.off()

# seurat_obj <- FindClusters(seurat_obj,
#                        resolution = 0.1) # 分辨率值越大，cluster 越多0~1之间

pdf("umap_tsne.pdf", width = 8, height = 7)
seurat_obj <- RunTSNE(seurat_obj,  dims = 1:15, reduction = "harmony",check_duplicates = FALSE)
DimPlot(seurat_obj,reduction = "tsne",label = TRUE,group.by="group",pt.size = 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15, reduction = "harmony",check_duplicates = FALSE)
DimPlot(seurat_obj,reduction = "umap",label = TRUE,group.by="group",pt.size = 1)
dev.off()

##到这里就可以完成批次校正

saveRDS(seurat_obj, file = "ALL_harmony_无MDD.rds")

