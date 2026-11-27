##改代码通过harmony消除异质性
#install.packages("harmony")#未安装需运行
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(harmony)

seurat_obj <- readRDS("~/脑课题/result/Seurat/merged_after_Vln_无MDD.rds")

sort(table(seurat_obj$group), decreasing = T)
table(seurat_obj$group)
seurat_obj_list <- SplitObject(seurat_obj, split.by = 'group')


seurat_obj_list <- lapply(seurat_obj_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 3000)
  return(x)
})
features <- SelectIntegrationFeatures(seurat_obj_list,
                                      nfeatures = 3000)
head(features, n = 20)


exclude_genes <- read.table('~/脑课题/result/Seurat/exclude_genes.txt') |> pull()
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
scale_data <- Reduce(cbind, scale_data_list)
seurat_obj@assays$RNA@scale.data <- scale_data
VariableFeatures(seurat_obj) <- HVGs



seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunHarmony(seurat_obj, "group")#可以指定按组还是按样本

pdf("elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(seurat_obj, ndims = 50)
dev.off()

seurat_obj <- FindNeighbors(seurat_obj,reduction = "harmony",dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) 
library(clustree)
pdf("clustree.pdf", width = 8, height = 10)
clustree(seurat_obj@meta.data,prefix = "RNA_snn_res.")
dev.off()


pdf("umap_tsne.pdf", width = 8, height = 7)
seurat_obj <- RunTSNE(seurat_obj,  dims = 1:15, reduction = "harmony",check_duplicates = FALSE)
DimPlot(seurat_obj,reduction = "tsne",label = TRUE,group.by="group",pt.size = 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15, reduction = "harmony",check_duplicates = FALSE)
DimPlot(seurat_obj,reduction = "umap",label = TRUE,group.by="group",pt.size = 1)
dev.off()


saveRDS(seurat_obj, file = "ALL_harmony_无MDD.rds")

