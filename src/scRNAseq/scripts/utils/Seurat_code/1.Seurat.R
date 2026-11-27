library(Seurat)
##———— ——————————————————Seurat——————————————————————————###

setwd("/home/lsg/Data/glioblastoma") 
























merged <- readRDS("/home/lsg/Data/glioblastoma/data/Seurat/merged_only_无MDD.rds")
str(merged)
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

all.genes <- rownames(merged)
print(all.genes)
merged <- ScaleData(merged, features = all.genes)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))

print(sum(!is.finite(merged@assays$RNA@data)))

merged <- RunUMAP(merged, dims = 1:10)
pdf("umap_无MDD.pdf", width = 7, height = 6)
DimPlot(merged, reduction = "umap", label = T)
dev.off()






merged <- readRDS("/home/lsg/Data/glioblastoma/data/Seurat/merged_only_无MDD.rds")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "group")

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

merged <- subset(merged,
                 (group == "CON" & nFeature_RNA < 3000) |
                   (group == "ASD" & nFeature_RNA < 5000) |
                   (group == "BD" & nFeature_RNA < 2500) |
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


merged <- readRDS("~/脑课题/result/Seurat/merged_after_Vln.rds")

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(merged, file = "~/脑课题/result/Seurat/merged_Normalized.rds")






library(glmGamPoi)
options(future.globals.maxSize = 1 * 1024^2 * 1024)


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


merged <- RunPCA(merged)
ElbowPlot(merged, ndims = 50)
merged <- RunTSNE(merged, dims = 1:10)
merged <- RunUMAP(merged, dims = 1:10)

DimPlot(merged,
        reduction = "umap",
        label = T)

DimPlot(merged,
        reduction = "tsne", # tsne, umap, pca
        label = T)
DimPlot(merged,reduction = "umap",label = TRUE,pt.size = 1.0,group.by = "group")
DimPlot(merged,reduction = "umap",label = TRUE,split.by="group",pt.size = 1,ncol = 3)















saveRDS(merged,
        file = "~/精神疾病/data/388/Seurat/data/merged_SCTed.rds")


