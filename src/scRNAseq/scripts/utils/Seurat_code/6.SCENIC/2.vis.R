#### 导入包 ----
library(SCENIC)
library(AUCell)
library(Seurat)

seurat.obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/SZ.rds")
regulonAUC <- importAUCfromText("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/auc_g_mtx.csv")
regulonAUC <-
  regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
cellInfo <- seurat.obj@meta.data

a <- getAUC(regulonAUC)
a <- t(a)


seurat.obj[["scenic"]] <-
  CreateAssayObject(counts = a)

scaled <-
  t(scale(
    t(getAUC(regulonAUC)),
    center = T,
    scale = T
  ))
tf.seleted <- c("ATF4", "E2F1", "CEBPE","EGR1",
                "ELK1", "ETV5", "NFKB1","TEF")

mismatched_tf <- tf.seleted[!tf.seleted %in% rownames(a)]
print(mismatched_tf)

regulonActivity_byCellType <-
  sapply(split(rownames(cellInfo), Idents(seurat.obj)),
         function(cells)
           rowMeans(getAUC(regulonAUC)[, cells]))



pheatmap(
  regulonActivity_byCellType[tf.seleted,],
  scale = "row",
  color =  colorRampPalette(c("navy", "white", "firebrick3"))(50),
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = "white", #单元格边框颜色
  cellheight = 5,
  cellwidth = 10,
  fontsize = 5,
  angle_col = "45", #列名旋转角度
  cluster_cols = F,
  cluster_rows = F,
)


regulonActivity_byCell <- a
seurat.obj$cell_id <- colnames(seurat.obj)
annotation_col <-
  as.data.frame(seurat.obj@meta.data[, c("group")])
rownames(annotation_col) <- rownames(seurat.obj@meta.data)

scaled <-
  t(scale(
    t(regulonActivity_byCell),
    center = T,
    scale = T
  ))

scaled[scaled >= 2] <- 2
scaled[scaled <= -2] <- -2
pheatmap(
  scaled[tf.seleted, ],
  annotation_col = annotation_col,
  color =  colorRampPalette(c("navy", "white", "firebrick3"))(50),
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = "white",
  cellwidth = 0.2,
  cellheight = 5,
  fontsize = 5,
  angle_col = "45",
  cluster_rows = F,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = F,
  height = 1,
  filename = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/SZ_heatmap_by_cell.pdf"
)


DefaultAssay(seurat.obj) <- "RNA"
tf <- "ETV6"
p1 <- FeaturePlot(
  seurat.obj,
  features = tf,
  reduction = "umap",
  pt.size = 2,
  order = T,
  min.cutoff = 0
)
p1

DefaultAssay(seurat.obj) <- "scenic"
tf <- "ETV6"
p2 <- FeaturePlot(
  seurat.obj,
  features = tf,
  reduction = 'umap',
  min.cutoff = 0,
  cols = c("lightgrey" , "#DE1F1F"),
  pt.size = 2,
  order = T
)
p2

p <- p1 + p2 + plot_layout(ncol = 2)
p
ggsave("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/feature_plot1.pdf",plot = p,width = 11,height = 5)
