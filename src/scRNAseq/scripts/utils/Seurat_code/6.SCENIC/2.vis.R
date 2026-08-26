#### 导入包 ----
library(SCENIC)
library(AUCell)
library(Seurat)

#### 导入数据 ----
# auc得分矩阵
seurat.obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/SZ.rds")
# 导入auc数据
regulonAUC <- importAUCfromText("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/auc_g_mtx.csv")
regulonAUC <-
  regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
cellInfo <- seurat.obj@meta.data

a <- getAUC(regulonAUC)
a <- t(a)


# 将得分矩阵添加进seurat
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

# 检查不匹配的转录因子
mismatched_tf <- tf.seleted[!tf.seleted %in% rownames(a)]
print(mismatched_tf)

#### 1. 按照细胞类型进行热图可视化 ----
regulonActivity_byCellType <-
  sapply(split(rownames(cellInfo), Idents(seurat.obj)),
         function(cells)
           rowMeans(getAUC(regulonAUC)[, cells]))

# regulon <- regulonActivity_byCellType[grep("MAFK(+)|FOXD3(+)|ALX3(+)|ETV5(+)|RARA(+)|FOXP1(+)|CTCF(+)|NR2F1(+)",x = rownames(regulonActivity_byCellType),ignore.case = T),]

# tf.seleted <- names(regulon) #指定转录因子
# tf.seleted <- c("ARNTL(+)","ATF4(+)","ELK4(+)","NFKB1(+)","IRF7(+)","THRB(+)")

pheatmap(
  regulonActivity_byCellType[tf.seleted,],
  scale = "row",
  # annotation_col=meta,
  color =  colorRampPalette(c("navy", "white", "firebrick3"))(50),
  #breaks=seq(-value.max, value.max, length.out = 100),
  # breaks = seq(1,100,length.out = 10),   #加了一行
  treeheight_row = 10, #聚类树行高
  treeheight_col = 10, #聚类树列高
  border_color = "white", #单元格边框颜色
  cellheight = 5, #单元格高度
  cellwidth = 10, #单元格宽度
  fontsize = 5, #字体大小
  angle_col = "45", #列名旋转角度
  cluster_cols = F,
  cluster_rows = F, #是否出现系统发生树
  #  filename = "/home/maxr/maxr/exercise/results/SMC_heatmap_by_celltype1.pdf"
)

#dev.off() #在plot区显示图片

#### 2. 按照细胞可视化 ----
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
  # scale = "row",
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


#### 3. 转录因子表达水平可视化 ----
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

#### 4. 转录因子活性在umap或tsne投影可视化 ----
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
