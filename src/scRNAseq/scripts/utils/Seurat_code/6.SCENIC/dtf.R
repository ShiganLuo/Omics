setwd("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/VS")

CON_TFS<- read.table("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/CON/tfs_targer.tsv", header=T, sep=",")
ASD_TFS<- read.table("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/ASD/tfs_targer.tsv", header=T, sep=",")
BD_TFS<- read.table("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/BD/tfs_targer.tsv", header=T, sep=",")
MDD_TFS<- read.table("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/MDD/tfs_targer.tsv", header=T, sep=",")
PTSD_TFS<- read.table("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/PTSD/tfs_targer.tsv", header=T, sep=",")
SZ_TFS<- read.table("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/tfs_targer.tsv", header=T, sep=",")

tf_to_keep0 <- CON_TFS$tf
tf_to_keep1 <- ASD_TFS$tf
tf_to_keep2 <- BD_TFS$tf
tf_to_keep3 <- MDD_TFS$tf
tf_to_keep4 <- PTSD_TFS$tf
tf_to_keep5 <- SZ_TFS$tf

# 取并集
union_tfs01 <- union(tf_to_keep0, tf_to_keep1)
union_tfs02 <- union(tf_to_keep0, tf_to_keep2)
union_tfs03 <- union(tf_to_keep0, tf_to_keep3)
union_tfs04 <- union(tf_to_keep0, tf_to_keep4)
union_tfs05 <- union(tf_to_keep0, tf_to_keep5)

scRNA0<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_CON.rds")#读取rds文件
scRNA1<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_ASD.rds")#读取rds文件
scRNA_filtered0 <- subset(scRNA0, features = union_tfs01)
scRNA_filtered1 <- subset(scRNA1, features = union_tfs01)
 
# 加载Seurat包
library(Seurat)

# 假设scRNA_filtered0和scRNA_filtered1都是Seurat对象
# 使用Seurat的`merge`函数来合并这两个对象

combined_seurat <- merge(scRNA_filtered0, scRNA_filtered1)

# 查看合并后的对象
print(combined_seurat)

saveRDS(combined_seurat,"ASDVSCON.rds")

seurat.obj<-combined_seurat

table(seurat.obj$group)#查看疾病与正常分组(IVL,Normal)
table(seurat.obj$cell_type)#查看细胞类型分组()

#### 寻找差异基因 ----
# 将细胞类型及刺激状态作为分组
Idents(seurat.obj) <- "group"

# 查看当前分组
table(Idents(seurat.obj))

#使用FindMarkers函数寻找差异表达基因
DefaultAssay(seurat.obj) <- "RNA"

##_____________________________________________________________________________________________

markers <- FindMarkers(seurat.obj,
                       ident.1 = "ASD",
                       ident.2 = "CON",
                       logfc.threshold = 0)



markers$celltype <- "ASDVSCON"
markers$gene <- rownames(markers)
#### 保存结果 ----
# 差异基因
write.table(
  markers,
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/ASDVSCON_deg.txt",
  quote = F,
  sep = ",",
  row.names = F
)

##————————————————————————火山图——————————————————————————————————————————————————————————————————————————————————————————————
#### 导入包 ----
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(dplyr)

#### 导入数据 ----
##markers <- read.table("~/brain/差异/output/deg/ASDVSCON_deg.txt", sep = ",", header = T)

#### 定义差异基因 ----
# 设置阈值
logFC.cutoff <- 0.25
pvalue.cutoff <- 0.05

markers$change <-
  as.factor(ifelse(
    markers$p_val < pvalue.cutoff &
      abs(markers$avg_log2FC) > logFC.cutoff,
    ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
    'NOT'
  ))
table(markers$change)

Volcano_data <- na.omit(markers)
Volcano_data$logP <- -log10(Volcano_data$p_val)

# 输入关注的基因
genes <- c("APOE")#####

Volcano_data$label <- ""
Volcano_data$label[Volcano_data$gene %in% genes] <-
  Volcano_data[Volcano_data$gene %in% genes, ]$gene

Volcano_paired <-
  rasterise(
    ggscatter(
      Volcano_data,
      x = "avg_log2FC",
      y = "logP",
      color = "change",
      palette = c("#4DBBD5FF", "#BBBBBB", "#E64B35FF"),#######
      size = 0.5,
      font.label = 8,
      repel = T,
      xlab = "log2 FoldChange",
      ylab = "-log10 (pvalue)"
    ),
    dpi = 600
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff),
             linetype = "dashed") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    plot.title = element_text(
      size = 8,
      color = "black",
      hjust = 0.5
    ),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank()
  ) +
  geom_text_repel(
    data = Volcano_data,
    aes(x = avg_log2FC,
        y = logP,
        label = label),
    size = 3,
    fontface = "italic",
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE,
    max.overlaps = 100
  ) +ggtitle("ASD VS CON") + theme(plot.title = element_text(hjust = 0.5, size = 20, color = "black", face = "bold")) +  scale_x_continuous(limits = c(-max(abs(
    Volcano_data$avg_log2FC
  )),
  max(abs(
    Volcano_data$avg_log2FC
  ))))


Volcano_paired

#### 保存结果 ----
# 保存图片
ggsave(
  "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/ASDVSCON_volcano.pdf",
  plot = Volcano_paired,
  height = 3,
  width = 3,
  dpi = 600
)

# 保存上调基因
write.table(
  Volcano_data[Volcano_data$change=="UP",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/ASDVSCON_up_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

# 保存下调基因
write.table(
  Volcano_data[Volcano_data$change=="DOWN",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/ASDVSCON_down_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

####——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
setwd("/home/chenzh/brain/SCENIC/by_group/downsrtream_analysis/VS")

CON_TFS<- read.table("~/brain/SCENIC/by_group/CON/tfs_targer.tsv", header=T, sep=",")
ASD_TFS<- read.table("~/brain/SCENIC/by_group/ASD/tfs_targer.tsv", header=T, sep=",")
BD_TFS<- read.table("~/brain/SCENIC/by_group/BD/tfs_targer.tsv", header=T, sep=",")
MDD_TFS<- read.table("~/brain/SCENIC/by_group/MDD/tfs_targer.tsv", header=T, sep=",")
PTSD_TFS<- read.table("~/brain/SCENIC/by_group/PTSD/tfs_targer.tsv", header=T, sep=",")
SZ_TFS<- read.table("~/brain/SCENIC/by_group/SZ/tfs_targer.tsv", header=T, sep=",")

tf_to_keep0 <- CON_TFS$tf
tf_to_keep1 <- ASD_TFS$tf
tf_to_keep2 <- BD_TFS$tf
tf_to_keep3 <- MDD_TFS$tf
tf_to_keep4 <- PTSD_TFS$tf
tf_to_keep5 <- SZ_TFS$tf

# 取并集
union_tfs01 <- union(tf_to_keep0, tf_to_keep1)
union_tfs02 <- union(tf_to_keep0, tf_to_keep2)
union_tfs03 <- union(tf_to_keep0, tf_to_keep3)
union_tfs04 <- union(tf_to_keep0, tf_to_keep4)
union_tfs05 <- union(tf_to_keep0, tf_to_keep5)

scRNA0<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_CON.rds")#读取rds文件
scRNA2<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_BD.rds")#读取rds文件
scRNA_filtered0 <- subset(scRNA0, features = union_tfs02)
scRNA_filtered2 <- subset(scRNA2, features = union_tfs02)

# 加载Seurat包
library(Seurat)

# 假设scRNA_filtered0和scRNA_filtered1都是Seurat对象
# 使用Seurat的`merge`函数来合并这两个对象

combined_seurat <- merge(scRNA_filtered0, scRNA_filtered2)

# 查看合并后的对象
print(combined_seurat)

saveRDS(combined_seurat,"BDVSCON.rds")

seurat.obj<-combined_seurat

table(seurat.obj$group)#查看疾病与正常分组(IVL,Normal)
table(seurat.obj$cell_type)#查看细胞类型分组()

#### 寻找差异基因 ----
# 将细胞类型及刺激状态作为分组
Idents(seurat.obj) <- "group"

# 查看当前分组
table(Idents(seurat.obj))

#使用FindMarkers函数寻找差异表达基因
DefaultAssay(seurat.obj) <- "RNA"

##_____________________________________________________________________________________________

markers <- FindMarkers(seurat.obj,
                       ident.1 = "BD",
                       ident.2 = "CON",
                       logfc.threshold = 0)



markers$celltype <- "BDVSCON"
markers$gene <- rownames(markers)
#### 保存结果 ----
# 差异基因
write.table(
  markers,
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/BDVSCON_deg.txt",
  quote = F,
  sep = ",",
  row.names = F
)

##————————————————————————火山图——————————————————————————————————————————————————————————————————————————————————————————————
#### 导入包 ----
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(dplyr)

#### 导入数据 ----
##markers <- read.table("~/brain/差异/output/deg/ASDVSCON_deg.txt", sep = ",", header = T)

#### 定义差异基因 ----
# 设置阈值
logFC.cutoff <- 0.25
pvalue.cutoff <- 0.05

markers$change <-
  as.factor(ifelse(
    markers$p_val < pvalue.cutoff &
      abs(markers$avg_log2FC) > logFC.cutoff,
    ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
    'NOT'
  ))
table(markers$change)

Volcano_data <- na.omit(markers)
Volcano_data$logP <- -log10(Volcano_data$p_val)

# 输入关注的基因
genes <- c("APOE")#####

Volcano_data$label <- ""
Volcano_data$label[Volcano_data$gene %in% genes] <-
  Volcano_data[Volcano_data$gene %in% genes, ]$gene

Volcano_paired <-
  rasterise(
    ggscatter(
      Volcano_data,
      x = "avg_log2FC",
      y = "logP",
      color = "change",
      palette = c("#4DBBD5FF", "#BBBBBB", "#E64B35FF"),#######
      size = 0.5,
      font.label = 8,
      repel = T,
      xlab = "log2 FoldChange",
      ylab = "-log10 (pvalue)"
    ),
    dpi = 600
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff),
             linetype = "dashed") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    plot.title = element_text(
      size = 8,
      color = "black",
      hjust = 0.5
    ),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank()
  ) +
  geom_text_repel(
    data = Volcano_data,
    aes(x = avg_log2FC,
        y = logP,
        label = label),
    size = 3,
    fontface = "italic",
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE,
    max.overlaps = 100
  ) +ggtitle("ASD VS CON") + theme(plot.title = element_text(hjust = 0.5, size = 20, color = "black", face = "bold")) +  scale_x_continuous(limits = c(-max(abs(
    Volcano_data$avg_log2FC
  )),
  max(abs(
    Volcano_data$avg_log2FC
  ))))


Volcano_paired

#### 保存结果 ----
# 保存图片
ggsave(
  "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/BDVSCON_volcano.pdf",
  plot = Volcano_paired,
  height = 3,
  width = 3,
  dpi = 600
)

# 保存上调基因
write.table(
  Volcano_data[Volcano_data$change=="UP",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/BDVSCON_up_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

# 保存下调基因
write.table(
  Volcano_data[Volcano_data$change=="DOWN",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/BDVSCON_down_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

####——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
####——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
setwd("/home/chenzh/brain/SCENIC/by_group/downsrtream_analysis/VS")

CON_TFS<- read.table("~/brain/SCENIC/by_group/CON/tfs_targer.tsv", header=T, sep=",")
ASD_TFS<- read.table("~/brain/SCENIC/by_group/ASD/tfs_targer.tsv", header=T, sep=",")
BD_TFS<- read.table("~/brain/SCENIC/by_group/BD/tfs_targer.tsv", header=T, sep=",")
MDD_TFS<- read.table("~/brain/SCENIC/by_group/MDD/tfs_targer.tsv", header=T, sep=",")
PTSD_TFS<- read.table("~/brain/SCENIC/by_group/PTSD/tfs_targer.tsv", header=T, sep=",")
SZ_TFS<- read.table("~/brain/SCENIC/by_group/SZ/tfs_targer.tsv", header=T, sep=",")

tf_to_keep0 <- CON_TFS$tf
tf_to_keep1 <- ASD_TFS$tf
tf_to_keep2 <- BD_TFS$tf
tf_to_keep3 <- MDD_TFS$tf
tf_to_keep4 <- PTSD_TFS$tf
tf_to_keep5 <- SZ_TFS$tf

# 取并集
union_tfs01 <- union(tf_to_keep0, tf_to_keep1)
union_tfs02 <- union(tf_to_keep0, tf_to_keep2)
union_tfs03 <- union(tf_to_keep0, tf_to_keep3)
union_tfs04 <- union(tf_to_keep0, tf_to_keep4)
union_tfs05 <- union(tf_to_keep0, tf_to_keep5)

scRNA0<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_CON.rds")#读取rds文件
scRNA3<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_MDD.rds")#读取rds文件
scRNA_filtered0 <- subset(scRNA0, features = union_tfs03)
scRNA_filtered3 <- subset(scRNA3, features = union_tfs03)

# 加载Seurat包
library(Seurat)

# 假设scRNA_filtered0和scRNA_filtered1都是Seurat对象
# 使用Seurat的`merge`函数来合并这两个对象

combined_seurat <- merge(scRNA_filtered0, scRNA_filtered3)

# 查看合并后的对象
print(combined_seurat)

saveRDS(combined_seurat,"MDDVSCON.rds")

seurat.obj<-combined_seurat

table(seurat.obj$group)#查看疾病与正常分组(IVL,Normal)
table(seurat.obj$cell_type)#查看细胞类型分组()

#### 寻找差异基因 ----
# 将细胞类型及刺激状态作为分组
Idents(seurat.obj) <- "group"

# 查看当前分组
table(Idents(seurat.obj))

#使用FindMarkers函数寻找差异表达基因
DefaultAssay(seurat.obj) <- "RNA"

##_____________________________________________________________________________________________

markers <- FindMarkers(seurat.obj,
                       ident.1 = "MDD",
                       ident.2 = "CON",
                       logfc.threshold = 0)



markers$celltype <- "MDDVSCON"
markers$gene <- rownames(markers)
#### 保存结果 ----
# 差异基因
write.table(
  markers,
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/MDDVSCON_deg.txt",
  quote = F,
  sep = ",",
  row.names = F
)

##————————————————————————火山图——————————————————————————————————————————————————————————————————————————————————————————————
#### 导入包 ----
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(dplyr)

#### 导入数据 ----
##markers <- read.table("~/brain/差异/output/deg/ASDVSCON_deg.txt", sep = ",", header = T)

#### 定义差异基因 ----
# 设置阈值
logFC.cutoff <- 0.25
pvalue.cutoff <- 0.05

markers$change <-
  as.factor(ifelse(
    markers$p_val < pvalue.cutoff &
      abs(markers$avg_log2FC) > logFC.cutoff,
    ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
    'NOT'
  ))
table(markers$change)

Volcano_data <- na.omit(markers)
Volcano_data$logP <- -log10(Volcano_data$p_val)

# 输入关注的基因
genes <- c("APOE")#####

Volcano_data$label <- ""
Volcano_data$label[Volcano_data$gene %in% genes] <-
  Volcano_data[Volcano_data$gene %in% genes, ]$gene

Volcano_paired <-
  rasterise(
    ggscatter(
      Volcano_data,
      x = "avg_log2FC",
      y = "logP",
      color = "change",
      palette = c("#4DBBD5FF", "#BBBBBB", "#E64B35FF"),#######
      size = 0.5,
      font.label = 8,
      repel = T,
      xlab = "log2 FoldChange",
      ylab = "-log10 (pvalue)"
    ),
    dpi = 600
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff),
             linetype = "dashed") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    plot.title = element_text(
      size = 8,
      color = "black",
      hjust = 0.5
    ),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank()
  ) +
  geom_text_repel(
    data = Volcano_data,
    aes(x = avg_log2FC,
        y = logP,
        label = label),
    size = 3,
    fontface = "italic",
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE,
    max.overlaps = 100
  ) +ggtitle("ASD VS CON") + theme(plot.title = element_text(hjust = 0.5, size = 20, color = "black", face = "bold")) +  scale_x_continuous(limits = c(-max(abs(
    Volcano_data$avg_log2FC
  )),
  max(abs(
    Volcano_data$avg_log2FC
  ))))


Volcano_paired

#### 保存结果 ----
# 保存图片
ggsave(
  "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/MDDVSCON_volcano.pdf",
  plot = Volcano_paired,
  height = 3,
  width = 3,
  dpi = 600
)

# 保存上调基因
write.table(
  Volcano_data[Volcano_data$change=="UP",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/MDDVSCON_up_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

# 保存下调基因
write.table(
  Volcano_data[Volcano_data$change=="DOWN",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/MDDVSCON_down_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

####——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
####——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
setwd("/home/chenzh/brain/SCENIC/by_group/downsrtream_analysis/VS")

CON_TFS<- read.table("~/brain/SCENIC/by_group/CON/tfs_targer.tsv", header=T, sep=",")
ASD_TFS<- read.table("~/brain/SCENIC/by_group/ASD/tfs_targer.tsv", header=T, sep=",")
BD_TFS<- read.table("~/brain/SCENIC/by_group/BD/tfs_targer.tsv", header=T, sep=",")
MDD_TFS<- read.table("~/brain/SCENIC/by_group/MDD/tfs_targer.tsv", header=T, sep=",")
PTSD_TFS<- read.table("~/brain/SCENIC/by_group/PTSD/tfs_targer.tsv", header=T, sep=",")
SZ_TFS<- read.table("~/brain/SCENIC/by_group/SZ/tfs_targer.tsv", header=T, sep=",")

tf_to_keep0 <- CON_TFS$tf
tf_to_keep1 <- ASD_TFS$tf
tf_to_keep2 <- BD_TFS$tf
tf_to_keep3 <- MDD_TFS$tf
tf_to_keep4 <- PTSD_TFS$tf
tf_to_keep5 <- SZ_TFS$tf

# 取并集
union_tfs01 <- union(tf_to_keep0, tf_to_keep1)
union_tfs02 <- union(tf_to_keep0, tf_to_keep2)
union_tfs03 <- union(tf_to_keep0, tf_to_keep3)
union_tfs04 <- union(tf_to_keep0, tf_to_keep4)
union_tfs05 <- union(tf_to_keep0, tf_to_keep5)

scRNA0<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_CON.rds")#读取rds文件
scRNA4<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_PTSD.rds")#读取rds文件
scRNA_filtered0 <- subset(scRNA0, features = union_tfs04)
scRNA_filtered4 <- subset(scRNA4, features = union_tfs04)

# 加载Seurat包
library(Seurat)

# 假设scRNA_filtered0和scRNA_filtered1都是Seurat对象
# 使用Seurat的`merge`函数来合并这两个对象

combined_seurat <- merge(scRNA_filtered0, scRNA_filtered4)

# 查看合并后的对象
print(combined_seurat)

saveRDS(combined_seurat,"PTSDVSCON.rds")

seurat.obj<-combined_seurat

table(seurat.obj$group)#查看疾病与正常分组(IVL,Normal)
table(seurat.obj$cell_type)#查看细胞类型分组()

#### 寻找差异基因 ----
# 将细胞类型及刺激状态作为分组
Idents(seurat.obj) <- "group"

# 查看当前分组
table(Idents(seurat.obj))

#使用FindMarkers函数寻找差异表达基因
DefaultAssay(seurat.obj) <- "RNA"

##_____________________________________________________________________________________________

markers <- FindMarkers(seurat.obj,
                       ident.1 = "PTSD",
                       ident.2 = "CON",
                       logfc.threshold = 0)



markers$celltype <- "PTSDVSCON"
markers$gene <- rownames(markers)
#### 保存结果 ----
# 差异基因
write.table(
  markers,
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/PTSDVSCON_deg.txt",
  quote = F,
  sep = ",",
  row.names = F
)

##————————————————————————火山图——————————————————————————————————————————————————————————————————————————————————————————————
#### 导入包 ----
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(dplyr)

#### 导入数据 ----
##markers <- read.table("~/brain/差异/output/deg/ASDVSCON_deg.txt", sep = ",", header = T)

#### 定义差异基因 ----
# 设置阈值
logFC.cutoff <- 0.25
pvalue.cutoff <- 0.05

markers$change <-
  as.factor(ifelse(
    markers$p_val < pvalue.cutoff &
      abs(markers$avg_log2FC) > logFC.cutoff,
    ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
    'NOT'
  ))
table(markers$change)

Volcano_data <- na.omit(markers)
Volcano_data$logP <- -log10(Volcano_data$p_val)

# 输入关注的基因
genes <- c("APOE")#####

Volcano_data$label <- ""
Volcano_data$label[Volcano_data$gene %in% genes] <-
  Volcano_data[Volcano_data$gene %in% genes, ]$gene

Volcano_paired <-
  rasterise(
    ggscatter(
      Volcano_data,
      x = "avg_log2FC",
      y = "logP",
      color = "change",
      palette = c("#4DBBD5FF", "#BBBBBB", "#E64B35FF"),#######
      size = 0.5,
      font.label = 8,
      repel = T,
      xlab = "log2 FoldChange",
      ylab = "-log10 (pvalue)"
    ),
    dpi = 600
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff),
             linetype = "dashed") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    plot.title = element_text(
      size = 8,
      color = "black",
      hjust = 0.5
    ),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank()
  ) +
  geom_text_repel(
    data = Volcano_data,
    aes(x = avg_log2FC,
        y = logP,
        label = label),
    size = 3,
    fontface = "italic",
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE,
    max.overlaps = 100
  ) +ggtitle("ASD VS CON") + theme(plot.title = element_text(hjust = 0.5, size = 20, color = "black", face = "bold")) +  scale_x_continuous(limits = c(-max(abs(
    Volcano_data$avg_log2FC
  )),
  max(abs(
    Volcano_data$avg_log2FC
  ))))


Volcano_paired

#### 保存结果 ----
# 保存图片
ggsave(
  "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/PTSDVSCON_volcano.pdf",
  plot = Volcano_paired,
  height = 3,
  width = 3,
  dpi = 600
)

# 保存上调基因
write.table(
  Volcano_data[Volcano_data$change=="UP",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/PTSDVSCON_up_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

# 保存下调基因
write.table(
  Volcano_data[Volcano_data$change=="DOWN",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/PTSDVSCON_down_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

####——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
####——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
setwd("/home/chenzh/brain/SCENIC/by_group/downsrtream_analysis/VS")

CON_TFS<- read.table("~/brain/SCENIC/by_group/CON/tfs_targer.tsv", header=T, sep=",")
ASD_TFS<- read.table("~/brain/SCENIC/by_group/ASD/tfs_targer.tsv", header=T, sep=",")
BD_TFS<- read.table("~/brain/SCENIC/by_group/BD/tfs_targer.tsv", header=T, sep=",")
MDD_TFS<- read.table("~/brain/SCENIC/by_group/MDD/tfs_targer.tsv", header=T, sep=",")
PTSD_TFS<- read.table("~/brain/SCENIC/by_group/PTSD/tfs_targer.tsv", header=T, sep=",")
SZ_TFS<- read.table("~/brain/SCENIC/by_group/SZ/tfs_targer.tsv", header=T, sep=",")

tf_to_keep0 <- CON_TFS$tf
tf_to_keep1 <- ASD_TFS$tf
tf_to_keep2 <- BD_TFS$tf
tf_to_keep3 <- MDD_TFS$tf
tf_to_keep4 <- PTSD_TFS$tf
tf_to_keep5 <- SZ_TFS$tf

# 取并集
union_tfs01 <- union(tf_to_keep0, tf_to_keep1)
union_tfs02 <- union(tf_to_keep0, tf_to_keep2)
union_tfs03 <- union(tf_to_keep0, tf_to_keep3)
union_tfs04 <- union(tf_to_keep0, tf_to_keep4)
union_tfs05 <- union(tf_to_keep0, tf_to_keep5)

scRNA0<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_CON.rds")#读取rds文件
scRNA5<-readRDS("~/brain/SCENIC/by_group/sub/scRNAsub_SZ.rds")#读取rds文件
scRNA_filtered0 <- subset(scRNA0, features = union_tfs05)
scRNA_filtered5 <- subset(scRNA5, features = union_tfs05)

## 加载Seurat包
library(Seurat)

# 假设scRNA_filtered0和scRNA_filtered1都是Seurat对象
# 使用Seurat的`merge`函数来合并这两个对象

combined_seurat <- merge(scRNA_filtered0, scRNA_filtered5)

# 查看合并后的对象
print(combined_seurat)

saveRDS(combined_seurat,"SZVSCON.rds")

seurat.obj<-combined_seurat

table(seurat.obj$group)#查看疾病与正常分组(IVL,Normal)
table(seurat.obj$cell_type)#查看细胞类型分组()

#### 寻找差异基因 ----
# 将细胞类型及刺激状态作为分组
Idents(seurat.obj) <- "group"

# 查看当前分组
table(Idents(seurat.obj))

#使用FindMarkers函数寻找差异表达基因
DefaultAssay(seurat.obj) <- "RNA"

##_____________________________________________________________________________________________

markers <- FindMarkers(seurat.obj,
                       ident.1 = "SZ",
                       ident.2 = "CON",
                       logfc.threshold = 0)



markers$celltype <- "SZVSCON"
markers$gene <- rownames(markers)
#### 保存结果 ----
# 差异基因
write.table(
  markers,
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/SZVSCON_deg.txt",
  quote = F,
  sep = ",",
  row.names = F
)

##————————————————————————火山图——————————————————————————————————————————————————————————————————————————————————————————————
#### 导入包 ----
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(dplyr)

#### 导入数据 ----
##markers <- read.table("~/brain/差异/output/deg/ASDVSCON_deg.txt", sep = ",", header = T)

#### 定义差异基因 ----
# 设置阈值
logFC.cutoff <- 0.25
pvalue.cutoff <- 0.05

markers$change <-
  as.factor(ifelse(
    markers$p_val < pvalue.cutoff &
      abs(markers$avg_log2FC) > logFC.cutoff,
    ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
    'NOT'
  ))
table(markers$change)

Volcano_data <- na.omit(markers)
Volcano_data$logP <- -log10(Volcano_data$p_val)

# 输入关注的基因
genes <- c("APOE")#####

Volcano_data$label <- ""
Volcano_data$label[Volcano_data$gene %in% genes] <-
  Volcano_data[Volcano_data$gene %in% genes, ]$gene

Volcano_paired <-
  rasterise(
    ggscatter(
      Volcano_data,
      x = "avg_log2FC",
      y = "logP",
      color = "change",
      palette = c("#4DBBD5FF", "#BBBBBB", "#E64B35FF"),#######
      size = 0.5,
      font.label = 8,
      repel = T,
      xlab = "log2 FoldChange",
      ylab = "-log10 (pvalue)"
    ),
    dpi = 600
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff),
             linetype = "dashed") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    plot.title = element_text(
      size = 8,
      color = "black",
      hjust = 0.5
    ),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank()
  ) +
  geom_text_repel(
    data = Volcano_data,
    aes(x = avg_log2FC,
        y = logP,
        label = label),
    size = 3,
    fontface = "italic",
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE,
    max.overlaps = 100
  ) +ggtitle("ASD VS CON") + theme(plot.title = element_text(hjust = 0.5, size = 20, color = "black", face = "bold")) +  scale_x_continuous(limits = c(-max(abs(
    Volcano_data$avg_log2FC
  )),
  max(abs(
    Volcano_data$avg_log2FC
  ))))


Volcano_paired

#### 保存结果 ----
# 保存图片
ggsave(
  "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/SZVSCON_volcano.pdf",
  plot = Volcano_paired,
  height = 3,
  width = 3,
  dpi = 600
)

# 保存上调基因
write.table(
  Volcano_data[Volcano_data$change=="UP",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/SZVSCON_up_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

# 保存下调基因
write.table(
  Volcano_data[Volcano_data$change=="DOWN",],
  file = "~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/SZVSCON_down_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)
