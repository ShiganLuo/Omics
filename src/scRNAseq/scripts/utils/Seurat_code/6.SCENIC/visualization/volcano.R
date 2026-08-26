#### 导入包 ----
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)

#### 导入数据 ----
markers <- read.table("./results/deg/tc.txt", sep = ",", header = T)

#### 定义差异基因 ----
# 设置阈值
logFC.cutoff <- 0.5
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
genes <- c("TM4SF1", "COL6A3", "APOE", "DDX5")
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
      palette = c("#4DBBD5FF", "#BBBBBB", "#E64B35FF"),
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
  ) +
  scale_x_continuous(limits = c(-max(abs(
    Volcano_data$avg_log2FC
  )),
  max(abs(
    Volcano_data$avg_log2FC
  ))))
Volcano_paired

#### 保存结果 ----
# 保存图片
ggsave(
  "./results/deg/tc_volcano.pdf",
  plot = Volcano_paired,
  height = 3,
  width = 3,
  dpi = 600
)

# 保存上调基因
write.table(
  Volcano_data[Volcano_data$change=="UP",],
  file = "./results/deg/tc_up_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

# 保存下调基因
write.table(
  Volcano_data[Volcano_data$change=="DOWN",],
  file = "./results/deg/tc_down_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)
