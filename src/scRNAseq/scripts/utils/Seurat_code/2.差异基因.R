######### 导入包 ----
library(Seurat)

#### 配置 ----
set.seed(717)
seurat.obj<-readRDS("~/脑课题/result/Seurat/ALL_harmony_无MDD.rds")#读取rds文件
table(seurat.obj$group)#查看各个疾病细胞数
# ASD    BD   CON  PTSD    SZ 
# 16326  5524 10536  3523  5983 

#####DEG分析
# 将细胞类型及刺激状态作为分组
Idents(seurat.obj) <- "group"
table(Idents(seurat.obj))
# CON   ASD    BD  PTSD    SZ 
# 10536 16326  5524  3523  5983

#使用FindMarkers函数寻找差异表达基因
DefaultAssay(seurat.obj) <- "RNA"
markers <- FindMarkers(seurat.obj,
                       ident.1 = "ASD",
                       ident.2 = "CON",
                       logfc.threshold = 0)
markers$disease <- "ASD"
markers$gene <- rownames(markers)

#### 保存结果 ----
# 差异基因
write.table(markers, file = "~/脑课题/result/DEG/harmony/ASD_deg1.txt", quote = F, sep = ",", row.names = F)


############火山图（DEG结果可视化）
#### 导入包 ----
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(dplyr)

markers<-read.table("~/脑课题/result/DEG/harmony/SZ_deg.txt", header = T,sep=",")

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



# ------绘制火山图-------
Volcano_paired <- rasterise(
  ggscatter(
    Volcano_data,
    x = "avg_log2FC",
    y = "logP",
    color = "change",
    palette = c("#4DBBD5FF", "#BBBBBB", "#E64B35FF"),  # 设置颜色调色板
    size = 0.5,
    font.label = 8,
    repel = TRUE,
    xlab = "log2 FoldChange",
    ylab = "-log10 (pvalue)"
  ),
  dpi = 600
)
Volcano_paired

#---------------------------------------------------
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
      repel = TRUE,
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
    plot.title = element_text(size = 8, color = "black", hjust = 0.5),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank()
  ) +
  scale_x_continuous(limits = c(-max(abs(Volcano_data$avg_log2FC)), 
                                max(abs(Volcano_data$avg_log2FC))))


#-------- 输入关注的基因-----------
genes <- c("ADIRF","IFI27","COL1A2")

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
  "~/脑课题/result/DEG/SZ_volcano.pdf",
  plot = Volcano_paired,
  height = 6,
  width = 6,
  dpi = 600
)

# 保存上调基因
write.table(
  Volcano_data[Volcano_data$change=="UP",],
  file = "~/脑课题/result/DEG/SZ_up_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)

# 保存下调基因
write.table(
  Volcano_data[Volcano_data$change=="DOWN",],
  file = "~/脑课题/result/DEG/SZ_down_deg.csv",
  quote = F,
  sep = ",",
  row.names = F
)







#### GSEA分析 ----
# gsea分析输入文件准备 ----
gsea.input <- FindMarkers(
  seurat.obj,
  ident.1 = "SZ",
  ident.2 = "CON",
  min.pct = 0,
  logfc.threshold = 0
)
gsea.input$celltype <- "SZ"
gsea.input$gene <- rownames(gsea.input)

# GSEA输入文件保存
write.table(
  gsea.input,
  file = "~/精神疾病/data/388/Seurat/data/GSEA/SZ_gsea.txt",
  quote = F,
  sep = ",",
  row.names = F
)
#saveRDS(seurat.obj,"/home/zhaor/exercise/data/merge.rds")


# 设置基因集
# H, C2, C5, C6, C7
# H: 标志基因集是通过聚集许多MSigDB基因集来表达明确的生物状态或过程的一致表达的标签。
# C2: 来自在线路径数据库、PubMed出版物和领域专家知识的策展基因集。
# C5: 本体基因集由用相同本体术语注释的基因组成。
# C6: 直接从来自癌症基因扰动的微阵列基因表达数据定义的致癌特征基因集。
# C7: 免疫特征基因集代表免疫系统内的细胞状态和扰动。
category <- "C5"
genesets <- msigdbr(species = "Homo sapiens",
                    category = category)
genesets <- subset(genesets,
                   select = c("gs_name", "gene_symbol"))

rownames(gsea.input) <- gsea.input$gene
gsea.input <-
  gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
genelist <-
  structure(gsea.input$avg_log2FC, names = rownames(gsea.input))

write.table(genelist, quote = F, sep = "\t", col.names = F,
            file = "~/精神疾病/data/388/Seurat/data/GSEA/gsea_SZ.rnk")

res <- GSEA(genelist, TERM2GENE = genesets, eps = 0)

#### 保存结果 ----
# 分析结果
saveRDS(res, "~/精神疾病/data/388/Seurat/data/GSEA/gsea_SZ.rds")

# 导出表格
write.csv(res, "~/精神疾病/data/388/Seurat/data/GSEA/gsea_SZ-table.csv",
        row.names = F)

res<-readRDS("~/精神疾病/data/388/Seurat/data/GSEA/gsea_PTSD.rds")

# 导出图形
library(enrichplot)
p <- gseaplot2(res, 13,
               color = "firebrick",
               pvalue_table = T)
p
ggsave(filename = "~/精神疾病/data/388/Seurat/picture/GSEA/PTSD4.pdf",
       p,
       width = 15,
       height = 8)










#### GO分析------
library(clusterProfiler)
library(org.Hs.eg.db)
outdir<-("~/精神疾病/data/388/Seurat/data/GO")


#### 2.1 LM vs N -----
case <- 'MDD_VS_CON'
markers<-read.table("~/精神疾病/data/388/Seurat/data/DEG/1/MDD_deg.txt",
                    header = T,sep=",")

deg_up <- filter(markers,
                 markers$p_val < 0.05 & markers$avg_log2FC > 0.25)
# markers[markers$change %in% 'UP', ]$gene
# head(deg)
# #重新排序
# top.genes <-deg
# deg[order(deg$avg_log2FC, decreasing = T),]$gene
# class(deg)
# deg<-as.data.frame(deg)
bp <-
  enrichGO(
    deg_up$gene,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

term <- bp@result
#保存结果
write.table(term,
            file = file.path(outdir, paste0(case,'_UP_GOBP.csv')),
            quote = F,
            sep = ",",
            row.names = F
)
#绘制TOP20
df1 <- term[1:20,]####提取上调通路结果文件中的前二十条
df1$labelx=rep(0,nrow(df1))
df1$labely=seq(nrow(df1),1)
p <- ggplot(data = df1, 
            aes(x = -log10(pvalue),
                y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity",
           alpha=1,
           fill= "#CD3333",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = df1$Description),
            size=3.5,
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  ggtitle(case)+
  scale_x_continuous(expand = c(0,0))
p
#保存
ggsave(file.path(outdir, paste0(case, '_UP_Barplot.pdf')),#这里也可以重新设置一个outdir2存到不同文件夹中
       p,
       height=7,
       width=5)




#### Down 
deg_down <- filter(markers,
                   markers$p_val < 0.05 & markers$avg_log2FC < -0.25)

bp <-
  enrichGO(
    deg_down$gene,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

term <- bp@result
#保存结果
write.table(term,
            file = file.path(outdir, paste0(case,'_DOWN_GOBP.csv')),
            quote = F,
            sep = ",",
            row.names = F
)
#绘制TOP20
df2 <- term[1:20,]
df2$labelx=rep(0,nrow(df2))
df2$labely=seq(nrow(df2),1)
p <- ggplot(data = df2, 
            aes(x = -log10(pvalue),
                y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity",
           alpha=1,
           fill= "#2166AC",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = df2$Description),
            size=3.5,
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  ggtitle(case)+
  scale_x_continuous(expand = c(0,0))
p
#保存
ggsave(file.path(outdir, paste0(case, '_DOWN_Barplot.pdf')),
       p,
       height=7,
       width=5)


