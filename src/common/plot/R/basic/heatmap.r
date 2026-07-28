heatmap_tag <- function() {
    # 加载必要的包
library(ComplexHeatmap)
library(circlize)
library(readr)

# 读取数据文件
data <- read_csv("data.csv")
group <- read_csv("group.csv", col_names = c("Sample", "Group"))
label_genes <- read_csv("label.csv", col_names = "Gene")

# 数据预处理
data_matrix <- as.matrix(data[, -1])
rownames(data_matrix) <- data$gene
zscore_matrix <- t(scale(t(data_matrix)))

# 创建分组注释
group_colors <- c(control = "#c77cff", patient = "#ff9900")
ha_column <- HeatmapAnnotation(
  Group = group$Group,
  col = list(Group = group_colors),
  annotation_name_side = "right")
# 颜色映射函数
color_fun <- colorRamp2(c(-2, 0, 2), c("#c77cff", "#FFFFFF", "#ff9900"))

# 定义需要高亮的基因行（带指引线）
highlight_idx <- which(rownames(zscore_matrix) %in% label_genes$Gene)
row_anno <- rowAnnotation(
  link = anno_mark(
    at = highlight_idx,
    labels = rownames(zscore_matrix)[highlight_idx],
    labels_gp = gpar(fontsize = 6),
    link_gp = gpar(lwd = 0.5)),
  width = unit(0.2, "cm"))

# 绘制热图主体
ht <- Heatmap(
  zscore_matrix,
  name = "Z-Score",
  col = color_fun,
  top_annotation = ha_column,
  show_row_names = FALSE, # 关闭默认行标签
  show_column_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  right_annotation = row_anno, # 在右侧添加指引线注释
  heatmap_legend_param = list(
    direction = "vertical",
    title_position = "topleft",
    legend_height = unit(3, "cm")
  ))

# 绘制图形并调整图例布局
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE,
    padding = unit(c(5, 5, 5, 10), "mm"),
    gap = unit(5, "mm"))
}
heatmap <- function() {
ann="/home/lsg/GWAS/Betta_splendens240901/haplotype/data/ID/Splendens615_Pop.txt"
input="/home/lsg/GWAS/Betta_splendens240901/haplotype/data/heatmap/sort_haplot.csv"
#加载所需的软件包
library(pheatmap)
#加载注释内容
annotation_row <- read.table(ann,header=T,row.names = 1,stringsAsFactors = T)
head(annotation_row)
#加载haplotype文件
hap <- read.table(input, header=T,row.names = 1, stringsAsFactors = F)
# # # color # 根据自己的调整

group = list(Pop = c("Stiktos" = "#735380","Smaragdina_guitar" = "#660000","Smaragdina" = "#8A2BE2","Siamorientalis" = "#20B2AA",
"Mahachaiensis" = "#FA8072","Imbellis" = "#E1FF00","Splendens" = "#BC80BD","Yellow" = "#E6D214","White" = "#EFDFBB",
"Veiltail" = "#AED0CD","Turquoise" = "#2ED29F","Steel" = "#4682B4","Royal" = "#0F55F1","Red" = "#C63731","Orange" = "#F7953B",
"Koi" = "#DE77AE","Halfmoon" = "#A2B61E" ,"Giant" = "#6F5CCC","Fighter" = "#907C62","Dumbo_HMPK" = "#F6D5A4",
"Dumbo_Halfmoon" = "#C47664","Dragon" = "#177D7A","Crowntail" = "#91FF00","Copper" = "#00FFF7","Black" = "#BEBFC5"))


p <- pheatmap(hap, show_colnames=FALSE, 
                legend_breaks = -1:2, 
                legend_labels = c("./.", "Reference", "Heterozygote", "Altemative"), 
                # color = c("#FFFFFF","#F9EFE5","#F6C8E0","#F472B5"),
                color = c("#FFFFFF","#21E706","#0C13D7","#D7270C"),
                #cutree_cols=2, 
                cluster_row = FALSE, 
                cluster_col = FALSE,
                annotation_row = annotation_row, 
                show_rownames = FALSE,
                annotation_colors = group, 
                # annotation_names_row = TRUE,
                # annotation_legend = FALSE,
                # labels_row = labels_row,
                main="chr3:0-5.21Mb",
                fontsize = 20,
                fointsize= 30,
                )

#pdf(file = "AAA.pdf", width = 10, height = 5)
png(file = paste0("chr3:0-5.21Mb",".png"), width = 6000, height = 4000,res = 300)
p
dev.off()
print("结束")
}
