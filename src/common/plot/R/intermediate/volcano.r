volcano <- function (res, outjpeg,mode = "gene",geneAnnotation = "", nlabel = 10, label.by = "padj"){
  # add gene_name column to the results table
  if ( mode == 'TE'){
    res = res %>% 
    rownames_to_column(var = "index") %>%
    mutate(gene_name = word(index, 1, sep = ":"))
  } else if ( mode == 'gene' && geneAnnotation != ""){
    df <- read.csv(geneAnnotation, sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
    # it is gene mode as long as it contains ensembl_id ,this can resolve that both str1:str2:…… ensembl_id
    res <- res %>%
      rownames_to_column(var = "index") %>%
      mutate(gene_name = ifelse(
      str_count(index, ":") >= 2, 
      str_extract(index, "^[^:]+"),  # else extract gene_name from df
      df$gene_name[match(index, rownames(df))]))
  } else if ( mode == 'gene' && geneAnnotation == ""){
    stop("Invalid geneAnnotation argument. Please provide a gene annotation file, because you choose mode gene.")
  } else {
    stop("Invalid mode argument. Choose either gene or TE.")
  }

    # print(head(res))
  # assign significance to results based on padj
  res <- res %>%
  mutate(
    gene_status = case_when(
      log2FoldChange > 0.58 & -log10(padj) > 1.3 ~ "Upregulated genes",
      log2FoldChange < -0.58 & -log10(padj) > 1.3 ~ "Downregulated genes",
      TRUE ~ "Other genes"
    )
  )
  res = res[!is.na(res$padj),]
  # Calculate the number of each gene_status
  status_counts <- res %>%
    group_by(gene_status) %>%
    summarise(count = n())

  # Create new legend labels containing quantities
  new_labels <- setNames(
    sapply(status_counts$gene_status, function(status) {
      count <- status_counts$count[status_counts$gene_status == status]
      paste0(status, " (n = ", count, ")")
    }),
    status_counts$gene_status
  )
  # print(head(res))
  
  # get labels for the highest or lowest genes according to either padj or log2FoldChange
  if (label.by == "padj") {
    up_genes <- res %>%
      filter(gene_status == "Upregulated genes") %>%
      arrange(padj) %>%
      head(nlabel)
    down_genes <- res %>%
      filter(gene_status == "Downregulated genes") %>%
      arrange(padj) %>%
      head(nlabel)
  } else if (label.by == "log2FoldChange") {
    up_genes <- head(arrange(res[res[,"gene_status"] == "Upregulated genes" ,], desc(log2FoldChange)),nlabel)
    down_genes <- head(arrange(res[res[,"gene_status"] == "Upregulated genes" ,], log2FoldChange),nlabel)
  } else
    stop ("Invalid label.by argument. Choose either padj or log2FoldChange.")
  
    p = ggplot(res, aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(col=gene_status)) + 
      # scale_fill_manual(values=c("Upregulated genes" ="red","Downregulated genes" = "blue","Other genes" = "black"), 
      #                       name = "gene classification", 
      #                       labels = c("Upregulated genes" = uplegendLabel, "Downregulated genes" = downlegendLabel,"Other genes" = otherlegendLabel))+
      scale_color_manual(values=c("Upregulated genes" ="red","Downregulated genes" = "blue","Other genes" = "black"),
        labels = new_labels) + 
      ggrepel::geom_text_repel(data=up_genes, aes(label=head(gene_name,nlabel)), color = "#F59494", size = 3)+
      ggrepel::geom_text_repel(data=down_genes, aes(label=head(gene_name,nlabel)), color = "#93ACF6", size = 3)+
      labs ( x = expression(log[2]("FoldChange")), y = expression(-log[10]("adjusted p-value")))+
      geom_vline(xintercept = 0.58, linetype = "dotted")+
      geom_vline(xintercept = -0.58, linetype = "dotted")+
      geom_hline(yintercept = 1.3, linetype = "dotted")+
      theme_classic() +  
      theme(
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA,linewidth = 0.2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()   
      )
  jpeg(width = 200, height = 150, units = "mm", res = 300, file = outjpeg)
    print(p)
  dev.off()
}

volcano1 <- function() {
  # 加载必要的包
library(ggplot2)
library(ggrepel)
library(dplyr)
# 读取数据
data <- read.csv("patient-control_diff_new.csv")
# 数据处理：计算-log10(P_Value)，添加显著性标记
data <- data %>%
  mutate(
    logP = -log10(P_Value),
    color_group = case_when(
      P_Value >= 0.05 ~ "not_sig",
      log2FC > 1 ~ "up", # 上调（使用阈值1）
      log2FC < -1 ~ "down", # 下调（使用阈值-1）
      TRUE ~ "not_sig"            # 中间区域
    )
  )

# 提取log2FC绝对值最大且P<0.05的10个基因（按绝对值排序）
top_genes <- data %>%
  filter(P_Value < 0.05) %>%
  group_by(GENE_SYMBOL) %>% # 按基因符号去重
  slice_max(abs(log2FC), n = 1) %>% # 取每个基因的最大绝对值
  ungroup() %>%
  arrange(desc(abs(log2FC))) %>%
  head(10) # 取前10个

# 定义新的配色方案
point_colors <- c("up" = "#CB6CB1", "down" = "#13A3A6", "not_sig" = "gray80")
label_fill <- c("up" = "#ffa500", "down" = "#B3A9EB") # 橙色和浅紫色

# 绘制火山图
ggplot(data, aes(x = log2FC, y = logP)) +
  geom_point(aes(color = color_group), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = point_colors) +

# 添加辅助线
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", linewidth = 0.5) +

# 添加基因标签
  geom_label_repel(
    data = top_genes,
    aes(label = GENE_SYMBOL, fill = ifelse(log2FC > 0, "up", "down")),
    color = "white",
    box.padding = 1.2, # 增加标签间距
    segment.color = "grey30", # 连接线颜色
    nudge_x = ifelse(top_genes$log2FC > 0, 3, -3), # 更强制偏移量
    direction = "both", # 允许双向调整
    max.overlaps = 20, # 增加最大重叠容忍度
    size = 3.5                        # 标签文字大小
  ) +

# 颜色填充映射
  scale_fill_manual(values = label_fill) +

# 坐标轴标签
  labs(x = "log2(Fold Change)", y = "-log10(P-Value)") +

# 主题设置
  theme_classic() + # 经典主题
  theme(
    legend.position = "none", # 隐藏图例
    axis.title = element_text(size = 12, face = "bold"), # 加粗坐标轴标题
    axis.text = element_text(size = 10) # 坐标轴文字大小
  ) +
# 坐标轴扩展范围（确保标签可见）
  expand_limits(x = c(max(abs(data$log2FC)) * 1.2 * c(-1, 1)))
}