###mitoGSEA计算-----
####27mo与2mo的GSEA棒棒糖图-----
library(msigdbr)
library(tidyverse)
library(readxl)
set.seed(717)
setwd("~/精神疾病/data/388/Seurat/data/GSEA/")
output <- "~/精神疾病/data/388/Seurat/data/GSEA/CSV/"
tissues <- c("ASD", "BD", "MDD", "PTSD", "SZ")  
library(tidyverse)  

adata <- NULL  

for (i in tissues) {  
  file_path <- paste0("~/精神疾病/data/388/Seurat/data/GSEA/CSV/gsea_", i, "-table.csv")  
  
  data <- read.csv(file_path)  
  
  if (!all(c("pvalue", "NES", "ID") %in% names(data))) {  
    stop("One or more of the required columns (pvalue, NES, ID) is missing in the CSV file.")  
  }  
  
  data$a <- -log10(data$pvalue)  
  data$tissue <- i  
  data$change <- ifelse(data$NES > 0, "up", "down")  
  
  if (is.null(adata)) {  
    adata <- data  
  } else {  
    adata <- rbind(adata, data)  
  }  
}  
write.csv(adata, file = "~/精神疾病/data/388/Seurat/data/GSEA/CSV/allgsea_mental_diseases.csv", row.names = FALSE)


up_data <- subset(adata,subset = adata$change %in% "up")
a <- table(up_data$ID)
a <- as.data.frame(a)
a1 <- subset(a,subset = a$Freq > 3)
up_data1 <- subset(up_data,subset = ID %in% a1$Var1)
library(forcats)
up_data1$Description <- as.factor(up_data1$Description)
up_data1$Description <- fct_inorder(up_data1$Description)
up_data1 <- subset(up_data1,subset = pvalue < 0.05)

p1 <- ggplot(up_data1, aes(tissue, Description)) +
  geom_point(aes(color=NES, size=-log10(pvalue)))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black")
  )+
  scale_color_gradient(low = "#F4A582", high="#D6604D")+
  labs(x = NULL, y = NULL)+
  guides(size = guide_legend(order = 1))

p1
ggsave(paste0(output,"4个及以上器官线粒体功能上调气泡图.pdf"), plot = p1, width = 6.75, height = 6.69)
a1 <- subset(a,subset = a$Freq %in% 1)
up_data1 <- subset(up_data,subset = ID %in% a1$Var1)
library(forcats)
up_data1$Description <- as.factor(up_data1$Description)
up_data1$Description <- fct_inorder(up_data1$Description)
up_data1 <- subset(up_data1,subset = pvalue < 0.05)

p2 <- ggplot(up_data1, aes(tissue, Description)) +
  geom_point(aes(color=NES, size=-log10(pvalue)))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black")
  )+
  scale_color_gradient(low = "#F4A582", high="#D6604D")+
  labs(x = NULL, y = NULL)+
  guides(size = guide_legend(order = 1))

p2
ggsave(paste0(output,"特异性下调气泡图.pdf"), plot = p2, width = 6.75, height = 10)
data <- read.csv(paste0("/home/lixy/SP-hope/results/deg/files/",i,"/48_vs_0-deg.csv"))
data$tissue <- i
adata <- data.frame(matrix(nrow = 1, ncol = 10))
colnames(adata) <- colnames(data)
tissues <- c(tissues <- c("Bl","Br","C", "H","K","Li","LN","Lu","P","SI","Sk","Sp","Th"))
library(tidyverse)
i<-"Bl"
for (i in tissues) {
  deg <- read.csv(paste0("/home/lixy/SP-hope/results/deg/files/",i,"/48_vs_0-deg.csv"))
  deg$tissue <- i
  dbs <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/Mouse.MitoCarta3.0.xls")
  gene_df <- dbs %>%
    mutate(gene = strsplit(Genes, ", ")) %>%
    unnest(gene)
  deg <- subset(deg,subset = deg$change %in% "Upregulated")
  deg <- subset(deg,subset = gene %in% gene_df$gene)
  adata <- rbind(adata,deg)
}
adata <- adata[-1,]
library(gdata)
library(reshape2)
library(UpSetR)
up_wide<-dcast(adata, gene~adata$tissue,
                 value.var = 'gene')
up_wide <- up_wide[,-1]
output <- "~/RNAseq/results/gsea/72_0/"
sets <- c(tissues <- c("Bl","Br","C", "H","K","Li","LN","Lu","P","SI","Sk","Sp","Th"))
pdf(file.path(paste0(output,"上调线粒体基因upset图.pdf")), width = 13.64, height = 9.54)

p1 <- upset(fromList(up_wide), 
            order.by = "freq",
            nsets = 100,
            nintersects = 30,
            mb.ratio = c(0.7, 0.3),#柱状图与矩阵点图之间的比例大小
            show.numbers = "yes",
            point.size = 1.8,#矩阵中圆圈的大小
            line.size = 0.8,
            sets.x.label = "Set Size", #柱状图的轴标签
            main.bar.color = "#CD5C5C",
            sets.bar.color = c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#F0027F","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#6495ED",
                               "cyan1", "royalblue4","darksalmon"),
            matrix.color = "#CD5C5C",
            mainbar.y.label = "Gene Intersections", 
            text.scale = 1.75,
            shade.color = "#CD5C5C",
            shade.alpha = 0.1)
p1
dev.off()
setwd("/home/caojt/RNAseq/")
output <- "./figure/results/gsea/"
dir.create(output, recursive = TRUE)
dir.create(output, recursive = T)
counts <-
  read.csv(
    file = "/home/gongfengcz/RNAseq/data/多器官/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv",
    row.names = 1
  )
library(DESeq2)
library(readr)
library(tidyverse)
library(BiocParallel)
metadata <- read.csv("/home/gongfengcz/RNAseq/data/多器官/GSE132040_MACA_Bulk_metadata.csv")
tissues <- c("Bl","Br","C", "H","K","Li","LN","Lu","P","SI","Sk","Sp","Th")
ages <- c("6","12","24","48","72","96")
i <- "Bl"
for (i in tissues) {
  output1 <- paste0("./figure/lanscape/",i)
  dir.create(output1, recursive = TRUE)
  metadata1 <- metadata[grep(i, metadata$source.name), ]
  colnames(counts) <- gsub(".gencode.vM19", "", colnames(counts))
  for (n in ages) {
    positive_group <- n
    negative_group <- "3"
    
    filtered_metadata <- metadata1 %>%
      filter(characteristics..age %in% c(positive_group, negative_group)) %>%
      mutate(group = as_factor(characteristics..age))
    
    filtered_count_data <- counts %>%
      dplyr::select(pull(filtered_metadata, Sample.name))
    register(MulticoreParam(10))
    if (all(pull(filtered_metadata, Sample.name) %in% colnames(filtered_count_data))) {
      dds <- DESeqDataSetFromMatrix(countData = filtered_count_data,
                                    colData = filtered_metadata,
                                    design = ~ group)
      dds
      dds <- DESeq(dds, parallel = T)
      saveRDS(dds, file.path(
        output,i,
        str_c(positive_group, "_vs_", negative_group, "-dds.rds")
      ))
      res <-
        results(dds, contrast = c("group", positive_group, negative_group))
      res_ordered <- res[order(res$padj),]
      deg <- na.omit(as.data.frame(res_ordered))
      deg <- deg %>% mutate(change = ifelse(
        pvalue <= 0.05 & abs(log2FoldChange) >= 1,
        ifelse(log2FoldChange > 1, "Upregulated", "Downregulated"),
        "Stable"
      )) %>%
        rownames_to_column("gene")
      deg %>% dplyr::count(change)
      
      write.csv(deg, file.path(
        output,i,
        str_c(positive_group, "_vs_", negative_group, "-deg.csv")
      ))
    }}}
setwd("~/RNAseq/")
library(readr)
library(readxl)
library(UCell)
library(Seurat)
library(ggplot2)
vision <- read.csv("~/RNAseq/results/vision_mito.csv")
metadata <- read.csv("/home/lixy/SP-hope/data/expr-merged/20220725_metadata.csv")
split_vector <- strsplit(metadata$group, "-") # 使用下划线分割字符串

first_part <- sapply(split_vector, "[", 1)

metadata$group <- first_part
metadata$group <- paste0(metadata$group,"-",metadata$age)
vision1 <- merge(metadata, vision, by.x = "sample", by.y = "X", all = F)
vision1 <- vision1[,-c(2:4)] 
colnames(vision1)
vision1 <- aggregate(. ~ sample, data = vision1, mean)
vision1_long<-melt(vision1,
                   id.vars = c('sample'),#需要保留不参与聚合的变量,
                   variable.name='Pathways',
                   value.name='Score')
vision1_long$tissue <- sapply(strsplit(vision1_long$sample, split = "-"), "[", 1)
vision1_long$time <- sapply(strsplit(vision1_long$sample, split = "-"), "[", 2)


vision2_long <- subset(vision1_long,subset = vision1_long$Pathways %in% "OXPHOS")
vision2_long <- vision2_long[!grepl("NA", vision2_long$sample), ]
vision2_long <- subset(vision2_long,subset = vision2_long$tissue %in% c("C","K","Li","LN","SI","Sk"))

vision2_long$time <- factor(vision2_long$time, levels = c('0', '6', '12', '24',"48","72","96"))
vision2_long <- vision2_long[order(vision2_long$time), ]

ggplot(vision2_long, aes(vision2_long$time, vision2_long$Score, group=tissue, color=tissue, shape=tissue))+
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#FDBE85", "#D94701", "#8856A7", "#9EBCDA", "#FAA75B", "#882E72", "#5289C7", "#5F4690", "#F04590", "#C63F6D", "#913C3E", "#566E3D", "#1D9183", "#2F65A5", "#243475")) + #自定义颜色
  theme_bw()+
  labs(x = "Time", # 定义x轴文本
       y = "OXPHOS")# 定义y轴文本

vision2_long <- subset(vision1_long,subset = vision1_long$Pathways %in% "OXPHOS")
vision2_long <- vision2_long[!grepl("NA", vision2_long$sample), ]
vision2_long <- subset(vision2_long,subset = vision2_long$tissue %in% c("P","SI","Sk","Sp"))

vision2_long$time <- factor(vision2_long$time, levels = c("0","6", '12', '24',"48","72","96"))
vision2_long <- vision2_long[order(vision2_long$time), ]

ggplot(vision2_long, aes(vision2_long$time, vision2_long$Score, group=tissue, color=tissue, shape=tissue))+
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_color_manual(values = c("grey", "grey", "grey", "firebrick", "#882E72", "#5289C7", "#5F4690", "#F04590", "#C63F6D", "#913C3E", "#566E3D", "#1D9183", "#2F65A5", "#243475")) + #自定义颜色
  theme_bw()+
  labs(x = "Time", # 定义x轴文本
       y = "Mitophagy")# 定义y轴文本

