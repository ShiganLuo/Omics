###mitoGSEA计算-----
####27mo与2mo的GSEA棒棒糖图-----
library(msigdbr)
library(tidyverse)
library(readxl)
set.seed(717)
#### 导入数据 ----
# tissues <- c(tissues <- c("Bl","BM","Br","C", "H","K","Li","LN","Lu","P","SI","Sk","Sp","Th"))
# output <- "/home/RNAseq"
# i<-"Bl"
# for (i in tissues) {
#   deg <- read.csv(paste0("/home/lixy/SP-hope/results/deg/files/",i,"/48_vs_0-deg.csv"))
#   output <- paste0("/home/caojt/RNAseq/figure/lanscape/多器官/",i,"/mito_gsea/")
#   dir.create(output, recursive = T)
#   #### GSEA分析 ----
#   # 设置基因集
#   # H, C2, C5, C6, C7
#   # H: hallmark gene sets  are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.
#   # C2: curated gene sets  from online pathway databases, publications in PubMed, and knnowledge of domain experts.
#   # C5: ontology gene sets  consist of genes annotated by the same ontology term.
#   # C6: oncogenic signature gene sets  defined directly from microarray gene expression data from cancer gene perturbations.
#   # C7: immunologic signature gene sets  represent cell states and perturbations within the immune system.
#   dbs <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/Mouse.MitoCarta3.0.xls")
#   gene_df <- dbs %>%
#     mutate(gene = strsplit(Genes, ", ")) %>%
#     unnest(gene)
#   genesets <- gene_df[,c(1,3)]
#   colnames(genesets) <- c("gs_name", "gene_symbol")
#   deg <- deg %>% column_to_rownames("gene")
#   deg <-
#     deg[order(deg$log2FoldChange, decreasing = T), ]
#   genelist <-
#     structure(deg$log2FoldChange, names = rownames(deg))
#   res <- clusterProfiler::GSEA(genelist,
#                                TERM2GENE = genesets,pvalueCutoff = 1,eps = 1e-1000,minGSSize = 1,
#                                maxGSSize = 50000,
#                                seed = 717)
#   
#   #### 保存结果 ---
#   saveRDS(res,
#           file = paste0("/home/caojt/RNAseq/figure/lanscape/多器官/",i,"/mito_gsea/","48_0_gsea.rds"))
#   write_csv(res@result %>% filter(pvalue <= 0.05),
#             file = paste0("/home/caojt/RNAseq/figure/lanscape/多器官/",i,"/mito_gsea/","48_0_gsea.csv")
#   )}
###设置工作路径------
setwd("~/精神疾病/data/388/Seurat/data/GSEA/")
output <- "~/精神疾病/data/388/Seurat/data/GSEA/CSV/"
###读取GSEA文件-----
tissues <- c("ASD", "BD", "MDD", "PTSD", "SZ")  
library(tidyverse)  

# 初始化adata为空，第一个CSV文件包含所有必要的列  
adata <- NULL  

# 循环读取和处理CSV文件  
for (i in tissues) {  
  file_path <- paste0("~/精神疾病/data/388/Seurat/data/GSEA/CSV/gsea_", i, "-table.csv")  
  
  # 读取CSV文件  
  data <- read.csv(file_path)  
  
  # 确保数据框中包含必要的列  
  if (!all(c("pvalue", "NES", "ID") %in% names(data))) {  
    stop("One or more of the required columns (pvalue, NES, ID) is missing in the CSV file.")  
  }  
  
  # 添加额外的列  
  data$a <- -log10(data$pvalue)  
  data$tissue <- i  
  data$change <- ifelse(data$NES > 0, "up", "down")  
  
  # 如果adata是空的，用当前的数据框初始化它  
  if (is.null(adata)) {  
    adata <- data  
  } else {  
    # 否则，将当前数据框添加到adata中  
    adata <- rbind(adata, data)  
  }  
}  
# 保存adata到CSV文件  
write.csv(adata, file = "~/精神疾病/data/388/Seurat/data/GSEA/CSV/allgsea_mental_diseases.csv", row.names = FALSE)


#####GSEA气泡图绘制-----
#####普遍上调-----
up_data <- subset(adata,subset = adata$change %in% "up")
a <- table(up_data$ID)
a <- as.data.frame(a)
a1 <- subset(a,subset = a$Freq > 3)
up_data1 <- subset(up_data,subset = ID %in% a1$Var1)
library(forcats)
up_data1$Description <- as.factor(up_data1$Description)
up_data1$Description <- fct_inorder(up_data1$Description)
up_data1 <- subset(up_data1,subset = pvalue < 0.05)
# down_data1 <- subset(down_data1,subset = Description %in% c("response to oxygen levels","dicarboxylic acid metabolic process","glucose metabolic process","fatty acid metabolic process","response to fatty acid","cellular carbohydrate metabolic process","energy derivation by oxidation of organic compounds","mitochondrial transmembrane transport","fatty acid oxidation","lipid oxidation","pyruvate metabolic process"))

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
###特异性下调-----
a1 <- subset(a,subset = a$Freq %in% 1)
up_data1 <- subset(up_data,subset = ID %in% a1$Var1)
library(forcats)
up_data1$Description <- as.factor(up_data1$Description)
up_data1$Description <- fct_inorder(up_data1$Description)
up_data1 <- subset(up_data1,subset = pvalue < 0.05)
# down_data1 <- subset(down_data1,subset = Description %in% c("response to oxygen levels","dicarboxylic acid metabolic process","glucose metabolic process","fatty acid metabolic process","response to fatty acid","cellular carbohydrate metabolic process","energy derivation by oxidation of organic compounds","mitochondrial transmembrane transport","fatty acid oxidation","lipid oxidation","pyruvate metabolic process"))

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
####失调线粒体基因upset图----
###读取差异基因-----
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
#默认
output <- "~/RNAseq/results/gsea/72_0/"
sets <- c(tissues <- c("Bl","Br","C", "H","K","Li","LN","Lu","P","SI","Sk","Sp","Th"))
pdf(file.path(paste0(output,"上调线粒体基因upset图.pdf")), width = 13.64, height = 9.54)

p1 <- upset(fromList(up_wide), 
            order.by = "freq",
            nsets = 100,
            nintersects = 30, #需要绘制的交集数目
            mb.ratio = c(0.7, 0.3),#柱状图与矩阵点图之间的比例大小
            # number.angles = 0,#柱子上方数字倾斜角度
            show.numbers = "yes",
            point.size = 1.8,#矩阵中圆圈的大小
            line.size = 0.8, #矩阵点图中点和线的大小
            sets.x.label = "Set Size", #柱状图的轴标签
            #main.bar.color = c("#73BAD6","#73BAD6","#73BAD6","#73BAD6","#73BAD6","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf"), #柱状图柱子颜色  #"#7dc3fe",
            main.bar.color = "#CD5C5C",
            sets.bar.color = c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#F0027F","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#6495ED",
                               "cyan1", "royalblue4","darksalmon"),
            # matrix.color = c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#F0027F"),  #交点的颜色
            matrix.color = "#CD5C5C",
            mainbar.y.label = "Gene Intersections", 
            text.scale = 1.75,
            shade.color = "#CD5C5C",
            shade.alpha = 0.1)
p1
dev.off()
####多组deseq2的热图------
####差异基因的计算------
###切换工作路径----
setwd("/home/caojt/RNAseq/")
output <- "./figure/results/gsea/"
dir.create(output, recursive = TRUE)
####差异基因的计算-----
#### 配置
dir.create(output, recursive = T)
#### 导入数据
# 导入counts文件
counts <-
  read.csv(
    file = "/home/gongfengcz/RNAseq/data/多器官/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv",
    row.names = 1
  )
# 导入样本信息
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
  #### 差异分析
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
      # 提取差异基因 提取实验组Pvs对照组N
      res <-
        results(dds, contrast = c("group", positive_group, negative_group))
      # 按照调整后P值排序
      res_ordered <- res[order(res$padj),]
      # 去掉NA
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
####折线图绘制-----
setwd("~/RNAseq/")
library(readr)
library(readxl)
library(UCell)
library(Seurat)
library(ggplot2)
vision <- read.csv("~/RNAseq/results/vision_mito.csv")
metadata <- read.csv("/home/lixy/SP-hope/data/expr-merged/20220725_metadata.csv")
#vision$X <- gsub(".gencode.vM19", "", vision$X)
# 假设你的矩阵名字是my_matrix，你想要操作的列名为"BAT_24"
split_vector <- strsplit(metadata$group, "-") # 使用下划线分割字符串

# 获取分割后的第一部分（即下划线前面的部分）
first_part <- sapply(split_vector, "[", 1)

# 将这个新的列添加到你的矩阵中
metadata$group <- first_part
# 假设你的矩阵名字是my_matrix，你想要修改的列是"group"
#metadata$group <- sub("Limb", "Limb_Muscle", metadata$group)
#metadata$group <- sub("Small", "Small_Intestine", metadata$group)
metadata$group <- paste0(metadata$group,"-",metadata$age)
vision1 <- merge(metadata, vision, by.x = "sample", by.y = "X", all = F)
vision1 <- vision1[,-c(2:4)] 
colnames(vision1)
vision1 <- aggregate(. ~ sample, data = vision1, mean)
vision1_long<-melt(vision1,
                   id.vars = c('sample'),#需要保留不参与聚合的变量,
                   variable.name='Pathways',
                   value.name='Score')
#vision1_long$group <- gsub("Limb_Muscle","LimbMuscle",vision1_long$group )
#ision1_long$group <- gsub("Small_Intestine","SmallIntestine",vision1_long$group )
vision1_long$tissue <- sapply(strsplit(vision1_long$sample, split = "-"), "[", 1)
vision1_long$time <- sapply(strsplit(vision1_long$sample, split = "-"), "[", 2)


###common通路的折线图-----
vision2_long <- subset(vision1_long,subset = vision1_long$Pathways %in% "OXPHOS")
vision2_long <- vision2_long[!grepl("NA", vision2_long$sample), ]
vision2_long <- subset(vision2_long,subset = vision2_long$tissue %in% c("C","K","Li","LN","SI","Sk"))

# 根据排序后的行索引对矩阵进行重新排列
vision2_long$time <- factor(vision2_long$time, levels = c('0', '6', '12', '24',"48","72","96"))
# 按指定顺序排序数据框
vision2_long <- vision2_long[order(vision2_long$time), ]

ggplot(vision2_long, aes(vision2_long$time, vision2_long$Score, group=tissue, color=tissue, shape=tissue))+
  geom_point(size = 2) +
  geom_line(size = 1) +
  # scale_linetype_manual(values = c(3, 1)) + #自定义线条形状
  scale_color_manual(values = c("#FDBE85", "#D94701", "#8856A7", "#9EBCDA", "#FAA75B", "#882E72", "#5289C7", "#5F4690", "#F04590", "#C63F6D", "#913C3E", "#566E3D", "#1D9183", "#2F65A5", "#243475")) + #自定义颜色
  # scale_shape_manual(values = c(17,18,19))+  #自定义点的形状
  theme_bw()+  #修改主题
  # theme(legend.title = element_blank(),#图例标题去除
  #       legend.text = element_text(family = 'Arial'),#字体
  #       legend.position = c(0.9,0.4),#位置
  #       legend.direction = "vertical")  +#水平或垂直
  labs(x = "Time", # 定义x轴文本
       y = "OXPHOS")# 定义y轴文本

###特异通路的折线图-----
vision2_long <- subset(vision1_long,subset = vision1_long$Pathways %in% "OXPHOS")
vision2_long <- vision2_long[!grepl("NA", vision2_long$sample), ]
vision2_long <- subset(vision2_long,subset = vision2_long$tissue %in% c("P","SI","Sk","Sp"))

# 根据排序后的行索引对矩阵进行重新排列
vision2_long$time <- factor(vision2_long$time, levels = c("0","6", '12', '24',"48","72","96"))
# 按指定顺序排序数据框
vision2_long <- vision2_long[order(vision2_long$time), ]

ggplot(vision2_long, aes(vision2_long$time, vision2_long$Score, group=tissue, color=tissue, shape=tissue))+
  geom_point(size = 2) +
  geom_line(size = 1) +
  # scale_linetype_manual(values = c(3, 1)) + #自定义线条形状
  scale_color_manual(values = c("grey", "grey", "grey", "firebrick", "#882E72", "#5289C7", "#5F4690", "#F04590", "#C63F6D", "#913C3E", "#566E3D", "#1D9183", "#2F65A5", "#243475")) + #自定义颜色
  # scale_shape_manual(values = c(17,18,19))+  #自定义点的形状
  theme_bw()+  #修改主题
  # theme(legend.title = element_blank(),#图例标题去除
  #       legend.text = element_text(family = 'Arial'),#字体
  #       legend.position = c(0.9,0.4),#位置
  #       legend.direction = "vertical")  +#水平或垂直
  labs(x = "Time", # 定义x轴文本
       y = "Mitophagy")# 定义y轴文本

