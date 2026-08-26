#### 加载包 ----
library(Seurat)
library(GSVA)
library(GSEABase) # 用于读取gmt格式文件
library(limma)
library(tidyverse)
library(ggplot2)
# library(ComplexHeatmap)


#### 配置 ----
setwd('~/精神疾病/data/388/Seurat/data/GSVA')
outdir_data <- paste0('./score/')
dir.create(outdir_data,recursive = T)
outdir_plots <- paste0('./plots/')
dir.create(outdir_plots,recursive = T)

#### 导入数据 ----
# 读取基因集数据库
s.sets <- getGmt("/DATA/public/MSigDB/7.5.1/h.all.v7.5.1.symbols.gmt")

sub.seurat.obj <- readRDS('~/精神疾病/data/388/Seurat/data/merged/merged_all_scted.rds')
sub.seurat.obj2 <- readRDS('~/精神疾病/data/388/Seurat/data/SCTed_data.rds')

# DefaultAssay(seurat.obj) <- "RNA"
# Idents(seurat.obj) <- 'big_cell_type'

# Cell <- c('Macrophages') # 'Macrophages'
# table(sub.seurat.obj$big_cell_type)
# Macrophages 
# 2994

# data是处理后的数据，counts是原始数据。提取为矩阵
expr <- as.matrix(sub.seurat.obj@assays$RNA@data)
## 进行gsva打分计算
# 结果中，行是不同通路，列是不同barcode
es.matrix <- gsva(
  expr,
  s.sets,
  min.sz = 10,# 基因集中最小基因数，默认10，表示分析中仅考虑包含大于10个基因的基因集
  max.sz = Inf,# 基因集中最大基因数，默认inf，无限
  tau = 1, # 计算权重  默认情况下，当method=“gsva”时，这个tau=1；
  #当method=“ssgsea”时，tau=0.25
  method = "gsva",
  abs.ranking = FALSE,#是否按基因表达绝对值排序
  verbose = TRUE,#是否输出结果
  parallel.sz = 1 # 并行计算的核数 ####可以改成1,1更快。原先是4.
)  # 默认情况下，kcdf=“Gaussian”;
# 当输入表达式值是整数计数时（RNA@counts），则应 kcdf="Poisson"
# saveRDS(es.matrix, paste0(outdir_data,"GSVA_hallmarker_matrix.rds"))
es.matrix <- readRDS(paste0(outdir_data,"GSVA_hallmarker_matrix.rds"))

# file_path <- file.path(path.expand("~/精神疾病/data/388/Seurat/data/GSVA/score"), "GSVA_hallmarker_matrix.rds")  
# # 使用saveRDS函数保存对象  
# saveRDS(es.matrix, file_path)

#### limma差异分析 ----
# 查看组别
table(sub.seurat.obj$group)
# CON  ASD   BD  MDD PTSD   SZ 
# 2719 2378 2868 1551 2513 2838 


# sub.seurat.obj$group <- sub.seurat.obj$orig.ident
# 
# # 重新设置'group'的水平顺序  
# new_levels <- c("CON", "ASD", "BD", "MDD","PTSD","SZ") 
# sub.seurat.obj$group <- factor(sub.seurat.obj$group, levels = new_levels)  
# 
# print(levels(seurat_obj$group))
# table(sub.seurat.obj$group)

# 设置分组
case <- "SZ"
control <- "CON"
# 提取指定的两列数据，行名是barcode
meta <-
  sub.seurat.obj@meta.data[, c("orig.ident",
                               "group")]
# 添加cell_id列，列的内容是barcode
meta$cell_id <- rownames(meta)

## 根据不同分组，提取各自包括的cell_id（barcode）
cell_id.1 <-
  meta[meta$group == case, ]$cell_id
cell_id.2 <-
  meta[meta$group == control, ]$cell_id

## 根据提取的不同组cell_id，将不同组的gsva得分信息提取出来
es.matrix.1 <-
  as.data.frame(es.matrix[, cell_id.1],
                row.names = row.names(es.matrix))
es.matrix.2 <-
  as.data.frame(es.matrix[, cell_id.2],
                row.names = row.names(es.matrix))
# 合并两组的得分信息，为一个矩阵
es.matrix.f <- cbind(es.matrix.1, es.matrix.2)

# 相比于gsva后的结果初始数据，合并所得的数据类型发生改变,并且列是按不同分组来排列的
class(es.matrix) # gsva后的结果初始数据
class(es.matrix.f) # 合并所得的数据


### 分组设计
# 提取
grouP <-
  c(rep("case", dim(es.matrix.1)[2]), rep("control", dim(es.matrix.2)[2])) %>% as.factor()

## 分组。1-表示属于该组；0-表示不属于该组
desigN <- model.matrix( ~ grouP + 0)
# 添加行名
rownames(desigN) <-
  c(colnames(es.matrix.1), colnames(es.matrix.2))
# desigN 比较
comparE <-
  makeContrasts(grouPcase - grouPcontrol, levels = desigN)
# 拟合
fiT <- lmFit(es.matrix.f, desigN)
# 降维
fiT2 <- contrasts.fit(fiT, comparE)
# 统计分析
fiT3 <- eBayes(fiT2)
# 差异分析
diff <- topTable(fiT3, coef = 1, number = 200)####做gobp打分时可以改一下参数，改大点，目前太小

## 得到差异分析的结果
t_results <-
  as.data.frame(diff$t, row.names = rownames(es.matrix.f))
head(t_results)
colnames(t_results) <- c("t_value") # 调整列名
# 保存t值结果
file_path <- file.path(path.expand("~/精神疾病/data/388/Seurat/data/GSVA/score"), "SZ_GSVA_hallmarker_t_value.rds")  
saveRDS(t_results, file_path)



######### 可视化 ----
Cell <- c('SZ') 
#t_results <- readRDS(paste0(outdir_data,"ASD_GSVA_hallmarker_t_value.rds"))

### 一一分散条形图展示 ----
## 对通路的名字格式进行调整
rownames(t_results) <- gsub("HALLMARK_", "", rownames(t_results))
rownames(t_results) <- gsub("_", " ", rownames(t_results))
rownames(t_results) <- tolower(rownames(t_results))
## 增加列，并改名
t_results$hallmark <- rownames(t_results)
colnames(t_results) <- c("t", "hallmark")
head(t_results)

## 根据不同t值，分不同组，给定不同色
t_results$hallmark = with(t_results, reorder(hallmark, t))
t_results$fill <- ""
t_results[t_results$t >= 2.58,]$fill <-  # t值绝对值大于2.58，上调
  "up"
t_results[t_results$t <= -2.58,]$fill <-  # 下调
  "down"
t_results[abs(t_results$t) < 2.58,]$fill <-  # 不明显
  "no"
t_results$color <- ""
t_results[abs(t_results$t) < 2.58,]$color <-
  "n"
t_results[abs(t_results$t) >= 2.58,]$color <- # t值绝对值大于2.58，具有显著性
  "y"

## 柱形图
p <-
  ggplot(t_results, aes(x = hallmark, y = t, fill = fill)) +
  geom_bar(stat = "identity", width = .8) +
  scale_fill_manual(
    values = c(
      "down" = "#6fa7d3",
      "up" = "#da535a",
      "no" = "#cccccc"
    ),
    guide = F
  ) +
  geom_hline(
    yintercept = c(-2.58, 2.58),
    color = "white",
    linetype = "dashed",
    size = 0.5
  ) +
  coord_flip() +
  xlab("") +
  ylab("t value of GSVA score") +
  geom_text(
    data = subset(t_results, t < 0), ## 下调组
    aes(
      x = hallmark,
      y = 0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 3,
    hjust = "inward" # outward 在条内  inward 在条外
  ) +
  geom_text(
    data = subset(t_results, t > 0), ## 上调组
    aes(
      x = hallmark,
      y = -0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 3,
    hjust = "outward" # inward 在条内  outward 在条外
  ) +
  scale_colour_manual(values = c("y" = "black", "n" = "#cccccc"), # y-
                      guide = FALSE) +
  labs(title=paste0(Cell,' vs. CON'))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust=0.5,size = 12),
    axis.text.y = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),
  )
p

ggsave(
  filename = paste0(outdir_plots,"SZ_GSVA_bar.pdf"),
  plot = p,
  width = 12, # tu 5.3  mac 4.4
  height = 6
)
