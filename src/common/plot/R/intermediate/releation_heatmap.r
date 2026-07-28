# 加载必要的包
library(psych)
library(ComplexHeatmap)
library(circlize)

# 读取并预处理数据（与之前相同）
env <- read.csv("env.csv", header = TRUE)
env_factors <- env$Factor
env_data <- env[, -1]
colnames(env_data) <- gsub("-", ".", colnames(env_data))
env_data_t <- as.data.frame(t(env_data))
colnames(env_data_t) <- env_factors

expressions <- read.csv("expressions.csv", header = TRUE)
gene_names <- expressions$gene
expressions_data <- expressions[, -1]
expressions_data_t <- as.data.frame(t(expressions_data))
colnames(expressions_data_t) <- gene_names

common_samples <- intersect(rownames(env_data_t), rownames(expressions_data_t))
env_data_t <- env_data_t[common_samples, ]
expressions_data_t <- expressions_data_t[common_samples, ]

# 计算相关性矩阵和p值矩阵
cor_results <- corr.test(env_data_t, expressions_data_t, method = "pearson")
r_matrix <- cor_results$r
p_matrix <- cor_results$p

# 设置颜色映射和热图参数
max_abs <- max(abs(r_matrix))
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("#0f86a9", "white", "#FC8452"))
col_fun2 = colorRamp2(c(-max_abs, 0, max_abs), c("#A5CC26", "white", "#FF7BAC"))
col_fun3 = colorRamp2(c(-max_abs, 0, max_abs), c("#B3A9EB", "white", "#ffa500"))
# 定义星号标注函数（优化显示位置）
add_stars <- function(j, i, x, y, width, height, fill){
  p_value <- p_matrix[i, j]
if (is.na(p_value)) return()
if (p_value < 0.001) {
    grid.text("***", x, y, vjust=0.7,gp = gpar(fontsize = 12))
  } elseif (p_value < 0.01) {
    grid.text("**", x, y, vjust=0.7,gp = gpar(fontsize = 12))
  } elseif (p_value < 0.05) {
    grid.text("*", x, y, vjust=0.7,gp = gpar(fontsize = 12))
  }
}

# 绘制热图（修正参数问题）
ht <- Heatmap(
  matrix = r_matrix,
  name = "Cor",
  col = col_fun2,
  rect_gp = gpar(col = "white", lwd = 1.5),
  row_dend_width = unit(1.5, "cm"),
  column_dend_height = unit(1.5, "cm"),
  column_dend_gp = gpar(col = "gray30",lwd = 1.2),
  row_dend_gp = gpar(col = "gray30",lwd = 1.2),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  row_split = 2,
  column_split = 2,
  row_gap = unit(1, "mm"), # 行间隙大小
  column_gap = unit(1, "mm"), 
  row_title = NULL,column_title = NULL,
# 修正单元格尺寸设置方式
  heatmap_width = unit(ncol(r_matrix)*0.62, "cm"), 
  heatmap_height = unit(nrow(r_matrix)*0.8, "cm"), 
  cell_fun = add_stars)

# 输出图形
draw(ht)