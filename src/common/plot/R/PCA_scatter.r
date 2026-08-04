# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取分组信息
group_info <- read.csv("Group_info.csv", stringsAsFactors =  TRUE)

# 读取基因表达数据
expression_data <- read.csv("PCA_testdata.csv", row.names =  1)
expression_matrix <- t(expression_data)  # 转置矩阵（样本为行，基因为列）

# 执行PCA分析（自动进行中心化和标准化）
pca_result <- prcomp(expression_matrix, scale. =  TRUE)
pca_scores <-  as.data.frame(pca_result$x[,  1:2])
colnames(pca_scores) <- c("PC1",  "PC2")

# 添加样本分组信息
pca_scores$Sample <- rownames(pca_scores)
pca_data <- merge(pca_scores, group_info, by =  "Sample")
pca_data$Group <- factor(pca_data$Group)

# 设置分组颜色方案
group_colors <- c("GP1"  =  "#02afca",  "GP2"  =  "#b379b4",  "GP3"  =  "#4d97cd")
group_colors <- c("GP1"  =  "#efca72",  "GP2"  =  "#93cc82",  "GP3"  =  "#88c4e8")

# 预先计算主成分的方差解释百分比
pc_var <- summary(pca_result)$importance
pc1_var <- round(pc_var[2,  1] *  100,  1)
pc2_var <- round(pc_var[2,  2] *  100,  1)

# 创建坐标轴标签
x_label <- paste0("PC1 (", pc1_var,  "%)")
y_label <- paste0("PC2 (", pc2_var,  "%)")

# 样式1：添加置信椭圆的PCA图（填充颜色与描边相同，透明度0.6）
p_ellipse <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
   geom_point(size =  3, alpha =  0.8) +
   stat_ellipse(
      geom =  "polygon",  # 绘制填充多边形
      alpha =  0.6,  # 设置透明度
      level =  0.95,  # 95%置信区间
      size =  0.8,  # 边界线粗细
      show.legend =  FALSE
   ) +
   scale_color_manual(values = group_colors) +
   scale_fill_manual(values = group_colors) +
   labs(
      title =  "PCA with Confidence Ellipses (95%)",
      x = x_label,
      y = y_label
   ) +
   theme_bw(base_size =  12) +
   theme(
      panel.grid.minor = element_blank(),
      legend.position =  "right",
      plot.title = element_text(hjust =  0.5)
   )

# 样式2：添加多边形轮廓线的PCA图（填充颜色与描边相同，透明度0.6）
hull_data <- pca_data %>%
   group_by(Group) %>%
   slice(chull(PC1, PC2))  # 计算每组凸包点

p_hull <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
   geom_point(size =  3, alpha =  0.8) +
   geom_polygon(
      data = hull_data,
      alpha =  0.6,  # 填充透明度
      size =  0.8,  # 轮廓线粗细
      show.legend =  FALSE
   ) +
   scale_color_manual(values = group_colors) +
   scale_fill_manual(values = group_colors) +
   labs(
      title =  "PCA with Convex Hull Polygons",
      x = x_label,
      y = y_label
   ) +
   theme_bw(base_size =  12) +
   theme(
      panel.grid.minor = element_blank(),
      legend.position =  "right",
      plot.title = element_text(hjust =  0.5)
   )

# 显示图形
print(p_ellipse)
print(p_hull)