#' Create a bar plot of hUSI scores and save as JPEG.
#'
#' @param vec Named numeric vector of hUSI scores.
#' @param outpng Path to save the output JPEG file.
#' @param base_width_mm Base device width in millimeters; grows with sample count.
#' @param base_height_mm Base device height in millimeters; grows with label length.
#' @param width_per_sample Extra width per bar to keep labels visible when many samples.
#' @param height_per_label Extra height per character of the longest label.
#' @param xlab_line Line offset for the x-axis title to avoid collisions with tick labels.
bar <- function(vec, outpng,
        base_width_mm = 150,
        base_height_mm = 120,
        width_per_sample = 20,
    height_per_label = 4,
    xlab_line = 3) {
    n <- length(vec)
    labels <- names(vec)
    max_label_len <- ifelse(length(labels) > 0, max(nchar(labels)), 5)

    # Scale device to fit labels: width grows with sample count, height with label length.
    width_mm <- max(base_width_mm, width_per_sample * n)
    height_mm <- max(base_height_mm, 80 + height_per_label * max_label_len)

    jpeg(width = width_mm, height = height_mm, units = "mm", res = 300, file = outpng)
    bottom_mar <- max(5, 2 + max_label_len * 0.3)  # add margin if labels are long
    par(mar = c(bottom_mar, 4, 2, 1))
    colors <- grDevices::rainbow(max(3, n))[seq_len(n)]
        barplot(vec, col = colors, xlab = "", ylab = "hUSI", cex.lab = 1.1, border = NA,
            las = 1, cex.names = min(0.9, 10 / max_label_len))
    mtext("sample", side = 1, line = xlab_line, cex = 0.9)
    dev.off()
}

stacked_bar <- function() {
    library(ggplot2)
library(dplyr)
# 读取并处理数据
data <- read.csv("test_data1.csv") %>%
  group_by(sample) %>%
  mutate(total = sum(cell_number),
         percentage = cell_number / total * 100) %>%
  ungroup()

# 获取第一个样本(CN_P21)的细胞类型顺序
first_sample_order <- data %>% 
  filter(sample == "CN_P21") %>% 
  pull(cell_type) %>% 
as.character()

# 设置细胞类型因子顺序（保持与第一个样本一致）
data$cell_type <- factor(data$cell_type, levels = first_sample_order)

# 自定义配色方案
custom_colors <- c(
"#bc90be", "#65a2d2", "#2bae9e", "#72b285", 
"#d4d98a", "#ffdf97", "#eb8b39", "#f58e87", "#e36146")

# 默认配色方案
ggplot(data, aes(x = sample, y = percentage, fill = cell_type)) +
  geom_col(position = "stack", width = 0.7,color = "white") +
  labs(x = "Sample", y = "Percentage (%)", fill = "Cell Type") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 自定义配色方案
ggplot(data, aes(x = sample, y = percentage, fill = cell_type)) +
  geom_col(position = "stack", width = 0.7,color = "white") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Sample", y = "Percentage (%)", fill = "Cell Type") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

sign_bar <- function() {
    # 读取数据；
dt <- read.csv("test.csv",header = T)
# 载入rstatix包；
library(rstatix)
library(dplyr)
# 按分组统计信息；
stats_summary <- dt %>% df_group_by(Group) %>% get_summary_stats()
# 对数据进行正态性检验(Shapiro-Wilk Normality Test)；
dt %>% shapiro_test(vars = "Expressions")
# 使用两种方法进行方差齐性检验；
# Levene’s test for homogeneity of variance；
dt$Group <- as.factor(dt$Group)
dt %>% levene_test(Expressions ~ Group)
# Bartlett test of homogeneity of variances；
bartlett.test(Expressions~Group,data = dt)

# 使用rstatix包进行方差分析；
dt %>% anova_test(Expressions~Group)
# 使用基础函数进行方差分析，整体来看差异显著;
oneway<-aov(Expressions~Group,data = dt)
anova(oneway)

# 使用Fisher LSD法进行均值比较；
# LSD法（Fisher’s Least Significant Difference）；
# LSD法检验微小的差异，比较方便的是直接得出显著性字母标记，不需人工标记；
# install.packages("agricolae")
library("agricolae")
out <- LSD.test(oneway,"Group",p.adj="bonferroni")
# 整理用于作图的数据框
rowname<-row.names(out$means)
mean<-out$means[,1]
sd<-out$means[,2]
marker<-out$groups$groups
plotdata<-data.frame(rowname,mean,sd,marker)
# ggplot2 绘制带显著性标记的柱状图
library("ggplot2")
p1<-ggplot(plotdata,aes(x=factor(rowname),y=mean))+geom_bar(position =position_dodge(0),fill="orange",width = 0.52,stat = "identity")
p2<-p1+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position=position_dodge(0.6),width=0.25)
p3<-p2+geom_text(aes(x=factor(rowname),y=mean+sd+0.1,label=marker),size=4,position = position_dodge(0.6))
p4<-p3+xlab("")+ylab("Lesion diameter (mm)")
p5<-p4+scale_y_continuous(limits = c(0, 4),expand=expansion(add = c(0, 0)))
# 更改y轴显示范围，这里的expand默认为TRUE。
mytheme<-theme_classic()+theme(axis.title = element_text(size = 12),
                               axis.text = element_text(size=12),
                               panel.grid.major = element_line(color = "white"),
                               panel.grid.minor = element_line(colour = "white"),
                               axis.text.x = element_text(size = 12,angle=45,vjust=0.7,hjust=0.8,color = "black"),
                               axis.text.y = element_text(size = 12,color = "black"),
                               legend.text = element_text(size = 12),legend.title = element_blank(),
                               plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
p5+mytheme
}