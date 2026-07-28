library(ggplot2)

# 示例数据（每个方向的频数）
data <- data.frame(
  direction = factor(c("N", "NE", "E", "SE", "S", "SW", "W", "NW"), 
                     levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")),
  count = c(12, 23, 18, 15, 10, 8, 5, 9)
)

# 绘制玫瑰图
jpeg(width = 300, height = 300, units = "mm", res = 300, file = "rose_diagram.jpeg")
ggplot(data, aes(x = direction, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_polar(start = 0) +  # 极坐标
  theme_minimal() +
  labs(title = "Rose Diagram", x = "", y = "")
dev.off()
