library(circlize)
# 读取数据
df <- read.csv("test_data2.csv")
# 创建邻接矩阵
nodes <- sort(unique(c(df$source, df$target)))
mat <- matrix(0, nrow = length(nodes), ncol = length(nodes),
              dimnames = list(nodes, nodes))

for(i in 1:nrow(df)) {
mat[as.character(df$source[i]), as.character(df$target[i])] <- df$count[i]
}

# 设置颜色方案
colors <- c("#7ac770", "#c8cd4c", "#eccf77", "#ffdf97", 
            "#f58e87", "#e18383", "#e36146", "#52bdb9", "#65a2d2")
# 分配节点颜色（循环使用9种颜色）
node_colors <- setNames(colors[(seq_along(nodes) - 1) %% length(colors) + 1], nodes)

# 绘制弦图
circos.clear()
circos.par(start.degree = 90, gap.degree = 2)

chordDiagram(
x = mat,
grid.col = node_colors,
directional = 1,
direction.type = "arrows",
link.arr.type = "big.arrow",
transparency = 0.25,
annotationTrack = c("grid", "name"),
preAllocateTracks = list(track.height = 0.1)
)

# 添加图例（可选）
legend("right", legend = names(node_colors), fill = node_colors, 
       title = "Nodes", cex = 0.6)