library(ggplot2)
library(aPEAR)
plotGK = function(df,outjpeg,mainTitle){
   p = ggplot(df, aes(x = reorder(df[,"Description"],  -log10(df[,"pvalue"])), y = -log10(df[,"pvalue"]))) +
        geom_point(aes(color = -log10(df[,"pvalue"]), size = Count )) +
        scale_color_gradient(low = "#E2A3A3", high = "#F54E4E")+
        coord_flip() +  # 翻转坐标轴，使GO_Term在纵轴上显示
        labs(title = mainTitle,
            x = "",
            y = expression(-log[10]("pvalue")),
            color = expression(-log[10]("pvalue")),
            size = "Count") +
        guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))+
    theme_minimal() +
    theme(
        panel.border = element_blank(),  # 去除整个面板边框
        axis.line = element_line(colour = "black"),  # 重新添加坐标轴线条
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(),  # 去除次要网格线
        axis.ticks = element_line(colour = "black"),  # 设置坐标轴刻度线颜色
        axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
        plot.title = element_text(hjust = 0.5),  # 标题居中
        axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
        axis.text.x = element_text(size = 15)
    )
    jpeg(width = 200, height = 150, units = "mm", res = 300, file = outjpeg)
        print(p)
    dev.off()
}
# df1 = read.csv("/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upGo_venn.csv",sep = "\t")
# plotGK(df1,"/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upGo.jpeg","")

df2 = read.csv("/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upKegg_venn.csv",sep = "\t")
# # # plotGK(df2,"/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upKegg.jpeg","")
p <- enrichmentNetwork(df2, drawEllipses = TRUE, fontSize = 7)
jpeg(width = 300, height = 321, units = "mm", res = 300, file = "/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upKegg_aPERA.jpeg")
    print(p)
dev.off()
df3 = read.csv("/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/downKegg.csv",sep = "\t")
# plotGK(df3,"/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/downKegg.jpeg","")
p <- enrichmentNetwork(df3, drawEllipses = TRUE, fontSize = 7)
jpeg(width = 321, height = 300, units = "mm", res = 300, file = "/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/downKegg_aPERA.jpeg")
    print(p)
dev.off()

