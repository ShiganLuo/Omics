# library(tidyverse)
library(dplyr)
library(ggplot2)
rm(list=ls())
args = commandArgs(trailingOnly = TRUE) 
base_dir = "/home/lsg/GWAS/Betta_splendens/overlap_gene/output"
infile = file.path(base_dir,args[1],"GO/data1_genelist.difgoall.filt.merg")
suffix = paste(args[1],"_go.jpeg",sep="")
outfile = file.path(base_dir,args[1],"plot",suffix)
print(outfile)
df = read.table(infile,header = TRUE,sep = "\t")
jpeg(width = 250, height = 300, units = "mm", res = 300, file = outfile)

# 绘制点图
ggplot(df, aes(x = reorder(df[,"GO_Term"], df[,"x1"]/df[,"n"]), y = df[,"AdjustedPv"])) +
  geom_point(aes(color = AdjustedPv, size = x1 )) +
  scale_color_gradient(low = "blue", high = "red")+
  coord_flip() +  # 翻转坐标轴，使GO_Term在纵轴上显示
  labs(title = "GO Enrichment Analysis",
       x = "GO Terms",
       y = "GeneRatio",
       color = "AdjustedPv",
       size = "Count") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20,hjust = 0.5),
    axis.title.x = element_text(size = 18),             # 调整x轴标签字体大小
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15))
# axis.text.x = element_text(angle = 45, hjust = 1)
dev.off()
