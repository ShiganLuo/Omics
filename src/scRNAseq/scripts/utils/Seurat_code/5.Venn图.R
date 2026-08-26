library (VennDiagram)  
library(openxlsx)

#数值导入，可对数值进行配对
ASD_up <- read_csv("~/精神疾病/data/388/Seurat/data/DEG/ASD_up_deg.csv")
BD_up <- read_csv("~/精神疾病/data/388/Seurat/data/DEG/BD_up_deg.csv")
MDD_up <- read_csv("~/精神疾病/data/388/Seurat/data/DEG/MDD_up_deg.csv")
PTSD_up <- read_csv("~/精神疾病/data/388/Seurat/data/DEG/PTSD_up_deg.csv")
SZ_up <- read_csv("~/精神疾病/data/388/Seurat/data/DEG/SZ_up_deg.csv")

#数据转置，如果不转后头函数venn.diagram对矩阵数值不识别
set1 = t(ASD_up)
set2 = t(BD_up)
set3 = t(MDD_up)
set4 = t(PTSD_up)
set5 = t(SZ_up)

head(up)

set1 <- na.omit(set1)  # 删除包含NA的行
set2 <- na.omit(set2)  # 删除包含NA的行
set3 <- na.omit(set3)  # 删除包含NA的行
set4 <- na.omit(set4)  # 删除包含NA的行
set5 <- na.omit(set5)  # 删除包含NA的行



set1 <- ASD_up[[gene_column]]
set2 <- BD_up[[gene_column]]
set3 <- MDD_up[[gene_column]]
set4 <- PTSD_up[[gene_column]]
set5 <- SZ_up[[gene_column]]

#五元#

venn.diagram(x=list(set1,set2,set3,set4,set5),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC"), # 填充色 配色https://www.58pic.com/
             
             category.names = c("ASD_down", "BD_down","MDD_down","PTSD_down","SZ_down") , #标签名
             
             cat.dist = c(0.2, 0.2, 0.2, 0.2, 0.2), # 标签距离圆圈的远近
             
             cat.pos = c(0, -10, 240, 120, 20), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col=c('black','black',"black","black", "black"),   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='~/精神疾病/data/388/Seurat/picture/Venn_down1.png',# 文件保存
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 300,  # 分辨率
             
             compression = "lzw"# 压缩算法
             
)

grid.draw(data)


#去除含label的列
# library(dplyr)
# ASD_down <- select(ASD_down, -label)
# BD_down <- select(BD_down, -label)
# MDD_down <- select(MDD_down, -label)
# PTSD_down <- select(PTSD_down, -label)
# SZ_down <- select(SZ_down, -label)


