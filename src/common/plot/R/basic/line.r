library(ggplot2)
BIM = read.csv("/home/lsg/GWAS/Betta_splendens/pi/BIM_20K10S_chr5.windowed.pi", header = TRUE, sep = "\t")
BMA = read.csv("/home/lsg/GWAS/Betta_splendens/pi/BMA_20K10S_chr5.windowed.pi", header = TRUE, sep = "\t")
BSM = read.csv("/home/lsg/GWAS/Betta_splendens/pi/BSM_20K10S_chr5.windowed.pi", header = TRUE, sep = "\t")
BSP = read.csv('/home/lsg/GWAS/Betta_splendens/pi/BSP_20K10S_chr5.windowed.pi', header = TRUE, sep = '\t')
BSS = read.csv("/home/lsg/GWAS/Betta_splendens/pi/BSS_20K10S_chr5.windowed.pi", header = TRUE, sep = "\t")
BIM$Species =  "Imbellis"
BMA$Species = "Mahachaiensis"
BSM$Species = "Smaragdina"
BSP$Species = "Splendens"
BSS$Species = "Siamorientalis"
jpeg(width = 300, height = 150, units = "mm", res = 300, file = "five.jpeg")
data = rbind(BIM, BMA, BSM, BSP, BSS)
ggplot(data, aes(x = BIN_START, y = PI, color = Species))+
    geom_smooth(method = "gam", se = FALSE, size = 0.8)+
    # geom_smooth(data = BIM, mapping = aes(x = BIN_START, y = PI, color = Species), method = "loess", se = FALSE, size = 0.8, color = "darkorange2")+
    # geom_smooth(data = BMA, mapping = aes(x = BIN_START, y = PI, color = Species), method = "loess", se = FALSE, size = 0.8, color = "salmon")+
    # geom_smooth(data = BSM, mapping = aes(x = BIN_START, y = PI, color = Species), method = "loess", se = FALSE, size = 0.8, color = "steelblue2")+
    # geom_smooth(data = BSP, mapping = aes(x = BIN_START, y = PI, color = Species), method = "loess", se = FALSE, size = 0.8, color = "#BC80BD")+
    # geom_smooth(data = BSS, mapping = aes(x = BIN_START, y = PI, color = Species), method = "loess", se = FALSE, size = 0.8, color = "lightseagreen")+
    theme_minimal()+
    theme(
        panel.background = element_rect(fill = "white", color = "black"), # 设置绘图区域背景为白色，边框为黑色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        axis.line = element_line(color = "black")  # 设置坐标轴线为黑色
    )+
    scale_color_manual(values = c("Imbellis" = "darkorange2" , "Mahachaiensis" = "salmon", "Smaragdina" = "steelblue2", "Splendens" = "#BC80BD", "Siamorientalis" = "lightseagreen"))+
    labs(x = "sequence", y = "Pi", color = "Species")
dev.off()
