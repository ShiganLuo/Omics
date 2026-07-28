rm(list=ls());options(stringsAsFactors=FALSE)
setwd("/home/lsg/GWAS/Betta_splendens/Het_LROH")

het = read.table("Splendens410_Het.het",sep = "",header = TRUE)
roh = read.table("Splendens410_ROH.hom.indiv",sep = "",header = TRUE)
fid = read.csv("Pop410_cluster_final.txt",sep = "\t",header = TRUE)
het_fid = merge(het, fid, by.x = "FID", by.y = "IID", all = FALSE)
roh_fid = merge(roh,fid,by.x = "FID",by.y = "IID",all = FALSE)
head(het_fid)
####plot Froh result as a vioplot####
library(ggplot2)
library(reshape2)
require(cowplot)
fig_order = c("Smaragdina_guitar","Mahachaiensis","Siamorientalis","Smaragdina","Stiktos","Imbellis","Splendens",
"Black","Copper","Crowntail","Dragon","Dumbo_Halfmoon","Dumbo_HMPK","Fighter","Giant","Halfmoon","Koi","Orange",
"Red","Royal","Steel","Turquoise","Veiltail","White","Yellow"
)
het_fid$variety = factor(het_fid$variety,levels = fig_order)
roh_fid$variety = factor(roh_fid$variety,levels = fig_order)
###Veiltail,Crowntail,Fighter,Halfmoon,HMPK#######
fill1=c("#F0E68C", "#E6E6FA", "#FF4500", "#DB7093", "#00CED1", 
                  "#FF7F50", "#D2691E", "#8A2BE2", "#A52A2A", "#A9A9A9", 
                  "#00FF00", "#ADD8E6", "#FFFF00", "#90EE90", "#FFA500", 
                  "#4682B4", "#EE82EE", "#D3D3D3", "#40E0D0", "#FF6347", 
                  "#FAFAD2", "#800000", "#FF00FF", "#32CD32", "#F08080")
p1 = ggplot(data=het_fid, aes(x=variety, y=F,colour=variety))+
	geom_violin(aes(fill=factor(variety)))+
	scale_colour_manual(values = fill1)+
	scale_fill_manual(values = fill1)+
	xlab("variety")+
        ylab("F")+
	theme(axis.text.x=element_text(angle=60,hjust=0.65,vjust=0.6,size=8),
	axis.title=element_text(size=10),axis.text.y=element_text(size=8))+
	theme(legend.position="none") 
p2 = ggplot(data = roh_fid,aes(x = variety,y = KBAVG,colour = variety))+
	geom_violin(aes(fill = factor(variety)))+
	scale_colour_manual(values = fill1)+
	scale_fill_manual(values = fill1)+
	xlab("variety")+
	ylab("KBAVG")+
	theme(axis.text.x=element_text(angle=60,hjust=0.65,vjust=0.6,size=8),
	axis.title=element_text(size=10),axis.text.y=element_text(size=8))+
	theme(legend.position="none")
png("F_LROH.violin.png", height=1800, width=1800, res=300)
ggdraw() + 
	draw_plot(p1,0,0.5,1,0.5) + 
	draw_plot(p2,0,0,1,0.5) 
dev.off()
##############################

