#!/usr/bin/env Rscript
rm(list=ls());options(stringsAsFactors=FALSE)
file<-dir(path="/home/lsg/GWAS/Betta_splendens240901/interpopulation/output/pi",pattern=".windowed.pi")

# print(file)
# bacoground = GWAS


###Veiltail,Crowntail,Fighter,Halfmoon,HMPK#######

data = vector("list",length(file))
Nd = vector("list",length(file))

for(i in 1:length(file)){
        mainplot=unlist(strsplit(file[i],".pi.windowed"))[1]
        # print(strsplit(file[i],".sites"))
        print(mainplot)
        file[i] = paste0("/home/lsg/GWAS/Betta_splendens240901/interpopulation/output/pi/",file[i]) 
        # print(file[i])
        dat = read.table(file[i],header=TRUE)
        # #dat = dat[dat$N_VARIANTS > 10,]
        # print(class(dat$PI))
        Nd[i] = mean(dat$PI)
        Outvals = boxplot(dat$PI)$out
        # print(Outvals)
        dat = dat[which(!dat$PI %in% Outvals),]

        data[[i]] = data.frame(Pop = rep(mainplot, length(dat$PI)),Pi = dat$PI)
}
print("--------------------------")
fill1=c("Imbellis" = "#E1FF00","Mahachaiensis" = "#FA8072",
"Smaragdina_guitar" = "#660000","Smaragdina" = "#8A2BE2","Splendens" = "#BC80BD","Siamorientalis" = "#20B2AA",
"Stiktos" = "#735380")
data_all = do.call("rbind",data)
data_all$Pop=factor(data_all$Pop,order=TRUE,levels=c("Splendens","Imbellis","Siamorientalis","Smaragdina","Smaragdina_guitar","Stiktos","Mahachaiensis"))

data_all = na.omit(data_all)
# data_all$Pop=factor(data_all$Pop,order=TRUE,levels=c("Veiltail","Crowntail","Fighter","Halfmoon","HMPK"))

library(ggplot2)
library(reshape2)
library(cowplot)
#require(cowplot)

# png("BS_nucleotide_diversity.png", height=1000, width=1000, res=300)
p1 = ggplot(data=data_all, aes(x=Pop,y=Pi,colour=Pop))+
        geom_boxplot(aes(fill=factor(Pop)),coef = 1.5)+
		geom_boxplot(outlier.shape = NA)+
        scale_colour_manual(values = fill1)+
        scale_fill_manual(values = fill1)+
        xlab("Population")+
        ylab("Nucleotide_diversity")+
	theme(axis.text.x=element_text(angle=60,hjust=0.55,vjust=0.6,size=8),
	axis.title=element_text(size=10),axis.text.y=element_text(size=8))+
	theme(legend.position="none",
    panel.background = element_rect(fill = "white", color = "black"),)
#p1
#dev.off()

########plot for sub pop########

###Red,Orange,Yellow,Turquoise,Royal,Steel,Black,White,Dumbo,Giant,Koi###
fill2=c("Halfmoon" = "#A2B61E" ,"Dragon" = "#177D7A","Fighter" = "#907C62","Steel" = "#4682B4",
"Orange" = "#F7953B","Yellow" = "#E6D214","White" = "#EFDFBB","Koi" = "#DE77AE","Turquoise" = "#2ED29F",
"Giant" = "#6F5CCC","Crowntail" = "#91FF00","Red" = "#C63731","Copper" = "#00FFF7","Dumbo_HMPK" = "#F6D5A4",
"Dumbo_Halfmoon" = "#C47664","Royal" = "#0F55F1","Veiltail" = "#AED0CD","Black" = "#BEBFC5")
newdata = do.call("rbind",data)
newdata$Pop = factor(newdata$Pop,order=TRUE,levels=c("Black","Copper","Crowntail","Dragon","Dumbo_Halfmoon","Dumbo_HMPK","Fighter","Giant","Halfmoon","Koi","Orange",
"Red","Royal","Steel","Turquoise","Veiltail","White","Yellow"))
newdata = na.omit(newdata)

# png("BS_nucleotide_diversity_subpop.png", height=1000, width=1000, res=300)
p2 = ggplot(data=newdata, aes(x=Pop,y=Pi,colour=Pop))+
        geom_boxplot(aes(fill=factor(Pop)),coef = 1.5)+
		geom_boxplot(outlier.shape = NA)+
        scale_colour_manual(values = fill2)+
        scale_fill_manual(values = fill2)+
        xlab("Sub-population")+
        ylab("")+
        ylim(0,0.006)+
	theme(axis.text.x=element_text(angle=60,hjust=0.55,vjust=0.6,size=8),
	axis.title=element_text(size=10),axis.text.y=element_text(size=8))+
	theme(legend.position="none",
    panel.background = element_rect(fill = "white", color = "black"),)
# dev.off()
###merge two plot ############################
png("BS_nucleotide_diversity_merge.png", height=1200, width=2000, res=300)
ggdraw() +
        draw_plot(p1,0,0,0.5,1) +
        draw_plot(p2,0.5,0,0.5,1)
################################
# p1 + p2
#p2
dev.off()
###############################



