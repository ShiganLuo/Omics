library(ggplot2)
library(ggpubr)
plotAletenedFeature <- function(df1,df2){

  Top_altered_RNA <- ggplot(df1,aes(x=df1$feature,y=1,fill = df1$number))+geom_bar(stat = "identity")+
    coord_flip()+
    scale_fill_gradient2(low="light blue",high="dark red",mid = "red",midpoint = 0.9,limits = c(0,1))+
    theme_bw()+
    theme(
      plot.margin = unit(x=c(15,0,10,0),units="pt"),
      legend.position="none",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face="bold", color="black", size=18),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_text(color="black",family = "Arial", size=10),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  Top_altered_DNA <- ggplot(df2,aes(x=df2$feature,y=1,fill = df2$number))+geom_bar(stat = "identity")+
    coord_flip()+theme_bw()+
    scale_fill_gradient2(low="light blue",high="dark red",mid = "red",midpoint = 0.9,limits = c(0,1))+
    theme(
      plot.margin = unit(x=c(15,0,10,0),units="pt"),
      legend.position="left",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face="bold", color="black", size=18),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_text(color="black",family = "Arial", size=10),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color="black", size=10, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  p <- ggarrange(Top_altered_DNA,Top_altered_RNA,align = "h",widths = c(7,1))
  ggsave(p,filename = "Top_altered_genes_20220918.png",width = 2.25,height = 4.01,device = "png")
}
infile = "/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/TopAlterFeature.csv"
df = read.table(infile,sep=',',header=TRUE)
df_TE = df[df['type'] == 'TE',]
df_gene = df[df['type'] == 'gene_id',]
plotAletenedFeature(df_TE,df_gene)