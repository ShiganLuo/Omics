setwd("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/VS")
library(clusterProfiler)
library(org.Hs.eg.db)


TF_MDD_UP<-read_csv("~/精神疾病/data/388/Seurat/data/DEG/MDD_up_deg.csv")$gene

MDD_TFS<- read.table("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/MDD/tfs_targer.tsv", header=T, sep=",")


filtered_MDD_TFS <- subset(MDD_TFS, tf %in% TF_MDD_UP)

target_gene_of_MDD_UP<-filtered_MDD_TFS$target_gene


#### 2.功能富集 ----

##_________________________________________________________________________________________________

case <- 'MDD '

deg_up<-target_gene_of_MDD_UP
bp <-
  enrichGO(
    deg_up,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

term <- bp@result


#保存结果
write.table(term,
            file = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/VS/target_gene_MDD_UP_GOBP.csv",
            quote = F,
            sep = ",",
            row.names = F
)
#绘制TOP20
df1 <- term[1:20,]####提取上调通路结果文件中的前二十条
df1$labelx=rep(0,nrow(df1))
df1$labely=seq(nrow(df1),1)
p <- ggplot(data = df1, 
            aes(x = -log10(pvalue),
                y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity",
           alpha=1,
           fill= "pink",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = df1$Description),
            size=3.5,
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  ggtitle(case)+
  scale_x_continuous(expand = c(0,0))
p
#保存
ggsave("~/brain/SCENIC/by_group/downsrtream_analysis/GO/target_gene_of_MDD_UP_GOBP.pdf", p,
       height=6,
       width=8)

###________________________________________________________________________________________________________

TF_BD_UP<-read_csv("~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/BDVSCON_up_deg.csv")$gene

BD_TFS<- read.table("~/brain/SCENIC/by_group/BD/tfs_targer.tsv", header=T, sep=",")


filtered_BD_TFS <- subset(BD_TFS, tf %in% TF_BD_UP)

target_gene_of_BD_UP<-filtered_BD_TFS$target_gene

#### 2.功能富集 ----

##_________________________________________________________________________________________________

case <- 'BD'

deg_up<-target_gene_of_BD_UP

# markers<-read.table("~/brain/差异/output/deg/ASDVSCON_deg.txt",
#                     header = T,sep=",")
# 
# deg_up <- filter(markers,
#                  markers$p_val < 0.05 & markers$avg_log2FC > 0.25)
# markers[markers$change %in% 'UP', ]$gene
# head(deg)
# #重新排序
# top.genes <-deg
# deg[order(deg$avg_log2FC, decreasing = T),]$gene
# class(deg)
# deg<-as.data.frame(deg)
bp <-
  enrichGO(
    deg_up,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

term <- bp@result
#保存结果
write.table(term,
            file = "~/brain/SCENIC/by_group/downsrtream_analysis/GO/target_gene_of_BD_UP_GOBP.csv",
            quote = F,
            sep = ",",
            row.names = F
)
#绘制TOP20
df1 <- term[1:20,]####提取上调通路结果文件中的前二十条
df1$labelx=rep(0,nrow(df1))
df1$labely=seq(nrow(df1),1)
p <- ggplot(data = df1, 
            aes(x = -log10(pvalue),
                y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity",
           alpha=1,
           fill= "pink",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = df1$Description),
            size=3.5,
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  ggtitle(case)+
  scale_x_continuous(expand = c(0,0))
p
#保存
ggsave("~/brain/SCENIC/by_group/downsrtream_analysis/GO/target_gene_of_BD_UP_GOBP.pdf", p,
       height=6,
       width=8)
###________________________________________________________________________________________________________

####__________________________________________________________________________________________________________________________

####__________________________________________________________________________________________________________________________

####__________________________________________________________________________________________________________________________

####__________________________________________________________________________________________________________________________

####__________________________________________________________________________________________________________________________

TF_ASD_DOWN<-read_csv("~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/ASDVSCON_down_deg.csv")$gene

ASD_TFS<- read.table("~/brain/SCENIC/by_group/ASD/tfs_targer.tsv", header=T, sep=",")


filtered_ASD_TFS <- subset(ASD_TFS, tf %in% TF_ASD_DOWN)

target_gene_of_ASD_DOWN<-filtered_ASD_TFS$target_gene

#### 2.功能富集 ----

##_________________________________________________________________________________________________

case <- 'ASD'

deg_down<-target_gene_of_ASD_DOWN

# markers<-read.table("~/brain/差异/output/deg/ASDVSCON_deg.txt",
#                     header = T,sep=",")
# 
# deg_up <- filter(markers,
#                  markers$p_val < 0.05 & markers$avg_log2FC > 0.25)
# markers[markers$change %in% 'UP', ]$gene
# head(deg)
# #重新排序
# top.genes <-deg
# deg[order(deg$avg_log2FC, decreasing = T),]$gene
# class(deg)
# deg<-as.data.frame(deg)
bp <-
  enrichGO(
    deg_down,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

term <- bp@result
#保存结果
write.table(term,
            file = "~/brain/SCENIC/by_group/downsrtream_analysis/GO/target_gene_of_ASD_DOWN_GOBP.csv",
            quote = F,
            sep = ",",
            row.names = F
)



#绘制TOP20
df2 <- term[1:20,]
df2$labelx=rep(0,nrow(df2))
df2$labely=seq(nrow(df2),1)
p <- ggplot(data = df2, 
            aes(x = -log10(pvalue),
                y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity",
           alpha=1,
           fill= "#87CEFA",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = df2$Description),
            size=3.5,
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  ggtitle(case)+
  scale_x_continuous(expand = c(0,0))
p
#保存
ggsave("~/brain/SCENIC/by_group/downsrtream_analysis/GO/target_gene_of_ASD_DOWN_GOBP.pdf", p,
       height=6,
       width=8)


####__________________________________________________________________________________________________________________________

TF_PTSD_DOWN<-read_csv("~/brain/SCENIC/by_group/downsrtream_analysis/VS/DEG/PTSDVSCON_down_deg.csv")$gene

PTSD_TFS<- read.table("~/brain/SCENIC/by_group/PTSD/tfs_targer.tsv", header=T, sep=",")


filtered_PTSD_TFS <- subset(PTSD_TFS, tf %in% TF_PTSD_DOWN)

target_gene_of_PTSD_DOWN<-filtered_PTSD_TFS$target_gene

#### 2.功能富集 ----

##_________________________________________________________________________________________________

case <- 'PTSD'

deg_down<-target_gene_of_PTSD_DOWN

# markers<-read.table("~/brain/差异/output/deg/ASDVSCON_deg.txt",
#                     header = T,sep=",")
# 
# deg_up <- filter(markers,
#                  markers$p_val < 0.05 & markers$avg_log2FC > 0.25)
# markers[markers$change %in% 'UP', ]$gene
# head(deg)
# #重新排序
# top.genes <-deg
# deg[order(deg$avg_log2FC, decreasing = T),]$gene
# class(deg)
# deg<-as.data.frame(deg)
bp <-
  enrichGO(
    deg_down,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

term <- bp@result
#保存结果
write.table(term,
            file = "~/brain/SCENIC/by_group/downsrtream_analysis/GO/target_gene_of_PTSD_DOWN_GOBP.csv",
            quote = F,
            sep = ",",
            row.names = F
)




#### Down 


bp <-
  enrichGO(
    deg_down,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

term <- bp@result
#保存结果
write.table(term,
            file = "~/brain/SCENIC/by_group/downsrtream_analysis/GO/target_gene_of_PTSD_DOWN_GOBP.csv",
            quote = F,
            sep = ",",
            row.names = F
)
#绘制TOP20
df2 <- term[1:20,]
df2$labelx=rep(0,nrow(df2))
df2$labely=seq(nrow(df2),1)
p <- ggplot(data = df2, 
            aes(x = -log10(pvalue),
                y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity",
           alpha=1,
           fill= "#87CEFA",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = df2$Description),
            size=3.5,
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  ggtitle(case)+
  scale_x_continuous(expand = c(0,0))
p
#保存
ggsave("~/brain/SCENIC/by_group/downsrtream_analysis/GO/target_gene_of_PTSD_DOWN_GOBP.pdf", p,
       height=6,
       width=8)

