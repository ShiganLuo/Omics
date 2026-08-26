#cd ~/LM/project2304/src/

####配置----
library(Seurat)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggplot2)
library(UpSetR)
library(qs)
library(clusterProfiler)
library(org.Hs.eg.db)

# #setwd('/home/chencx/project/PH/src')
# path <- '~/scrna/PH/results/F2/'
# setwd(path)
 outdir2 <- '~/精神疾病/data/388/Seurat/picture'
# dir.create(outdir2, recursive = T)

#DefaultAssay(seurat_obj) <- 'RNA'

ASD_up <- read_csv("DEG/ASD_up_deg.csv")
BD_up <- read_csv("DEG/BD_up_deg.csv")
MDD_up <- read_csv("DEG/MDD_up_deg.csv")
PTSD_up <- read_csv("DEG/PTSD_up_deg.csv")
SZ_up <- read_csv("DEG/SZ_up_deg.csv")

DEG_split <- list(ASD_down, BD_down, MDD_down, PTSD_down, SZ_down)  

saveRDS(DEG_split,"~/精神疾病/data/388/Seurat/data/DEG/DEG_split.rds")
# 
# color_ct <- c("#FB8072", "#80B1D3", "#FDB462", "#8DD3C7",
#                     "#FCCDE5", "#BC80BD", "#CCEBC5", "#377EB8", "#1B9E77") #,'#B3DE69'
# feature <- c("ACM","DCM","HCM","ICM",
#              "AMI","HF","CHD","CS","COV") #,'Normal'

#### upset图+Barplot ----
## 1全部差异基因 ----
# outdir <- './files/all/'
# dir.create(outdir, recursive = T)

#### 1.1 upset ----
# color_ct <- c("#FCCDE5", "#377EB8", "#8DD3C7", "#FDB462",
#               "#1B9E77", "#80B1D3", "#FB8072", "#CCEBC5", "#BC80BD")
# DEG_split <- df_list
# DEG_split <- split(DEG, DEG$celltype)
for(i in 1:5){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
upsetData=fromList(DEG_split)
#case = '1125_disease_all_'


#绘图
pdf(file=file.path(outdir2, paste0(case, 'upset.pdf')),onefile = FALSE,width=8,height=5)
upset(upsetData,
      nsets = length(DEG_split),               #展示多少个数据
      nintersects = 5,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "yes",                   #柱状图上方是否显示数值
      number.angles = 0,                     #字体角度
      point.size = 3,                         #点的大小
      matrix.color="#b0b9b8",                     #交集点颜色
      line.size = 0.8,                        #线条粗细
      mainbar.y.label = "Gene Intersections", #y轴文字
      sets.x.label = "Set Size",              #x柱状图文字
      sets.bar.color = "#b0b9b8", #左下方条形图颜色
      main.bar.color ="black"#, #右上方条形图颜色
      # queries = list(
      #   list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","HF","CHD","CS","COV"), color= "#1fab89", active=T),
      #   #active：如果为F，则在每一列上方显示一个三角
      #   list(query=intersects, params=list("ACM","DCM","HCM","ICM","CHD","CS","COV"), color= "#1fab89", active=T),
      #   list(query=intersects, params=list("AMI"), color="#eeb401", active=T),
      #   list(query=intersects, params=list('CHD'), color="#633372", active=T),
      #   list(query=intersects, params=list('CS'), color="#D55E00", active=T),
      #   list(query=intersects, params=list('COV'), color="#008080", active=T),
      #   list(query=intersects, params=list('ACM','DCM'), color="#c82d31", active=T),
      #   list(query=intersects, params=list('ICM','CS','AMI'), color="#0e2c82", active=T),
      #   list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CHD","CS","COV"), color="#1fab89", active=T),
      #   list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CS","COV"), color= "#1fab89", active=T),
      #   list(query=intersects, params=list("ACM","DCM","HF","ICM", "AMI","CS"), color= "#c82d31", active=T)
      #   )
      ) 
dev.off()


#保存交集文件
intersectGenes=Reduce(intersect,DEG_split)
write.table(file=file.path(outdir, paste0(case, 'intergenes.txt')),
            intersectGenes,sep=",",quote=F,col.names=T,row.names=F)

#### 1.2 特异基因 ----
type <- 'ICM'

unique_geneset <- setdiff(DEG_split[[9]],
                          union(DEG_split[[1]],
                                c(DEG_split[[2]],DEG_split[[3]],DEG_split[[4]],DEG_split[[5]],
                                  DEG_split[[6]],DEG_split[[7]],DEG_split[[8]])))

write.table(file=file.path(outdir, paste0(case, type, '_uniquegene.txt')),
            unique_geneset,sep=",",quote=F,col.names=T,row.names=F)

##### 2 上调差异基因 ----
outdir <- './files/up/'
dir.create(outdir, recursive = T)

DEG <- readRDS("/home/chencx/project/PH/src/results/DEG/disease/files/markers.BasedOncondition.rds")|>
  filter(avg_log2FC > 0.25 & p_val < 0.05)
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
upsetData=fromList(DEG_split)
case = '1125_disease_UP_'

#绘图
color_ct <- c("#FCCDE5", "#377EB8", "#8DD3C7", "#80B1D3",
              "#FDB462", "#FB8072", "#1B9E77", "#BC80BD", "#CCEBC5")

pdf(file=file.path(outdir2, paste0(case, 'upset.pdf')),onefile = FALSE,width=8,height=5)
upset(upsetData,
      nsets = length(DEG_split),               #展示多少个数据
      nintersects = 24,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "yes",                   #柱状图上方是否显示数值
      number.angles = 0,                     #字体角度
      point.size = 3,                         #点的大小
      matrix.color="#b0b9b8",                     #交集点颜色
      line.size = 0.8,                        #线条粗细
      mainbar.y.label = "Intersections of up-regulated genes", #y轴文字
      sets.x.label = "Set Size",              #x柱状图文字
      sets.bar.color = color_ct,#"#b0b9b8", #左下方条形图颜色
      main.bar.color ="black", #右上方条形图颜色
      queries = list(
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","HF","CHD","CS","COV"), color= "#1fab89", active=T),
        #active：如果为F，则在每一列上方显示一个三角
        list(query=intersects, params=list("ACM","DCM","HCM","ICM","HF","CS","COV"), color= "#1fab89", active=T),
        list(query=intersects, params=list("AMI"), color="#eeb401", active=T),
        list(query=intersects, params=list('HF'), color="#633372", active=T),
        list(query=intersects, params=list('CS'), color="#D55E00", active=T),
        list(query=intersects, params=list('COV'), color="#008080", active=T),
        list(query=intersects, params=list('ACM','DCM'), color="#c82d31", active=T),
        list(query=intersects, params=list('ICM','CS','AMI'), color="#0e2c82", active=T),
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CS","COV"), color= "#1fab89", active=T)
        )
      )  
dev.off()

#保存交集文件
intersectGenes=Reduce(intersect,DEG_split)
write.table(file=file.path(outdir, paste0(case, 'intergenes.txt')),
            intersectGenes,sep=",",quote=F,col.names=T,row.names=F)

#### 2.4特异的基因集合 ----
type <- 'MI' # 5
type <- 'TOF' # 6
type <- 'HF' # 3
# unique_geneset <- setdiff(DEG_split[[type]],
#                           union(DEG_split[[2]],c(DEG_split[[3]],DEG_split[[4]],DEG_split[[5]],DEG_split[[6]])))

type <- 'HCM DCM' # 9
# unique_geneset <- setdiff(union(DEG_split[[1]],DEG_split[[2]]), 
#                           union(DEG_split[[3]],c(DEG_split[[4]],DEG_split[[5]],DEG_split[[6]])))

#### 功能富集 
# top.genes <- unique_geneset
# bp <-
#   enrichGO(
#     top.genes,
#     OrgDb = org.Hs.eg.db,
#     keyType = 'SYMBOL',
#     ont = "BP",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.2)
# 
# term <- bp@result
# term$celltype <- type
# 
# write.table(
#   term,
#   file = file.path(outdir, paste0(case,type,'GO.txt')),
#   row.names = F,
#   quote = F,
#   sep = "\t")
# saveRDS(term, file.path(outdir, paste0(case,type,'GO.rds')))
#  
#### MI
term <- readRDS(file.path(outdir, paste0(case,type,'GO.rds')))
terms <- c(
  'cardiac muscle tissue development',
  'amino acid activation',
  'heart process',
  'heart contraction',
  'cardiac chamber development',
  'regulation of blood circulation',
  'actin filament-based movement',
  'cardiac conduction',
  'cardiac ventricle development',
  'heart morphogenesis'
)


#### TOF
term <- readRDS(file.path(outdir, paste0(case,type,'GO.rds')))
terms <- c(
  'organic anion transport',
  'carboxylic acid transport',
  'canonical Wnt signaling pathway',
  'L-arginine transmembrane transport',
  'regulation of canonical Wnt signaling pathway',
  'basic amino acid transmembrane transport',
  'aromatic amino acid transport',
  'regulation of Wnt signaling pathway',
  'organic acid transport',
  'azole transmembrane transport'
)
#### HF
term <- readRDS(file.path(outdir, paste0(case,type,'GO.rds')))
terms <- c(
  'regulation of blood pressure',
  'regulation of systemic arterial blood pressure',
  'intrinsic apoptotic signaling pathway in response to DNA damage',
  'positive regulation of endothelial cell proliferation',
  'mitotic DNA damage checkpoint signaling',
  'endothelial cell proliferation',
  'positive regulation of protein kinase B signaling',
  'mitotic G1 DNA damage checkpoint signaling',
  'positive regulation of angiogenesis',
  'positive regulation of vasculature development'
)


#### DCM, HCM
term <- readRDS(file.path(outdir, paste0(case,type,'GO.rds')))
terms <- c(
  'nucleocytoplasmic transport',
  'small nucleolar ribonucleoprotein complex assembly',
  'regulation of protein ubiquitination',
  'response to interferon-alpha',
  'histone H2A monoubiquitination',
  'activin receptor signaling pathway',
  'histone ubiquitination',
  'ventricular septum development',
  'vasodilation',
  'vascular process in circulatory system',
  'integrin-mediated signaling pathway'
)

colnames(term)
# "ID" "Description" "GeneRatio"   "BgRatio"   "pvalue"   "p.adjust"   "qvalue"   "geneID"   "Count" 
df <- term[term$Description %in% terms, ]
df$pvalue <- -log10(df$pvalue)
#df$pvalue[df$pvalue >400] = 400 
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)

ggplot(data = df,
       aes(Count, reorder(Description,Count), fill = pvalue)) +
  geom_bar(stat="identity",
           #fill='red',
           alpha=1,
           width = 0.8) + 
  # scale_fill_gradient(low = '#fff5a5', high = "#eeb401") + #MI
  # scale_fill_gradient(low = '#fdc7ff', high = "#633372") + #TOF
  # scale_fill_gradient(low = '#ffdfdf', high = "#CC79A7") + #HF
  scale_fill_gradient(low = '#c7f5fe', high = "#0072B2") + #HCM DCM
  geom_text(aes(x = labelx,
                y = Description,
                label = Description),
            size=3.5, 
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),plot.title = element_text(hjust = 0.5,vjust = -1),
        legend.title = element_text(colour = 'black',size = 10,vjust = 3),#lenged的大小颜色位置
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  labs(title = type, x="Count",fill = "-log10(pvalue)")+#设置legend 文字
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))#条目太长换行

ggsave(file.path(outdir2,paste0(case, type, 'Barplot.pdf')),
       ggplot2::last_plot(),#最后一张图
       height=3,
       width=5)



#### 3 下调差异基因 ----
outdir <- './files/down/'
dir.create(outdir, recursive = T)

#### 3.1 upset图 ----
DEG <- readRDS("/home/chencx/project/PH/src/results/DEG/disease/files/markers.BasedOncondition.rds")|>
  filter(avg_log2FC < -0.25 & p_val < 0.05)
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
upsetData=fromList(DEG_split)
case = '1125_disease_DOWN_'


#绘图
color_ct <- c("#CCEBC5", "#1B9E77", "#FDB462", "#377EB8",
              "#8DD3C7", "#FCCDE5", "#8DD3C7", "#FB8072", "#BC80BD")

pdf(file=file.path(outdir2, paste0(case, 'upset.pdf')),onefile = FALSE,width=8,height=5)
upset(upsetData,
      nsets = length(DEG_split),               #展示多少个数据
      nintersects = 24,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "yes",                   #柱状图上方是否显示数值
      number.angles = 0,                     #字体角度
      point.size = 3,                         #点的大小
      matrix.color="#b0b9b8",                     #交集点颜色
      line.size = 0.8,                        #线条粗细
      mainbar.y.label = "Intersections of down-regulated genes", #y轴文字
      sets.x.label = "Set Size",              #x柱状图文字
      sets.bar.color = color_ct,#"#b0b9b8", #左下方条形图颜色
      main.bar.color ="black", #右上方条形图颜色
      queries = list(
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","HF","CHD","CS","COV"), color= "#1fab89", active=T),
        #active：如果为F，则在每一列上方显示一个三角
        list(query=intersects, params=list("ACM","DCM","HCM","ICM","CHD","CS","COV"), color= "#1fab89", active=T),
        list(query=intersects, params=list('CHD'), color="#633372", active=T),
        list(query=intersects, params=list('COV'), color="#008080", active=T),
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CHD","COV"), color= "#1fab89", active=T)
        )
      ) 

dev.off()

#保存交集文件
intersectGenes=Reduce(intersect,DEG_split)
write.table(file=file.path(outdir, paste0(case, 'intergenes.txt')),
            intersectGenes,sep=",",quote=F,col.names=T,row.names=F)


#### 3.2功能富集 ----
library(clusterProfiler)
library(org.Hs.eg.db)
top.genes <- intersectGenes
#功能富集  #生物过程：BH; 细胞分组：CC;分子功能：MF
bp <-
  enrichGO(
    top.genes,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

term <- bp@result
term$celltype <- 'down-regulated terms'

write.table(
  term,
  file = file.path(outdir, paste0(case, 'GO.txt')),
  row.names = F,
  quote = F,
  sep = "\t")
saveRDS(term, file.path(outdir, paste0(case, 'GO.rds')))

#### 3.3柱状图 ----
terms <- c('oxidative phosphorylation',
           'cellular respiration',
           'mitochondrial ATP synthesis coupled electron transport',
           'proton motive force-driven ATP synthesis',
           'electron transport chain',
           'purine ribonucleoside triphosphate biosynthetic process',
           'ATP metabolic process',
           'nucleoside triphosphate metabolic process',
           'viral process',
           'intrinsic apoptotic signaling pathway',
           'proteasome-mediated ubiquitin-dependent protein catabolic process')
colnames(term)
# "ID" "Description" "GeneRatio"   "BgRatio"   "pvalue"   "p.adjust"   "qvalue"   "geneID"   "Count" 
df <- term[term$Description %in% terms, ]
df$pvalue <- -log10(df$pvalue)
#df$pvalue[df$pvalue >400] = 400 
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)

ggplot(data = df,
       aes(Count, reorder(Description,Count), fill = pvalue)) +
  geom_bar(stat="identity",
           #fill='red',
           alpha=1,
           width = 0.8) + 
  scale_fill_gradient(low = "#a3daff", high = "#4ea1d3") +
  geom_text(aes(x = labelx,
                y = Description,
                label = Description),
            size=3.5, 
            hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -1),
        legend.title = element_text(colour = 'black',size = 10,vjust = 3),#lenged的大小颜色位置
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  labs(fill = "-log10(pvalue)")+#设置legend 文字
  xlab("Count")+#主标题
  ggtitle("Down-regulated pathway")+#+scale_x_continuous(expand = c(0,0))
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))#条目太长换行


ggsave(file.path(outdir2,paste0(case, 'Barplot.pdf')),
       ggplot2::last_plot(),#最后一张图
       height=4,
       width=6)


