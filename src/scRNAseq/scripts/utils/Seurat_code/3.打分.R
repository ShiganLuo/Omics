library(Seurat)
library(msigdbr) # 基因集
library(stringr) # 字符串的操作
library(ggplot2)


#### 配置与数据导入
setwd('~/精神疾病/data/388/Seurat/data/')  ## 需要修改为自己的文件夹路径
outdir_data <- paste0('./score/')
dir.create(outdir_data,recursive = T)
outdir_plots <- paste0('./picture/')
dir.create(outdir_plots,recursive = T)

seurat_obj <- readRDS('~/精神疾病/data/388/Seurat/data/merged_all_scted.rds')

table(seurat_obj$orig.ident)
#  ASD   BD  CON  MDD PTSD   SZ 
#2378 2868 2719 1551 2513 2838 

## C5-GO
category <- "C5" # C5-GO  C2-KEGG  H-Hallmark
# 从MSigDB数据库中获取基因集
Dataset <- msigdbr(species = "Homo sapiens", category = category)  # 人-Homo sapiens  鼠-Mus musculus
Dataset <- subset(Dataset, gs_subcat == 'GO:BP')

# 筛选出关注的通路对应的基因集
terms_GO <- c('Regulation of dendritic cell antigen processing and presentation',
              'cytoplasmic translation',
              'immune response-activating signaling pathway',
              'regulation of antigen processing and presentation',
              'positive regulation of leukocyte activation',
              'microglial cell activation',
              'regulation of nervous system development',
              'myeloid cell homeostasis',
              'Activation of immune response',
              'oxoglutarate metabolic process',
              'Neuron death')


pathway.name = "cytoplasmic translation"
name_path <- "cytoplasmic translation"  
terms_GO <- toupper(terms_GO) |> str_replace_all(" ", "_")
terms_GO <- paste0('GOBP_',terms_GO)
dataset <- subset(Dataset, gs_name %in% terms_GO)
table(dataset$gs_name)
#GOBP_ACTIVATION_OF_IMMUNE_RESPONSE 
    #555 
#GOBP_CYTOPLASMIC_TRANSLATION 
    #170 
#GOBP_REGULATION_OF_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION 
    #11 



# 对数据集中 基因集 进行循环：对每个基因集，提取属于这个基因集的基因，并将结果存储在一个列表中
geneSets_GO <- lapply(unique(dataset$gs_name), 
                      function(x){dataset$gene_symbol[dataset$gs_name == x]}) # 有点慢
names(geneSets_GO) <- unique(dataset$gs_name)


# ## H-Hallmark
# category <- "H"
# Dataset <- msigdbr(species = "Homo sapiens", category = category)
# table(Dataset$gs_name)
# dataset <- Dataset
# geneSets_H <- lapply(unique(dataset$gs_name), 
#                      function(x){dataset$gene_symbol[dataset$gs_name == x]})
# names(geneSets_H) <- unique(dataset$gs_name)
# 

# AUCell ------------------------------------------------------------------
library(AUCell)
### GO通路打分 ----
# 排序：对每个细胞,基于其所有基因的表达进行排序,从高到低
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data)
# 计算每个细胞中每个基因集的AUC. 为应对大数据及去噪，只有每个细胞表达量最高的前5%的基因参与计算
cells_AUC <- AUCell_calcAUC(geneSets_GO, cells_rankings, aucMaxRank=ceiling(0.05 * nrow(cells_rankings)))  
# 提取结果
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
# 重命名列
colnames(AUCell_socre) <- colnames(AUCell_socre)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
score_GO <- AUCell_socre
#saveRDS(score_GO, paste0(outdir_data, "AUCell_score_GO_pathway.rds"))



# ### H-Hallmark通路打分 ----
# cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 
# cells_AUC <- AUCell_calcAUC(geneSets_H, cells_rankings, aucMaxRank=ceiling(0.05 * nrow(cells_rankings))) ##
# AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
# colnames(AUCell_socre) <- colnames(AUCell_socre)|>
#   str_replace_all("HALLMARK_", "") |>
#   str_replace_all("_", " ") |>
#   str_to_sentence()
# # saveRDS(AUCell_socre, paste0(outdir_data,"AUCell_score_Hallmark_pathway.rds"))
# 
# terms_H <- c(
#   'tgf beta signaling',
#   'epithelial mesenchymal transition'
# )
# terms_H <- str_to_sentence(terms_H) # 句首字母大写
# score_H <- AUCell_socre[ ,colnames(AUCell_socre) %in% terms_H]
# colnames(score_H)
# # [1] "Epithelial mesenchymal transition" "Tgf beta signaling"





### 相关性准备 ----
## 挑选关注的基因，整合其表达与打分结果，为后续相关性分析做准备
#score_selected <- cbind(score_GO, score_H)
seurat_obj_AUCell <- seurat_obj
# a <- t(as.data.frame(seurat_obj@assays[["RNA"]]@data[c("MRC1","SIGLEC10"),]))
# score_all_selected <- cbind(score_selected, a)

seurat_obj_AUCell@meta.data <- cbind(seurat_obj@meta.data, score_GO)
head(seurat_obj_AUCell@meta.data)

#saveRDS(seurat_obj_AUCell, paste0(outdir_data,"AUCell_score_all.rds"))


### 打分结果可视化 ----
#name_path <- 'regulation of dendritic cell antigen processing and presentation'
library(ggpubr)

seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/merged/merged_all_scted.rds")
# seurat_obj <- subset(seurat_obj,disease %in% c("CHD","AMI","ICM","CS","COVID-19","normal"))
C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)
terms <- c("GOBP_ACTIVATION_OF_IMMUNE_RESPONSE"
  #"GOBP_NEURON_DEATH"
  #"GOBP_MICROGLIAL_CELL_ACTIVATION"
  ##"GOBP_CYTOPLASMIC_TRANSLATION"
  #"GOBP_REGULATION_OF_ENDOCYTOSIS"
  #"GOBP_INTERLEUKIN-8_PRODUCTION"
)

selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% terms) #%>% pull(gene_symbol)
selected_gene_sets
cells_rankings <- AUCell_buildRankings(as.matrix(seurat_obj@assays$RNA@data))
cells_rankings
geneSets<-lapply(unique(selected_gene_sets$gs_name),
                 function(x){selected_gene_sets$gene_symbol[selected_gene_sets$gs_name==x]})
#View(geneSets)
names(geneSets) <- unique(selected_gene_sets$gs_name)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC
#saveRDS(cells_AUC,"./2/cells_AUC.rds")
for( i in terms){
  aucs <- as.numeric(getAUC(cells_AUC)[i, ])
  #将AUC的结果添加到seurat.object的meta.data中，并画图
  #seurat_obj$AUC <- aucs
  seurat_obj$i <- aucs
  df <- seurat_obj@meta.data
  colnames(df)[7] <- paste0(i)
  # p <- VlnPlot(seurat_obj,features = paste0(i), 
  #         pt.size = 0,group.by = "disease")
  colnames(df)[7] <- gsub("GOBP_", "", colnames(df)[7])
  colnames(df)[7] <- gsub("_", " ", colnames(df)[7])
  colnames(df)[7] <- tolower(colnames(df)[7])
  # 画图
  Data <- df
  mean.score <- df
  mean.score <- mean(mean.score[,7])
  p <- ggplot(Data, aes(x=orig.ident, y=Data[,7],fill=orig.ident)) + 
    geom_violin(trim=FALSE,color="white") + #绘制小提琴图
    geom_boxplot(width=0.2,position=position_dodge(0.9))+
    # scale_fill_manual(values = ggsci::pal_aaas()(10))+ #设置填充的颜色
    scale_fill_manual(values = c("#9370DB", "#96C3D8", "#5F9BBE","#4dbbd5ff","#7FC97F","#F19294"))+ #,"#4A9D47","darkgoldenrod1","darkorange","#F19294","#b5a6d3"))+
    theme_bw()+ #背景变为白色
    geom_hline(yintercept = mean.score,
               color = "black",
               linetype = "dashed") +
    theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(size=8,), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain family="Times", face="plain"
          axis.title.y=element_text(size = 8), #设置y轴标题的字体属性
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                   size=8), # face="italic", family="Arial", 
          legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                    size=8),
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+  #不显示网格线
    ylab("")+xlab("")+ #设置x轴和y轴的标题+
    ggtitle(colnames(df)[7])
  # outdir <- "~/20240415_heart_SMC/result/4-deg/v1/4-gene/diseaseDEG/5disease/intersect_termgene_commonupDEG_score/2/"
  # ggsave(paste0( outdir, i,".pdf"),p,height = 4,width =6)
}
p
outdir <- "~/精神疾病/data/388/Seurat/picture/score/"
ggsave(paste0( outdir, i,".pdf"),p,height = 4,width =6)

##___


###########作图时将CON放在第一位
seurat_obj$group <- seurat_obj$orig.ident

# 重新设置'group'的水平顺序  
new_levels <- c("CON", "ASD", "BD", "MDD","PTSD","SZ") 
seurat_obj$group <- factor(seurat_obj$group, levels = new_levels)  

print(levels(seurat_obj$group))
table(seurat_obj$group)

file_path <- file.path(path.expand("~/精神疾病/data/388/Seurat/data"), "SCTed_data.rds")  

# 使用saveRDS函数保存对象  
saveRDS(seurat_obj, file_path)







