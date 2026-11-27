library(Seurat)
library(msigdbr)
library(stringr)
library(ggplot2)


setwd('~/精神疾病/data/388/Seurat/data/')  ## 需要修改为自己的文件夹路径
outdir_data <- paste0('./score/')
dir.create(outdir_data,recursive = T)
outdir_plots <- paste0('./picture/')
dir.create(outdir_plots,recursive = T)

seurat_obj <- readRDS('~/精神疾病/data/388/Seurat/data/merged_all_scted.rds')

table(seurat_obj$orig.ident)

category <- "C5" # C5-GO  C2-KEGG  H-Hallmark
Dataset <- msigdbr(species = "Homo sapiens", category = category)  # 人-Homo sapiens  鼠-Mus musculus
Dataset <- subset(Dataset, gs_subcat == 'GO:BP')

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



geneSets_GO <- lapply(unique(dataset$gs_name), 
                      function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets_GO) <- unique(dataset$gs_name)



library(AUCell)
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data)
cells_AUC <- AUCell_calcAUC(geneSets_GO, cells_rankings, aucMaxRank=ceiling(0.05 * nrow(cells_rankings)))  
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
colnames(AUCell_socre) <- colnames(AUCell_socre)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
score_GO <- AUCell_socre








seurat_obj_AUCell <- seurat_obj

seurat_obj_AUCell@meta.data <- cbind(seurat_obj@meta.data, score_GO)
head(seurat_obj_AUCell@meta.data)



library(ggpubr)

seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/merged/merged_all_scted.rds")
C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)
terms <- c("GOBP_ACTIVATION_OF_IMMUNE_RESPONSE"
)

selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% terms)
selected_gene_sets
cells_rankings <- AUCell_buildRankings(as.matrix(seurat_obj@assays$RNA@data))
cells_rankings
geneSets<-lapply(unique(selected_gene_sets$gs_name),
                 function(x){selected_gene_sets$gene_symbol[selected_gene_sets$gs_name==x]})
names(geneSets) <- unique(selected_gene_sets$gs_name)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC
for( i in terms){
  aucs <- as.numeric(getAUC(cells_AUC)[i, ])
  seurat_obj$i <- aucs
  df <- seurat_obj@meta.data
  colnames(df)[7] <- paste0(i)
  colnames(df)[7] <- gsub("GOBP_", "", colnames(df)[7])
  colnames(df)[7] <- gsub("_", " ", colnames(df)[7])
  colnames(df)[7] <- tolower(colnames(df)[7])
  Data <- df
  mean.score <- df
  mean.score <- mean(mean.score[,7])
  p <- ggplot(Data, aes(x=orig.ident, y=Data[,7],fill=orig.ident)) + 
    geom_violin(trim=FALSE,color="white") + #绘制小提琴图
    geom_boxplot(width=0.2,position=position_dodge(0.9))+
    scale_fill_manual(values = c("#9370DB", "#96C3D8", "#5F9BBE","#4dbbd5ff","#7FC97F","#F19294"))+ #,"#4A9D47","darkgoldenrod1","darkorange","#F19294","#b5a6d3"))+
    theme_bw()+
    geom_hline(yintercept = mean.score,
               color = "black",
               linetype = "dashed") +
    theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(size=8,),
          axis.title.y=element_text(size = 8),
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                   size=8),
          legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                    size=8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ylab("")+xlab("")+ #设置x轴和y轴的标题+
    ggtitle(colnames(df)[7])
}
p
outdir <- "~/精神疾病/data/388/Seurat/picture/score/"
ggsave(paste0( outdir, i,".pdf"),p,height = 4,width =6)



seurat_obj$group <- seurat_obj$orig.ident

new_levels <- c("CON", "ASD", "BD", "MDD","PTSD","SZ") 
seurat_obj$group <- factor(seurat_obj$group, levels = new_levels)  

print(levels(seurat_obj$group))
table(seurat_obj$group)

file_path <- file.path(path.expand("~/精神疾病/data/388/Seurat/data"), "SCTed_data.rds")  

saveRDS(seurat_obj, file_path)







