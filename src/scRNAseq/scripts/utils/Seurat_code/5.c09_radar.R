library(ggradar)

pathways <- readRDS('/home/chencx/project/PH/angiogenesis/function/results/files/metascape_CAPFs.rds')
groups <-  as.character(unique(pathways$group))

#### c08 ----
sub_obj <- subset(seurat_obj0, cell_type %in% 'Cap.c08.TNFRSF4')  #,'Cap.c09.RGS6'

df_c08 <- data.frame(rownames(sub_obj@meta.data))
for (i in groups) {
  sub_df <- pathways[pathways$group == i,]
  for (k in 1:3) {
    s <- sub_df[k,8]
    assign(paste0('str',k), strsplit(s, '\\|')[[1]])
  }
  gene_list <- Reduce(union, list(str1,str2,str3))
  
  cells_ranking <- AUCell_buildRankings(sub_obj@assays$RNA@counts)
  cells_AUC <- AUCell_calcAUC(gene_list, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
  AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
  df_c08 <- cbind(df_c08,AUCell_socre)
}

df_c08 <- df_c08[,-1]
colnames(df_c08) <- groups
c08_score <- colMeans(df_c08)
c08_score <- (c08_score-min(c08_score))/(max(c08_score)-min(c08_score))
c08_score <- as.data.frame(c08_score)

ggradar(c08_score, background.circle.transparency = 0,
        group.colours = '#E64B35',
        values.radar = c('0','0.5','1'),
        plot.extent.x.sf = 1.2,
        base.size = 3,axis.label.size = 4,grid.label.size = 4,
        group.point.size = 2,
        group.line.width = 1,
        # plot.title = 'none',
        legend.position = 'top',
        legend.text.size = 10
)+
  theme(plot.title = element_text(hjust = 0.2,size = 12))
ggsave(paste0(outdir2,'radar_pathways_c08.pdf'),last_plot(),width = 6,height = 4)


#### c09 ----
sub_obj <- subset(seurat_obj0, cell_type %in% 'Cap.c09.RGS6')

df_c09 <- data.frame(rownames(sub_obj@meta.data))
for (i in groups) {
  sub_df <- pathways[pathways$group == i,]
  for (k in 1:3) {
    s <- sub_df[k,8]
    assign(paste0('str',k), strsplit(s, '\\|')[[1]])
  }
  gene_list <- Reduce(union, list(str1,str2,str3))
  
  cells_ranking <- AUCell_buildRankings(sub_obj@assays$RNA@counts)
  cells_AUC <- AUCell_calcAUC(gene_list, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
  AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
  df_c09 <- cbind(df_c09,AUCell_socre)
}

df_c09 <- df_c09[,-1]
colnames(df_c09) <- groups
c09_score <- colMeans(df_c09)
c09_score <- (c09_score-min(c09_score))/(max(c09_score)-min(c09_score))
c09_score <- t(c09_score)
c09_score <- as.data.frame(c09_score)

ggradar(c09_score, background.circle.transparency = 0,
        group.colours = '#4DBBD5',
        values.radar = c('0','0.5','1'),
        plot.extent.x.sf = 1.2,
        base.size = 3,axis.label.size = 4,grid.label.size = 4,
        group.point.size = 2,
        group.line.width = 1,
        # plot.title = 'none',
        legend.position = 'top',
        legend.text.size = 10
)+
  theme(plot.title = element_text(hjust = 0.2,size = 12))
ggsave(paste0(outdir2,'radar_pathways_c09.pdf'),last_plot(),width = 6,height = 4)

# adata <- data
# data <- as.data.frame(c08_score)
# data$max <- 0.35
# data$min <- 0.2
# data <- t(data)
# data <- as.data.frame(data[c(2,3,1),])
# radarchart(data,
#            pcol = "#E64B35",
#            pfcol =  scales::alpha("#E64B35", 0.5),
#            plty = "solid",
#            cglty = "solid",
#            cglcol = "black",
#            cglwd =0.5)
# 
# data <- as.data.frame(c09_score)
# data$max <- 0.35
# data$min <- 0.2
# data <- t(data)
# data <- as.data.frame(data[c(2,3,1),])
# radarchart(data,
#            pcol = "#4DBBD5",
#            pfcol =  scales::alpha("#4DBBD5", 0.5),
#            plty = "solid",
#            cglty = "solid",
#            cglcol = "black",
#            cglwd =0.5)

#### 通路打分 ----
func <- as.character(unique(pathways$Description))

##
obj_c08 <- subset(seurat_obj0, cell_type == 'Cap.c08.TNFRSF4')
df_c08 <- data.frame(rownames(obj_c08@meta.data))
# i <- "cell migration involved in sprouting angiogenesis"
for (i in func) {
  sub_df <- pathways[pathways$Description == i,]
  s <- sub_df[,8]
  gene_list <- strsplit(s, '\\|')[[1]]
  
  cells_ranking <- AUCell_buildRankings(obj_c08@assays$RNA@counts)
  cells_AUC <- AUCell_calcAUC(gene_list, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
  AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
  df_c08 <- cbind(df_c08,AUCell_socre)
}
df_c08 <- df_c08[,-1]
colnames(df_c08) <- func
c08_score <- colMeans(df_c08) |> as.data.frame()

##
obj_c09 <- subset(seurat_obj0, cell_type == 'Cap.c09.RGS6')
df_c09 <- data.frame(rownames(obj_c09@meta.data))
# i <- "cell migration involved in sprouting angiogenesis"
for (i in func) {
  sub_df <- pathways[pathways$Description == i,]
  s <- sub_df[,8]
  gene_list <- strsplit(s, '\\|')[[1]]
  
  cells_ranking <- AUCell_buildRankings(obj_c09@assays$RNA@counts)
  cells_AUC <- AUCell_calcAUC(gene_list, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
  AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
  df_c09 <- cbind(df_c09,AUCell_socre)
}
df_c09 <- df_c09[,-1]
colnames(df_c09) <- func
c09_score <- colMeans(df_c09) |> as.data.frame()

##
score <- cbind(c08_score,c09_score)
colnames(score) <- c('c08','c09')
write.csv(score, paste0(outdir,'c08_c09_score.csv'))

library(pheatmap)
p <- pheatmap(score,cluster_cols = F, 
         cluster_rows = F, scale="row",
         fontsize_number = 20, border="white", #设置边框为白色
         treeheight_row = 10, treeheight_col = 10,
         shape = "circle",cellwidth = 10,cellheight =10,#设置热图方块宽度和高度
         # color=colorRampPalette(c("darkkhaki",'white', "darkgreen"))(10000),#设置颜色梯度
         # color=colorRampPalette(c("darkblue",'yellow', "darkred"))(10000),#设置颜色梯度
         # legend = F,#main=s_m6a, # 显示图例 # 设置图形标题
         fontsize_row = 10,fontsize_col = 10 ) # 分别设置横向和纵向字体大小

