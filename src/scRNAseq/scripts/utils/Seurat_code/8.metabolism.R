#### 加载R包 ----
library(Seurat)
library(SeuratData)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(magrittr)
library(scMetabolism)
library(ggplot2)
library(rsvd)

#### 配置 ----
path <- '/home/maxr/精神疾病/data/388/Seurat/data/metabolism'
dir.create(path, recursive = T)
setwd(path)

# load dataset
ifnb <-
  readRDS("~/精神疾病/data/388/Seurat/data/SCTed_new.rds")
head(ifnb@meta.data)
#              orig.ident nCount_RNA nFeature_RNA group percent.mt nCount_SCT nFeature_SCT
# CON_CON1_1        CON      14969         5223   CON          0       1547         1015
# CON_CON1_2        CON      13853         5063   CON          0       1831         1185
# CON_CON1_3        CON       5634         3104   CON          0       1501         1043
# CON_CON1_4        CON       5134         3027   CON          0       1320          952
# CON_CON1_5        CON       4357         2435   CON          0       1226          788
# CON_CON1_6        CON       4065         2447   CON          0       1789         1378

table(ifnb@meta.data$group)
# ASD   BD  CON  MDD PTSD   SZ 
# 2378 2868 2719 1551 2513 2838

ifnb.list <- SplitObject(ifnb, split.by = "group")

ifnb.list = NormalizeData(object = ifnb.list, assay = "RNA")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca", k.anchor = 20)
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

# 聚类
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

DefaultAssay(immune.combined) <- "RNA"
Idents(immune.combined) <- "group"

# UMAP可视化

# p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "group")
# p2 <- DimPlot(immune.combined, reduction = "umap", group.by = "seurat_annotations",label = TRUE,repel = TRUE)

pdf("immune.combined.UMAP.pdf", height = 5,width = 12)
p1 
dev.off()

# 保存数据
saveRDS(immune.combined,"immune.combined.rds")

#### 
seurat_obj <- readRDS("immune.combined.rds")
Idents(seurat_obj) <- seurat_obj$group

###打分计算
countexp.Seurat <- sc.metabolism.Seurat(obj = seurat_obj, 
                                        method = "AUCell", # supports VISION, AUCell, ssgsea, and gsva, which VISION is the default method.
                                        imputation = F, ncores = 2, 
                                        metabolism.type = "REACTOME") # supports KEGG and REACTOME, where KEGG contains 85 metabolism pathways and REACTOME contains 82 metabolism pathways.
saveRDS(countexp.Seurat,"countexp.Seurat1.rds")



####-------------------------- UMAP展示通路分布-------------------------------####
countexp.Seurat <- readRDS("countexp.Seurat1.rds")

signature_exp <- countexp.Seurat@assays$METABOLISM$score
signature_meta <- countexp.Seurat@meta.data[,c("group")]

DimPlot.metabolism(obj = countexp.Seurat, pathway = "Integration of energy metabolism", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Metabolism of RNA", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)


####-------------------------- 选择通路，绘制气泡图-------------------------------####
input.pathway<-c("Metabolism of nitric oxide NOS3 activation and regulation",
                 "Selenoamino acid metabolism",
                 "Regulation of lipid metabolism by pparalpha",
                 "Foxo mediated transcription of oxidative stress metabolic and neuronal genes", 
                 "Sialic acid metabolism",
                 "Metabolism of amino acids and derivatives",
                 "Regulation of glycolysis by fructose 2 6 bisphosphate metabolism",
                 "Carnitine metabolism",
                 "TP53 regulates metabolic genes",
                 "Metabolism of steroid hormones"
                 )
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "group", norm = "y") + 
  scale_x_discrete(limits = c("CON", "ASD", "BD", "MDD", "PTSD", "SZ"))
p1 <- DotPlot.metabolism

####-------单独通路量化（箱形图）
input.pathway<-c("Metabolism of steroid hormones")
BoxPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "group", ncol = 2)+ 
  scale_x_discrete(limits = c("CON", "ASD", "BD", "MDD", "PTSD", "SZ"))



#### --------------------------------热图展示与control的差异----------------------------####
signature_df <- data.frame(signature_meta, t(signature_exp))
colnames(signature_df)[1:nrow(signature_exp)+1] <- rownames(signature_exp)

#提取各组数据
signature_df_CON <- subset(signature_df, signature_meta%in%"CON")[0:nrow(signature_exp)+1]
signature_df_ASD <- subset(signature_df, signature_meta%in%"ASD")[0:nrow(signature_exp)+1]
signature_df_BD <- subset(signature_df, signature_meta%in%"BD")[0:nrow(signature_exp)+1]
signature_df_MDD <- subset(signature_df, signature_meta%in%"MDD")[0:nrow(signature_exp)+1]
signature_df_PTSD <- subset(signature_df, signature_meta%in%"PTSD")[0:nrow(signature_exp)+1]
signature_df_SZ <- subset(signature_df, signature_meta%in%"SZ")[0:nrow(signature_exp)+1]

signature_df_CON <- signature_df_CON[, -1]
signature_df_ASD <- signature_df_ASD[, -1]
signature_df_BD <- signature_df_BD[, -1]
signature_df_MDD <- signature_df_MDD[, -1]
signature_df_PTSD <- signature_df_PTSD[, -1]
signature_df_SZ <- signature_df_SZ[, -1]

###对比p值
df_pvalue <- c()
for (i in 1:ncol(signature_df_SZ)) {
  pvalue <- wilcox.test(signature_df_CON[,i], signature_df_SZ[,i]) %>% .$p.value
  df <- data.frame(Pathway = colnames(signature_df_SZ)[i], p_val = pvalue)
  df_pvalue <- rbind(df_pvalue,df)
}

#
p.adj <- p.adjust(df_pvalue$p_val,method = "BH") 
df_pvalue$p.adj <- p.adj

ctrl_df <- apply(signature_df_CON, 2, median) %>% as.data.frame() 
colnames(ctrl_df) <- "CON"
ctrl_df$Pathway <- as.vector(rownames(ctrl_df))

stim_df <- apply(signature_df_SZ, 2, median) %>% as.data.frame()
colnames(stim_df) <- "SZ"
stim_df$Pathway <- as.vector(rownames(stim_df))

combine_df <- full_join(ctrl_df, stim_df, by="Pathway")
pathway_df <- full_join(df_pvalue, combine_df, by="Pathway")
rownames(pathway_df) <- as.vector(pathway_df$Pathway)
pathway_df$Pathway <- NULL

# metabolism pathway
significance_pathway <- subset(pathway_df, p.adj<0.05)
# write.table(significance_pathway, "SZ_output_metabolic pathway score.txt", sep = "\t", quote = F)


###可视化
# pathway_plot_df <- significance_pathway
pathway_plot_df <- read.table("~/精神疾病/data/388/Seurat/data/metabolism/PTSD_output_metabolic pathway score.txt", header = TRUE, sep = "\t")
pathway_plot_df$p_val <- NULL
pathway_plot_df$p.adj <- NULL
pathway_plot_df[,c(1:2)] <- as.numeric(unlist(pathway_plot_df[,c(1:2)]))
# pathway_plot_df <- pathway_plot_df[c(1:3,7:9),]

pdf("PTSD_significance_metabolism pathway_pheatmap.pdf",height = 20,width = 8)
pheatmap(pathway_plot_df,
         cluster_row = F,
         treeheight_row = "20",
         cluster_col = F,
         show_rownames=T,
         cellwidth = 40,
         cellheight = 20,
         border_color = NA,
         color = colorRampPalette(colors = c("#333189","#306faf","#efe92b"))(100),
         legend = T,
         scale="column") 
dev.off()


