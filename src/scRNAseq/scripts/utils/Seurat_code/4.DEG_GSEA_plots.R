
#导入包 ----
library(Seurat)
library(ggplot2)
library(tidyverse)
source('source("/home/xingwl/share/20220802_scrna-m6A/custom_plot_function.R")')
setwd('~/m6A/RNA_LAU/DEGs')

# file.forder <- c('L1_P2_vs_L1_P0', 'L1_P3_vs_L1_P0', 'L15_P2_vs_L15_P0', 'L15_P3_vs_L15_P0', 'P2_vs_P0', 'P3_vs_P0')

#### 1.P2_vs_P0 ----
k <- 'P2_vs_P0'
outdir <- paste0("./results/files/padj0.05_logfc1/FUN_GSEA/", k)
outdir2 <- paste0("./results/plots/padj0.05_logfc1/FUN_GSEA/", k)
#### C5 ----
res <- readRDS(paste0("./results/files/padj0.05_logfc1/FUN_GSEA/",k,"_C5_result.rds"))
res@result$Description <- res@result$Description|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("GOCC_", "") |>
  str_replace_all("GOMF_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
View(res@result)
selected_features <- c('GOBP_PANCREAS_DEVELOPMENT',# DOWN
                       'GOBP_MEIOTIC_CELL_CYCLE_PHASE_TRANSITION',
                       'GOBP_ADENOHYPOPHYSIS_DEVELOPMENT',
                       'GOBP_PITUITARY_GLAND_DEVELOPMENT',
                       'GOBP_REGULATION_OF_T_HELPER_1_TYPE_IMMUNE_RESPONSE',#UP
                       'GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY',
                       'GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_CYTOKINE_PRODUCTION',
                       'GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_4_PRODUCTION',
                       'GOBP_IMMUNOLOGICAL_SYNAPSE_FORMATION',
                       'GOMF_IMMUNE_RECEPTOR_ACTIVITY')
# i = 'GOBP_CANONICAL_WNT_SIGNALING_PATHWAY'
for(i in selected_features){
  p <- res |>
    cat_gseaplot(
      i,
      subplots = c(1, 2),
      pvalue_table = T,
      title = i|>
        str_replace_all("GOBP_", "") |>
        str_replace_all("_", " ") |>
        str_to_sentence())
  p
  ggsave(file.path(outdir2, paste0('C5_',i,'.pdf')),
         plot = p,
         height = 2.5,
         width = 3.5)
}

#### C2 ----
res <- readRDS(paste0("./results/files/padj0.05_logfc1/FUN_GSEA/",k,"_C2_result.rds"))
res@result$Description <- res@result$Description|>
  str_replace_all("KEGG_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
View(res@result)


#### C7 ----
res <- readRDS(paste0("./results/files/padj0.05_logfc1/FUN_GSEA/",k,"_C7_result.rds"))
res@result$Description <- res@result$Description|>
  str_replace_all("KEGG_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
View(res@result)


#### barplot ----
library(ggpubr)
res <- readRDS(paste0("./results/files/padj0.05_logfc1/FUN_GSEA/",k,"_C5_result.rds"))
data <- res@result[res@result$ID %in% selected_features,]
data <- data[order(data$NES,decreasing = T),]
# num <- 6
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
data$pvalue <- -log10(data$pvalue)
data$change <- c(rep('UP-regulated',6),rep('Down-regulated',4))
data$group <- c(rep(1,6),rep(-1,4))

p2 <- ggbarplot(data, x="Description", y="pvalue", fill = "change", color = "white",
          orientation = "horiz",  #横向显示
          palette = c("#00a6e1", "#ee6470"),    #配色方案
          legend = "top",    #图例位置
          sort.val = 'asc',#下降排序"dasc",    #上升排序"asc",
          sort.by.groups=TRUE)+#按照组排序
  geom_text(aes(x=labely,
                y=labelx,
                label = Description),
            size=4,
            hjust =0)+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 0.5),
        axis.title.x = element_text(colour = 'black', size = 12))+
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(0,0))+
  ylab("-log10(pvalue)")+xlab("")+
  labs(fill='')+
  ggtitle("")

ggsave(file.path(outdir2, 'C5_GSEA_barplot.pdf'), plot = p2, height = 6, width = 6)



#### 2.P3_vs_P0 ----
k <- 'P3_vs_P0'
# outdir <- paste0("./results/files/padj0.05_logfc1/FUN_GSEA/", k)
outdir2 <- paste0("./results/plots/padj0.05_logfc1/FUN_GSEA/", k)

res <- readRDS(paste0("./results/files/padj0.05_logfc1/FUN_GSEA/",k,"_C5_result.rds"))
res@result$Description <- res@result$Description|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("GOCC_", "") |>
  str_replace_all("GOMF_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
View(res@result)
selected_features <- c('GOBP_MEIOTIC_CELL_CYCLE_PHASE_TRANSITION', #DOWN
                       'GOCC_MEIOTIC_SPINDLE',
                       'GOBP_REGULATION_OF_T_HELPER_1_TYPE_IMMUNE_RESPONSE', #UP
                       'GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE',
                       'GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS',
                       'GOBP_POSITIVE_REGULATION_OF_NIK_NF_KAPPAB_SIGNALING',
                       # 'GOBP_IN_UTERO_EMBRYONIC_DEVELOPMENT',
                       'GOBP_REGULATION_OF_PRODUCTION_OF_SMALL_RNA_INVOLVED_IN_GENE_SILENCING_BY_RNA')
# i = 'GOBP_CANONICAL_WNT_SIGNALING_PATHWAY'
for(i in selected_features){
  p <- res |>
    cat_gseaplot(
      i,
      subplots = c(1, 2),
      pvalue_table = T,
      title = i|>
        str_replace_all("GOBP_", "") |>
        str_replace_all("_", " ") |>
        str_to_sentence())
  p
  ggsave(file.path(outdir2, paste0('C5_',i,'.pdf')),
         plot = p,
         height = 2.5,
         width = 3.5)
}

#### C2 ----
res <- readRDS(paste0("./results/files/padj0.05_logfc1/FUN_GSEA/",k,"_C2_result.rds"))
res@result$Description <- res@result$Description|>
  str_replace_all("KEGG_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
View(res@result)


#### C7 ----
res <- readRDS(paste0("./results/files/padj0.05_logfc1/FUN_GSEA/",k,"_C7_result.rds"))
res@result$Description <- res@result$Description|>
  str_replace_all("KEGG_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
View(res@result)


#### 2.4 barplot ----
library(ggpubr)
res <- readRDS(paste0("./results/files/padj0.05_logfc1/FUN_GSEA/",k,"_C5_result.rds"))
data <- res@result[res@result$ID %in% selected_features,]
data <- data[order(data$NES,decreasing = T),]
# num <- 6
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
data$pvalue <- -log10(data$pvalue)
data$change <- c(rep('UP-regulated',5),rep('Down-regulated',2))
data$group <- c(rep(1,5),rep(-1,2))

p2 <- ggbarplot(data, x="Description", y="pvalue", fill = "change", color = "white",
                orientation = "horiz",  #横向显示
                palette = c("#00a6e1", "#ee6470"),    #配色方案
                legend = "top",    #图例位置
                sort.val = 'asc',#下降排序"dasc",    #上升排序"asc",
                sort.by.groups=TRUE)+#按照组排序
  geom_text(aes(x=labely,
                y=labelx,
                label = Description),
            size=3,
            hjust =0)+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 0.5),
        axis.title.x = element_text(colour = 'black', size = 12))+
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(0,0))+
  ylab("-log10(pvalue)")+xlab("")+
  labs(fill='')+
  ggtitle("")

ggsave(file.path(outdir2, 'C5_GSEA_barplot.pdf'), plot = p2, height = 6, width = 6.3)

