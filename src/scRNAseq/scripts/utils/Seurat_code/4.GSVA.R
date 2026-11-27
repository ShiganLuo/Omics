#### 加载包 ----
library(Seurat)
library(GSVA)
library(GSEABase)
library(limma)
library(tidyverse)
library(ggplot2)


setwd('~/精神疾病/data/388/Seurat/data/GSVA')
outdir_data <- paste0('./score/')
dir.create(outdir_data,recursive = T)
outdir_plots <- paste0('./plots/')
dir.create(outdir_plots,recursive = T)

s.sets <- getGmt("/DATA/public/MSigDB/7.5.1/h.all.v7.5.1.symbols.gmt")

sub.seurat.obj <- readRDS('~/精神疾病/data/388/Seurat/data/merged/merged_all_scted.rds')
sub.seurat.obj2 <- readRDS('~/精神疾病/data/388/Seurat/data/SCTed_data.rds')



expr <- as.matrix(sub.seurat.obj@assays$RNA@data)
es.matrix <- gsva(
  expr,
  s.sets,
  min.sz = 10,# 基因集中最小基因数，默认10，表示分析中仅考虑包含大于10个基因的基因集
  max.sz = Inf,# 基因集中最大基因数，默认inf，无限
  tau = 1,
  method = "gsva",
  abs.ranking = FALSE,#是否按基因表达绝对值排序
  verbose = TRUE,#是否输出结果
  parallel.sz = 1
)
es.matrix <- readRDS(paste0(outdir_data,"GSVA_hallmarker_matrix.rds"))


table(sub.seurat.obj$group)



case <- "SZ"
control <- "CON"
meta <-
  sub.seurat.obj@meta.data[, c("orig.ident",
                               "group")]
meta$cell_id <- rownames(meta)

cell_id.1 <-
  meta[meta$group == case, ]$cell_id
cell_id.2 <-
  meta[meta$group == control, ]$cell_id

es.matrix.1 <-
  as.data.frame(es.matrix[, cell_id.1],
                row.names = row.names(es.matrix))
es.matrix.2 <-
  as.data.frame(es.matrix[, cell_id.2],
                row.names = row.names(es.matrix))
es.matrix.f <- cbind(es.matrix.1, es.matrix.2)

class(es.matrix)
class(es.matrix.f)


grouP <-
  c(rep("case", dim(es.matrix.1)[2]), rep("control", dim(es.matrix.2)[2])) %>% as.factor()

desigN <- model.matrix( ~ grouP + 0)
rownames(desigN) <-
  c(colnames(es.matrix.1), colnames(es.matrix.2))
comparE <-
  makeContrasts(grouPcase - grouPcontrol, levels = desigN)
fiT <- lmFit(es.matrix.f, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
diff <- topTable(fiT3, coef = 1, number = 200)####做gobp打分时可以改一下参数，改大点，目前太小

t_results <-
  as.data.frame(diff$t, row.names = rownames(es.matrix.f))
head(t_results)
colnames(t_results) <- c("t_value") # 调整列名
file_path <- file.path(path.expand("~/精神疾病/data/388/Seurat/data/GSVA/score"), "SZ_GSVA_hallmarker_t_value.rds")  
saveRDS(t_results, file_path)



Cell <- c('SZ') 

rownames(t_results) <- gsub("HALLMARK_", "", rownames(t_results))
rownames(t_results) <- gsub("_", " ", rownames(t_results))
rownames(t_results) <- tolower(rownames(t_results))
t_results$hallmark <- rownames(t_results)
colnames(t_results) <- c("t", "hallmark")
head(t_results)

t_results$hallmark = with(t_results, reorder(hallmark, t))
t_results$fill <- ""
t_results[t_results$t >= 2.58,]$fill <-
  "up"
t_results[t_results$t <= -2.58,]$fill <-
  "down"
t_results[abs(t_results$t) < 2.58,]$fill <-
  "no"
t_results$color <- ""
t_results[abs(t_results$t) < 2.58,]$color <-
  "n"
t_results[abs(t_results$t) >= 2.58,]$color <-
  "y"

p <-
  ggplot(t_results, aes(x = hallmark, y = t, fill = fill)) +
  geom_bar(stat = "identity", width = .8) +
  scale_fill_manual(
    values = c(
      "down" = "#6fa7d3",
      "up" = "#da535a",
      "no" = "#cccccc"
    ),
    guide = F
  ) +
  geom_hline(
    yintercept = c(-2.58, 2.58),
    color = "white",
    linetype = "dashed",
    size = 0.5
  ) +
  coord_flip() +
  xlab("") +
  ylab("t value of GSVA score") +
  geom_text(
    data = subset(t_results, t < 0),
    aes(
      x = hallmark,
      y = 0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 3,
    hjust = "inward" # outward 在条内  inward 在条外
  ) +
  geom_text(
    data = subset(t_results, t > 0),
    aes(
      x = hallmark,
      y = -0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 3,
    hjust = "outward" # inward 在条内  outward 在条外
  ) +
  scale_colour_manual(values = c("y" = "black", "n" = "#cccccc"), # y-
                      guide = FALSE) +
  labs(title=paste0(Cell,' vs. CON'))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust=0.5,size = 12),
    axis.text.y = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),
  )
p

ggsave(
  filename = paste0(outdir_plots,"SZ_GSVA_bar.pdf"),
  plot = p,
  width = 12,
  height = 6
)
