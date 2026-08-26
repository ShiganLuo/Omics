#### 导入包 ----
library(Seurat)

#### 配置 ----
dir.create("~/精神疾病/data/388/Seurat/data/SCENIC/ASD_celltype")

#### 导入数据 ----
txt5 <- read.table("~/精神疾病/data/388/Seurat/data/merged/SZ_merged_data_.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#### 提取某个细胞亚型 ----
# sub.seurat.obj <- subset(
#   seurat.obj1,
#   idents = "MDD"
# )

# sub.seurat.obj <- sub.seurat.obj[sample(rownames(sub.seurat.obj), 400), 
#                                  sample(colnames(sub.seurat.obj), 400)]

txt5 <- txt5[sample(colnames(txt5), 400)]
#### 保存细胞亚型 ----
saveRDS(txt5, file = "~/精神疾病/data/388/Seurat/data/SCENIC/disease/SZ_test.rds")
