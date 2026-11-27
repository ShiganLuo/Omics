#### 导入包 ----
library(Seurat)

dir.create("~/精神疾病/data/388/Seurat/data/SCENIC/ASD_celltype")

txt5 <- read.table("~/精神疾病/data/388/Seurat/data/merged/SZ_merged_data_.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)



txt5 <- txt5[sample(colnames(txt5), 400)]
saveRDS(txt5, file = "~/精神疾病/data/388/Seurat/data/SCENIC/disease/SZ_test.rds")
