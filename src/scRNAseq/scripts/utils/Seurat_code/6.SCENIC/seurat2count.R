#### Information ----
# Title   :   Seurat object is transformed into expression matrix
# File    :   seurat2count.R
# Author  :   Songqi Duan
# Contact :   songqi.duan@outlook.com
# License :   Copyright (C) 2017-2021 by Songqi Duan | 段松岐
# Created :   2021/03/07 15:40:10
# Updated :   2021/05/25 23:07:17

pacman::p_load(Seurat)
pacman::p_load(optparse)

option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = FALSE,
        action = "store", help = "Seurat Object"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = FALSE,
        action = "store", help = "Output Path"
    ),
    make_option(c("-g", "--group"),
        type = "character", default = FALSE,
        action = "store", help = "Group"
    )
)
opt <- parse_args(OptionParser(
    option_list = option_list,
    usage = "This Script is a test for arguments!"
))

seurat_obj <- readRDS(opt$input)
DefaultAssay(seurat_obj) <- "RNA"

if(opt$group != "all") {
seurat_obj <- seurat_obj[, seurat_obj[["group"]]==opt$group]
}
counts <- as.matrix(seurat_obj@assays[["RNA"]]@counts)

if (!dir.exists(paste0(opt$output))) {
    dir.create(paste0(opt$output), recursive = T)
}

write.table(cbind(rownames(counts), counts),
    file = paste0(opt$output, "/counts.txt"),
    sep = "\t", quote = F,
    row.names = F
)
