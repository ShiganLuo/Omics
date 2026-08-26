library(dplyr)
library(Seurat)
library(patchwork)
options(bitmapType = "cairo")
#标准格式：行是基因，列是细胞;此数据相反，且以,分隔，数据需要进行转换
scRNA_py = function(path,sample){
    raw_data = read.csv(path,header = TRUE,row.names = 1, stringsAsFactors = FALSE, sep = ",")
    # 转置数据框并保留原始行名和列名
    transposed_data <- t(raw_data)
    transposed_data <- as.data.frame(transposed_data)
    # 将列名和行名重新命名为基因和细胞名称
    colnames(transposed_data) <- rownames(raw_data)  # 原始数据的行名是细胞名
    rownames(transposed_data) <- colnames(raw_data)  # 原始数据的列名是基因名

    df_Seurat = CreateSeuratObject(counts = transposed_data, project = sample, min.cells = 3, min.features = 200)
    # df_Seurat = CreateSeuratObject(counts = transposed_data, project = sample)
    return (df_Seurat)
}
scRNA_10x = function(path,sample){
    Sat.data = Read10X(path)
    Sat = CreateSeuratObject(counts = Sat.data, project = sample, min.cells = 3, min.features = 200)
    # Sat = CreateSeuratObject(counts = Sat.data, project = sample)
    # Sat = RenameCells(Sat, new.names = paste0(sample,"-",colnames(Sat),sep=""))
    Sat[['sample']] = sample
    return (Sat)
}
Transpon = function(path,sample){
    #scTE output always need to be tranformed
    raw_data = read.csv(path,header = TRUE,row.names = 1, stringsAsFactors = FALSE, sep = ",")
    transposed_data <- t(raw_data)
    transposed_data <- as.data.frame(transposed_data)
    # 将列名和行名重新命名为基因和细胞名称
    colnames(transposed_data) <-  paste(sample, "-",rownames(raw_data),sep="")  # 原始数据的行名是细胞名
    rownames(transposed_data) <- colnames(raw_data)  # 原始数据的列名是基因名
    df_Seurat = CreateSeuratObject(counts = transposed_data, project = sample, min.cells = 3, min.features = 200)
    # df_Seurat = CreateSeuratObject(counts = transposed_data, project = sample)
    df_Seurat[['sample']] = sample
    return (df_Seurat)
}
samples = c("GBM27","GBM28","GBM29")

SC_sample = list()
TE_sample = list()

for(sample in samples){
    SC_file = paste("/home/lsg/Data/glioblastoma/output/cellranger/",sample,"/outs/filtered_feature_bc_matrix",sep = "")
    # read data and Create Seurat object
    SC_sample[[sample]] = scRNA_10x(SC_file,sample)
    cat("\nRead file:",SC_file)
    cat("\n",sample,"scRNA-seq Seurat对象维度:",dim(SC_sample[[sample]]),"\n")
    TE_file = paste("/home/lsg/Data/glioblastoma/output/scTE/cellranger_yes/",sample,".csv",sep = "")
    # read data and Create Seurat object
    TE_sample[[sample]] = Transpon(TE_file,sample)
    cat("\nRead file",TE_file)
    cat("\n",sample,"TE Seurat对象维度:",dim(TE_sample[[sample]]))
}

SC = merge(x = SC_sample[[1]],y = SC_sample[-1],add.cell.ids = names(SC_sample))
TE = merge(x = TE_sample[[1]],y = TE_sample[-1],add.cell.ids = names(TE_sample))
cat("\nSC Seurat对象维度",dim(SC))
cat("\nTE Seurat对象维度",dim(TE),"\n")

rm(SC_sample,TE_sample,sample)


saveRDS(SC, file = "/home/lsg/Data/glioblastoma/output/rds/merge/SC.rds")
saveRDS(TE, file = "/home/lsg/Data/glioblastoma/output/rds/merge/TE.rds")
print("SC rds已保存")
print("TE rds已保存")
GBM = merge(SC,TE)
saveRDS(GBM, file = "/home/lsg/Data/glioblastoma/output/rds/merge/GBM.rds") 
cat("\nGBM Seurat对象维度",dim(GBM),"\n")
print("GBM rds已保存")
