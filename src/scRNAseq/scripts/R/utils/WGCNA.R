library(sceasy)
library(reticulate)
sceasy::convertFormat("/home/lsg/Data/glioblastoma/output/new/h5ad/SC-bbknn.h5ad", from="anndata", to="seurat",
                       outFile='filename.rds')