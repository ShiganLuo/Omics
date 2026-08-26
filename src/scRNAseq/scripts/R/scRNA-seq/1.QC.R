library(tidyverse)
library(patchwork)
library(Seurat)
options(bitmapType = "cairo")


# add mitochondrial information to SeuratObject
QC_plot = function(Sat,out){
    Sat[["percent.mt"]] = PercentageFeatureSet(
        Sat,
        pattern = "^MT-"
    )
    # Sat = CellCycleScoring(
    #     Sat,
    #     s.features = cc.genes$c.genes,
    #     g2m.features = cc.genes$g2m.genes
    # )#Seurat自带cc.genes
    print(colnames(Sat@meta.data))#"orig.ident"   "nCount_RNA"   "nFeature_RNA"Seurat自带
    # print(str(Sat))

    p = VlnPlot(Sat,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        group.by = "sample",
        log = T,
        pt.size = 0.1    
    )

    ggsave(paste("./figures/R/",out,"_violin_plot.png",sep=""), plot = p, width = 10, height = 6, dpi = 300)
    p = RidgePlot(
        object = Sat,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        log = T,
        ncol = 1,
        group.by = "sample"
    )
    ggsave(paste("./figures/R/",out,"_ridge_plot.png",sep=""), plot = p, width = 10, height = 6, dpi = 300)
    p1 = FeatureScatter(Sat,
        feature1 = "nFeature_RNA",
        feature2 = "nCount_RNA",
        group.by = "sample"    
    )
    p2 = FeatureScatter(Sat,
        feature1 = "nCount_RNA",
        feature2 = "percent.mt",
        group.by = "sample"
    )

    ggsave(paste("./figures/R/",out,"_scatter_plot.png",sep=""), plot = p1+p2, width = 10, height = 6, dpi = 300)

}
QC_filter = function(Sat,nfeature_min,nfeature_max,percent_max){
    ### filter
    # you should caculate MT like QC_plot
    Sat[["percent.mt"]] = PercentageFeatureSet(
        Sat,
        pattern = "^MT-"
    )
    print(colnames(Sat@meta.data))
    Sat = subset(Sat,subset = nFeature_RNA > nfeature_min & nFeature_RNA < nfeature_max & percent.mt < percent_max)
    return (Sat)
}

execute = function(Sat,outid,nfeature_min,nfeature_max,percent_max){
    # plot for filter,so that you can specify the parments of QC_filter
    QC_plot(Sat,outid)
    #filter
    cat("\nbefore filter,the dimension of",outid,"is",dim(Sat),"\n")
    Sat = QC_filter(Sat,nfeature_min,nfeature_max,percent_max)
    cat("\nafter filter,the dimension of",outid,"is",dim(Sat))
    #save
    saveRDS(Sat, file = paste("/home/lsg/Data/glioblastoma/output/rds/QC/",outid,"_QC.rds",sep = ""))
    cat("\n",outid,"已保存\n")
}


SC = readRDS("/home/lsg/Data/glioblastoma/output/rds/merge/SC.rds")
TE = readRDS("/home/lsg/Data/glioblastoma/output/rds/merge/TE.rds")
GBM = readRDS("/home/lsg/Data/glioblastoma/output/rds/merge/GBM.rds") 

execute(SC,"SC",100,10000,10)
execute(TE,"TE",100,10000,3)
execute(GBM,"GBM",100,10000,10)
