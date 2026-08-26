library(tidyverse)
library(patchwork)
library(Seurat)
SC = readRDS("/home/lsg/Data/glioblastoma/output/rds/SC_QC.rds")
TE = readRDS("/home/lsg/Data/glioblastoma/output/rds/TE_QC.rds")
GBM = readRDS("/home/lsg/Data/glioblastoma/output/rds/GBM_QC.rds")
### the data need to be splited,whci include many samples
##v2
normalization1 = function(Sat){
    Sat_list = SplitObject(Sat,split.by = "sample")
    
    for(i in 1:length(Sat_list)){
        # the first normalization
        
        Sati = Sat_list[[i]]
        Sati = NormalizeData(Sati,normalization.method = "LogNormalize")
        # find high variable features
        Sati = FindVariableFeatures(Sati,selection.method = "vst",nfeatures = 3000)
        Sat_list[[i]] = Sati
        print(Sat_list[[i]])

    }
    # data Integration
    anchors = FindIntegrationAnchors(
        object.list = Sat_list,
        normalization.method = "LogNormalize",
        verbose = TRUE
    )
    Sat = IntegrateData(anchorset = anchors)
    DefaultAssay(Sat) = "integrated"
    # the fourth normalization
    Sat = ScaleData(Sat,features = rownames(Sat))
    return (Sat)

}
##version=v5,suggestd
normalization2 = function(Sat){
    Sat_list = SplitObject(Sat,split.by = "sample")
    for(i in 1:length(Sat_list)){
        Sat_list[[i]] = SCTransform(Sat_list[[i]],
            variable.features.n = 3000,
            verbose = FALSE)
    }
    features = SelectIntegrationFeatures(object.list = Sat_list,nfeatures = 3000)
    Sat_list = PrepSCTIntegration(object.list = Sat_list,anchor.features = features)
    anchors = FindIntegrationAnchors(object.list = Sat_list,
        normalization.method = "SCT",
        anchor.features = features
    )
    Sat = IntegrateData(anchorset = anchors,
        normalization.method = "SCT"
    )
    DefaultAssay(Sat) = "integrated"
    return (Sat)
}
SC = normalization1(SC)
print(str(SC))
# SC = normalization2(SC)
# print(str(SC))
# print(DefaultAssay(SC))
