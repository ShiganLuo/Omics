library(IRkernel)
library(Seurat)
library(tidyverse)
library(tictoc)
library(stringr)
library(data.table)


# ##------------------抽样（每个疾病400个细胞）------------------##
# seurat <- readRDS("~/精神疾病/data/388/Seurat/data/SCTed_new.rds")
# scRNA <- seurat
# # 
# scRNA1 <- subset(x = scRNA, group == "CON")  
# subcell1 <- sample(colnames(scRNA1),400)
# scRNA2 <- subset(x = seurat, group == "ASD")  
# subcell2 <- sample(colnames(scRNA2),400)
# scRNA3 <- subset(x = scRNA, group == "BD")  
# subcell3 <- sample(colnames(scRNA3),400)
# scRNA4 <- subset(x = scRNA, group == "MDD")  
# subcell4 <- sample(colnames(scRNA4),400)
# scRNA5 <- subset(x = scRNA, group == "PTSD")  
# subcell5 <- sample(colnames(scRNA5),400)
# scRNA6 <- subset(x = scRNA, group == "SZ")  
# subcell6 <- sample(colnames(scRNA6),400)
# # 
# # seurat <- scRNA[,subcell6]
# seurat <- scRNA[,c(subcell1,subcell2,subcell3,subcell4,subcell5,subcell6)]
# saveRDS(seurat,"~/精神疾病/data/388/Seurat/data/SCENIC/A_new/ALL/ALL.rds")




####----------------------------scenic分析----------------------------####
--------------------------------------------------------------------------
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/data.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/objects.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/rip-seq.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/rna-seq.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/scrna-seq.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/utilities.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/visualization.R")

seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/ALL/ALL.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/ALL")



# seurat_obj <- subset(seurat,idents = c("ASD"))
seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/ASD/ASD.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/ASD")

# seurat_obj <- subset(seurat,idents = c("BD"))
seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/BD/BD.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/BD")


# seurat_obj <- subset(seurat,idents = c("CON"))
seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/CON/CON.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/CON")


# seurat_obj <- subset(seurat,idents = c("MDD"))
seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/MDD/MDD.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/MDD")


# seurat_obj <- subset(seurat,idents = c("PTSD"))
seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/PTSD/PTSD.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/PTSD")


# seurat_obj <- subset(seurat,idents = c("SZ"))
seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/SZ.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ")



