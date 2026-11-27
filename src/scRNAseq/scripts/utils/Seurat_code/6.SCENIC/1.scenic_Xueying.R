library(IRkernel)
library(Seurat)
library(tidyverse)
library(tictoc)
library(stringr)
library(data.table)






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



seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/ASD/ASD.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/ASD")

seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/BD/BD.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/BD")


seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/CON/CON.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/CON")


seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/MDD/MDD.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/MDD")


seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/PTSD/PTSD.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/PTSD")


seurat_obj <- readRDS("~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ/SZ.rds")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "~/精神疾病/data/388/Seurat/data/SCENIC/A_new/SZ")



