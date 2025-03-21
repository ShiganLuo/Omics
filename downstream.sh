#!/bin/bash
set -e
Rscript workflow/scripts/DESeq2/DESeq2.r --mode TElocal \
    --matrix output/counts/mouseTEcount.cntTable \
    --group workflow/assets/group.csv \
    --pattern "wild type" "Eif2ak4-/- (GCN2 knock-out)" \
    --outdir output \
    --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    --figure pca heatmap volcano

Rscript workflow/scripts/DESeq2/DESeq2.r --mode TElocal \
    --matrix output/counts/mouseTElocal.cntTable \
    --group workflow/assets/group.csv \
    --pattern "wild type" "Eif2ak4-/- (GCN2 knock-out)" \
    --outdir output \
    --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    --figure heatmap volcano