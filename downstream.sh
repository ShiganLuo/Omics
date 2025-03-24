#!/bin/bash
#RNA-SNP
# after RNASNP-noRemov.smk
set -e
Rscript workflow/scripts/DESeq2/DESeq2.r --mode TEcount \
    --matrix output/counts/mouseTEcount.cntTable \
    --group workflow/assets/group.csv \
    --pattern "wild type" "Eif2ak4-/- (GCN2 knock-out)" \
    --outdir output/ \
    --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    --figure pca heatmap volcano \
    --TEcountMode all

Rscript workflow/scripts/DESeq2/DESeq2.r --mode TElocal \
    --matrix output/counts/mouseTElocal.cntTable \
    --group workflow/assets/group.csv \
    --pattern "wild type" "Eif2ak4-/- (GCN2 knock-out)" \
    --outdir output/ \
    --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    --figure heatmap volcano
Rscript wworkflow/scripts/go.r
Rscript workflow/scripts/kegg.r
Rscript workflow/scripts/gsea.r \
    -i output/DESeq2/TEcount_Gene.csv \
    -g waitingflow/data/mh.all.v2024.1.Mm.symbols.gmt \
    -o output/ \
    -a /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    -t "Hallmark pathways altered by GCN2 in untreated mouse"

gseapy prerank \
    -r output/DESeq2/GSEA/TEcount_Gene_GSEA.rnk \
    -g data/gsea/TwoGenes.gmt \
    -o output/DESeq2/GSEA/TwoGenes \
    -f jpeg
