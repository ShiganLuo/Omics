#!/bin/bash
# Rscript /disk5/luosg/Totipotent20251031/workflow/RNA-SNP/scripts/DESeq2/DESeq2_advance.r \
#     --mode TEcount \
#     --matrix /disk5/luosg/Totipotent20251031/output/counts/TEcount/mouse/all_TEcount.cntTable \
#     --group /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE185005.tsv \
#     --pattern mESC ciTotiSC \
#     --outdir output/result/ciTotiSC \
#     --annotation /disk5/luosg/Reference/GENCODE/mouse/GRCm39/geneIDAnnotation.csv \
#     --figure pca heatmap volcano \
#     --TEcountMode all

# Rscript /disk5/luosg/Totipotent20251031/workflow/RNA-SNP/scripts/DESeq2/DESeq2_advance.r \
#     --mode TEcount \
#     --matrix /disk5/luosg/Totipotent20251031/output/counts/TEcount/mouse/all_TEcount.cntTable \
#     --group /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE166216.tsv \
#     --pattern mESC TLSC \
#     --outdir output/result/TLSC \
#     --annotation /disk5/luosg/Reference/GENCODE/mouse/GRCm39/geneIDAnnotation.csv \
#     --figure pca heatmap volcano \
#     --TEcountMode all

# Rscript /disk5/luosg/Totipotent20251031/workflow/RNA-SNP/scripts/DESeq2/DESeq2_advance.r \
#     --mode TEcount \
#     --matrix /disk5/luosg/Totipotent20251031/output/counts/TEcount/human/all_TEcount.cntTable \
#     --group /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE224794.tsv \
#     --pattern hESC hTBLC \
#     --outdir output/result/hTBLC \
#     --annotation /disk5/luosg/Reference/GENCODE/human/GRCh38/geneIDAnnotation.csv \
#     --figure pca heatmap volcano \
#     --TEcountMode all

# Rscript /disk5/luosg/Totipotent20251031/workflow/RNA-SNP/scripts/DESeq2/DESeq2_advance.r \
#     --mode TEcount \
#     --matrix /disk5/luosg/Totipotent20251031/output/counts/TEcount/human/all_TEcount.cntTable \
#     --group /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE204801.tsv \
#     --pattern prEpiSC ci8CLC \
#     --outdir output/result/ci8CLC \
#     --annotation /disk5/luosg/Reference/GENCODE/human/GRCh38/geneIDAnnotation.csv \
#     --figure pca heatmap volcano \
#     --TEcountMode all
Rscript /disk5/luosg/Totipotent20251031/workflow/RNA-SNP/scripts/DESeq2/DESeq2_advance.r \
    --mode TEcount \
    --matrix /disk5/luosg/Totipotent20251031/PRJNA663159/TEcount/mouse/allTEcount.cntTable \
    --group /disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE204801.tsv \
    --pattern prEpiSC ci8CLC \
    --outdir output/result/ci8CLC \
    --annotation /disk5/luosg/Reference/GENCODE/human/GRCh38/geneIDAnnotation.csv \
    --figure pca heatmap volcano \
    --TEcountMode all