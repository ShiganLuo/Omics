#!/bin/bash
#RNA-SNP
# after RNASNP-noRemov.smk
# you should modify workflow/assets/group.csv to match your actual group
# also --pattern in the following command,the first is control, the second is experiment
set -e
function downstream_mouse(){
    pattern=$1
    outdir=$2
    group=$3
    title=$4 # gene kd/ko
    ###
    # Rscript workflow/scripts/DESeq2/DESeq2.r --mode TEcount \
    #     --matrix output/counts/mouseTEcount.cntTable \
    #     --group workflow/assets/group.csv \
    #     --pattern "siCtr" "siMat2a" \
    #     --outdir output/ \
    #     --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    #     --figure pca heatmap volcano \
    #     --TEcountMode all
    Rscript workflow/scripts/DESeq2/DESeq2.r --mode TEcount \
        --matrix output/${outdir}/counts/mouseTEcount.cntTable \
        --group ${group} \
        --pattern ${pattern} \
        --outdir output/${outdir}/ \
        --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
        --figure pca heatmap volcano \
        --TEcountMode all

    # Rscript workflow/scripts/DESeq2/DESeq2.r --mode TElocal \
    #     --matrix output/counts/mouseTElocal.cntTable \
    #     --group workflow/assets/group.csv \
    #     --pattern "siCtr" "siMat2a" \
    #     --outdir output/ \
    #     --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    #     --figure heatmap volcano
    Rscript workflow/scripts/DESeq2/DESeq2.r --mode TElocal \
        --matrix output/${outdir}/counts/mouseTElocal.cntTable \
        --group ${group} \
        --pattern ${pattern} \
        --outdir output/${outdir}/ \
        --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
        --figure heatmap volcano
    
    # Rscript workflow/scripts/DESeq2/go-kegg.r \
    #     --mode go kegg \
    #     --up output/DESeq2/upDown/TEcount_Gene_up.csv \
    #     --down output/DESeq2/upDown/TEcount_Gene_down.csv \
    #     --outdir output/ \
    #     --title "Mat2a knock-down"
    Rscript workflow/scripts/DESeq2/go-kegg.r \
        --mode go kegg \
        --up output/${outdir}/DESeq2/upDown/TEcount_Gene_up.csv \
        --down output/${outdir}/DESeq2/upDown/TEcount_Gene_down.csv \
        --outdir output/${outdir}/ \
        --title ${title}

    # Rscript workflow/scripts/DESeq2/gsea.r \
    #     -i output/DESeq2/TEcount_Gene.csv \
    #     -g data/gsea/mh.all.v2024.1.Mm.symbols.gmt \
    #     -o output/ \
    #     -a /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    #     -t "Hallmark pathways altered by Mat2a in mouse"
    mainTitle="Hallmark pathways altered after ${title} in mouse"
    Rscript workflow/scripts/DESeq2/gsea.r \
        -i output/${outdir}/DESeq2/TEcount_Gene.csv \
        -g data/gsea/mh.all.v2024.1.Mm.symbols.gmt \
        -o output/${outdir}/ \
        -a /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
        -t ${mainTitle}

    # gseapy prerank \
    #     -r output/DESeq2/GSEA/TEcount_Gene_GSEA.rnk \
    #     -g data/gsea/TwoGenes.gmt \
    #     -o output/DESeq2/GSEA/TwoGenes \
    #     -f jpeg
    gseapy prerank \
        -r output/${outdir}/DESeq2/GSEA/TEcount_Gene_GSEA.rnk \
        -g data/gsea/TwoGenes.gmt \
        -o output/${outdir}/DESeq2/GSEA/TwoGenes \
        -f jpeg

}
export -f downstream_mouse
pattern=""
downstream_mouse

