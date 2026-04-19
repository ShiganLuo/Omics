#!/bin/bash
outdir="/data/pub/zhousha/20260411_RNAseq/output/test"
ngs_disambiguate \
    -s B1_CoD1_Mouse_WT_rep1 \
    -o ${outdir} \
    -a hisat2 \
    /data/pub/zhousha/20260411_RNAseq/output/hisat2/GRCh38/B1_CoD1_Mouse_WT_rep1.bam \
    /data/pub/zhousha/20260411_RNAseq/output/hisat2/GRCm39/B1_CoD1_Mouse_WT_rep1.bam