#!/bin/bash
Rscript /data/pub/zhousha/20260411_MERIPseq/workflow/RNA-SNP/modules/exomePeak/bin/exomePeak.r \
    --gtf /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf \
    --ip_bams /data/pub/zhousha/20260411_MERIPseq/output/igv/dedup/SRR29583211.dedup.bam /data/pub/zhousha/20260411_MERIPseq/output/igv/dedup/SRR29583210.dedup.bam \
    --input_bams /data/pub/zhousha/20260411_MERIPseq/output/igv/dedup/SRR29583209.dedup.bam /data/pub/zhousha/20260411_MERIPseq/output/igv/dedup/SRR29583208.dedup.bam \
    --treated_ip_bams /data/pub/zhousha/20260411_MERIPseq/output/igv/dedup/SRR29583207.dedup.bam /data/pub/zhousha/20260411_MERIPseq/output/igv/dedup/SRR29583206.dedup.bam \
    --treated_input_bams /data/pub/zhousha/20260411_MERIPseq/output/igv/dedup/SRR29583205.dedup.bam /data/pub/zhousha/20260411_MERIPseq/output/igv/dedup/SRR29583204.dedup.bam \
    --outprefix /data/pub/zhousha/20260411_MERIPseq/output/exomePeak/Dox

