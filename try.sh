#!/bin/bash
function trimSingle(){
    fastq=$1
    outdir=$2
    sample=$(basename -s .fq.gz ${fastq})
    trim_galore=/opt/TrimGalore-0.6.10/trim_galore
    ${trim_galore} \
        --cores 6 \
        --quality 30 \
        -o ${outdir} \
        --basename ${sample} \
        ${fastq}
}
export -f trimSingle
fastq=/ChIP_seq_2/StemCells/RNASNP_PT/data/GSE166216/fq/SRR13633379.fq.gz
outdir=./
trimSingle ${fastq} ${outdir}