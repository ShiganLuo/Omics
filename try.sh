#!/bin/bash
function trimSingle(){
    fast1=$1
    outdir=$2
    trim_galore=/opt/TrimGalore-0.6.10/trim_galore
    ${trim_galore} --phred33 \
        --cores 6 \
        --quality 30 \
        -o ${outdir} \
        --basename ${sample} \
        ${fastq}
}
export -f trimSingle
fastq=