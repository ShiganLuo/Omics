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
# trimSingle ${fastq} ${outdir}
function xeno(){
    script=scripts/XenofilteR.r
    csvIn=../output/GSE204801/xenofilterR/xenofilterR_input.csv
    csvRe=../output/GSE204801/xenofilterR/xenofilterR_reName.csv
    outdir=../output/GSE204801/xenofilterR/bam
    /usr/bin/Rscript ${script} \
    --inputFile ${csvIn} \
    --outputDir ${outdir} \
    --MM 8 \
    --workers 3
}
# xeno
function stringTie(){
    gtf=$1
    bam=$2
    refGtf=$3
    stringtie -o ${gtf} ${bam} -G ${refGtf} -p 8
}
gtf=SRR13633379.gtf
bam=/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/2pass/SRR13633379/mouse/SRR13633379Aligned.sortedByCoord.out.bam
refGtf=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/exon.gtf
# stringTie ${gtf} ${bam} ${refGtf}
# gffcompare -r ${refGtf} -o comparison ${gtf}
