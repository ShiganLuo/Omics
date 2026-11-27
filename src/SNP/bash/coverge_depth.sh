#!/bin/bash
bam=$1
outfile=$2
genomeSize=$3
function depth(){
    bam=$1
    sample_id=$(basename -s .bam ${bam})
    outfile=$2
    genomeSize=$3
    if [ ! -f ${outfile} ]; then
        echo -e "sample\tgenome_depth\tcoverage_depth\tcoverage\tcoverage_ratio" > ${outfile}
    fi
    samtools depth ${bam} | \
        awk -v sample=${sample_id} -v genomeSize=${genomeSize} '{sum+=$3} END { print sample  "\t" sum/genomeSize "\t" sum/NR "\t" NR "\t" NR/genomeSize } ' \
        >> ${outfile}
}
export -f depth

