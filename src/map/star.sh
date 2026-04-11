#!/bin/bash
function star_index(){
    local genome_index=$1
    local genome_fa=$2
    local genome_gtf=$3
    local read_lenth=$4
    local index_core=$5
    local STAR=$6
    ${STAR} --runMode genomeGenerate \
    --runThreadN ${index_core} \
    --genomeDir ${genome_index} \
    --genomeFastaFiles ${genome_fa} \
    --sjdbGTFfile ${genome_gtf} \
    --sjdbOverhang ${read_lenth} # readLength-1
}
export -f star_index
#### example
# genome_index=/ChIP_seq_2/Data/F24A040011431_HOMnkgoT_0115/data/index/199
# genome_fa=/ChIP_seq_2/StemCells/data/genome/GRCh38.primary_assembly.genome.fa
# genome_gtf=/ChIP_seq_2/StemCells/data/genome/gencode.v47.primary_assembly.annotation.gtf
# read_lenth=199
# index_core=25
# star_index ${genome_index} ${genome_fa} ${genome_gtf} ${read_lenth} ${index_core}

### star_rRNAIndex function introduction
##don't need gtf,because what we care is map ,not transcripts annotation
function star_rRNAIndex(){
    STAR=/opt/STAR-2.7.11b/bin/Linux_x86_64/STAR
    rRNA_index=$1
    rRNA_fa=$2
    index_core=$3
    ${STAR} --runMode genomeGenerate \
     --genomeDir ${rRNA_index} \
     --genomeFastaFiles ${rRNA_fa} \
     --runThreadN ${index_core} \
     --genomeSAindexNbases 6
}
export -f star_rRNAIndex

function star_align(){
    STAR=/opt/STAR-2.7.11b/bin/Linux_x86_64/STAR
    genome_index=$1
    index_core=$2
    read1=$3
    read2=$4
    outPrefix=$5
    ${STAR} --runThreadN ${index_core} \
    --genomeDir ${genome_index} \
    --readFilesCommand zcat \
    --readFilesIn ${read1} ${read2} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${outPrefix} \
    --outReadsUnmapped Fastx \
    --quantMode GeneCounts
}
export -f star_align
#### example
# genome_index=/ChIP_seq_2/Data/F24A040011431_HOMnkgoT_0115/data/index/99
# index_core=25
# read1=/ChIP_seq_2/Data/F24A040011431_HOMnkgoT_0115/data/serumA/serumA_1_val_1.fq.gz
# read2=/ChIP_seq_2/Data/F24A040011431_HOMnkgoT_0115/data/serumA/serumA_2_val_2.fq.gz
# outPrefix=/ChIP_seq_2/Data/F24A040011431_HOMnkgoT_0115/output/star/serumA
# star_align ${genome_index} ${index_core} ${read1} ${read2} ${outPrefix}