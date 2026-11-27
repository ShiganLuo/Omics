#!/bin/bash
#建立注释库
retrieve_pl=/opt/annovar/retrieve_seq_from_fasta.pl
convert_pl=/opt/annovar/convert2annovar.pl
annotate_pl=/opt/annovar/annotate_variation.pl
table_pl=/opt/annovar/table_annovar.pl
set -e
function GenePred(){
    local gtf=$1
    local refGene=$2
    gtfToGenePred -genePredExt ${gtf} ${refGene}
}
export -f GenePred
# gtf=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/gencode.vM36.primary_assembly.annotation.gtf
# refGene=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/GRCm39_refGene.txt
# GenePred ${gtf} ${refGene}
# gtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/gencode.v47.primary_assembly.annotation.gtf
# refGene=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/annovar/GRCh38/GRCh38_refGene.txt
# GenePred ${gtf} ${refGene}
function retrieve(){
    local retrieve_pl=$1
    local genome=$2
    local refGene=$3
    local refGeneMrna=$4
    dir=$(dirname ${refGene})
    mkdir -p ${dir}
    perl ${retrieve_pl} \
    --format refGene \
    --seqfile ${genome} ${refGene} \
    --out ${refGeneMrna}
}
export -f retrieve

function convert(){
    local convert_pl=$1
    local vcf=$2
    local out=$3
    perl ${convert_pl} \
    -format vcf4 \
    -withfreq ${vcf} > ${out}
}
export -f convert
vcf=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/vcf/SRR17762738.bed
out=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/vcf/SRR17762738.avinput
# convert ${convert_pl} ${vcf} ${out}
function table(){
    local table_pl=$1
    local avinput=$2
    local db=$3
    local out=$4
    local buildver=$5
    perl ${table_pl} \
    ${avinput} ${db} \
    -buildver ${buildver} \
    -out ${out} \
    -remove -protocol refGene \
    -operation g \
    -nastring . \
    -csvout
}
export -f table
avinput=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/annovar/SRR17762738.avinput
db=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/
buildver=GRCm39
out=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/annovar/SRR17762738_table
