#!/bin/bash
#建立注释库
retrieve_pl=/opt/annovar/retrieve_seq_from_fasta.pl
convert_pl=/opt/annovar/convert2annovar.pl
annotate_pl=/opt/annovar/annotate_variation.pl
table_pl=/opt/annovar/table_annovar.pl
#conda:RNA-SNP
## notice: Again do not use annotate_variation unless you are an expert and know the intricate differences of the many arguments.
#############1
# echo "gffread BettaSplendens/genomic.gff -T -o BettaSplendens/BettaSplendens.gtf"
# gffread BettaSplendens/genomic.gff -T -o BettaSplendens/BettaSplendens.gtf
# #############2
# echo "gtfToGenePred -genePredExt BettaSplendens/BettaSplendens.gtf BettaSplendens/BettaSplendens_refGene.txt"
# gtfToGenePred -genePredExt BettaSplendens/BettaSplendens.gtf BettaSplendens/BettaSplendens_refGene.txt
# #############3
# echo "perl ${retrieve} \
#   --format refGene \
#   --seqfile BettaSplendens/GCF_900634795.4_fBetSpl5.4_genomic.fna BettaSplendens/BettaSplendens_refGene.txt \
#   --out BettaSplendens/BettaSplendens_refGeneMrna.fa"
# perl ${retrieve} \
#   --format refGene \
#   --seqfile BettaSplendens/GCF_900634795.4_fBetSpl5.4_genomic.fna BettaSplendens/BettaSplendens_refGene.txt \
#   --out BettaSplendens/BettaSplendens_refGeneMrna.fa
# ##############4
# echo "bgzip -d -c /home/lsg/GWAS/Betta_splendens240901/data/Splendens_615ind_maf0.02miss0.5minQ30allele2.onlyNumberChr.pass.vcf.gz > BettaSplendens/Splendens_615ind_maf0.02miss0.5minQ30allele2.onlyNumberChr.pass.vcf"
# bgzip -d -c /home/lsg/GWAS/Betta_splendens240901/data/Splendens_615ind_maf0.02miss0.5minQ30allele2.onlyNumberChr.pass.vcf.gz > BettaSplendens/Splendens_615ind_maf0.02miss0.5minQ30allele2.onlyNumberChr.pass.vcf
# #############5
# echo "perl ${convert} \
#     -format vcf4 -allsample -withfreq BettaSplendens/Splendens_615ind_maf0.02miss0.5minQ30allele2.onlyNumberChr.pass.vcf > BettaSplendens.avinput"
# perl ${convert} \
#     -format vcf4 -allsample -withfreq BettaSplendens/Splendens_615ind_maf0.02miss0.5minQ30allele2.onlyNumberChr.pass.vcf > BettaSplendens.avinput
# #############6
# echo "perl ${annotate} \
#     -geneanno -dbtype refGene \
#     -out Betta_Splendens.annotate \
#     -build BettaSplendens BettaSplendens.avinput BettaSplendens/"
# perl ${annotate} \
#     -geneanno -dbtype refGene \
#     -out Betta_Splendens.annotate \
#     -build BettaSplendens BettaSplendens.avinput BettaSplendens/ 
############8
# echo "perl ${table} \
#     BettaSplendens.avinput BettaSplendens/ \
#     -buildver BettaSplendens \
#     -out Splendens \
#     -remove -protocol refGene \
#     -operation g \
#     -nastring . \
#     -csvout"
# perl ${table} \
#     BettaSplendens.avinput BettaSplendens/ \
#     -buildver BettaSplendens \
#     -out Betta_Splendens.table \
#     -remove -protocol refGene \
#     -operation g \
#     -nastring . \
#     -csvout
set -e
function GenePred(){
    local gtf=$1
    local refGene=$2
    gtfToGenePred -genePredExt ${gtf} ${refGene}
    # mamba install bioconda::ucsc-gtftogenepred
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
# genome=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa
# refGene=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/GRCm39_refGene.txt
# refGeneMrna=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/GRCm39_refGeneMrna.fa
# retrieve ${retrieve_pl} ${genome} ${refGene} ${refGeneMrna}

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
# table ${table_pl} ${avinput} ${db} ${out} ${buildver} 