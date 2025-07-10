#!/bin/bash
#建立注释库
retrieve=/opt/annovar/retrieve_seq_from_fasta.pl
convert=/opt/annovar/convert2annovar.pl
annotate=/opt/annovar/annotate_variation.pl
table=/opt/annovar/table_annovar.pl
#conda:RNA-SNP
## notice: Again do not use annotate_variation unless you are an expert and know the intricate differences of the many arguments.
set -e
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
function GenePred(){
    gtf=$1
    refGene=$2
    gtfToGenePred -genePredExt ${gtf} ${refGene}
}
# gtf=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/gencode.vM36.primary_assembly.annotation.gtf
# refGene=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/GRCm39_refGene.txt
# GenePred ${gtf} ${refGene}
# gtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/gencode.v47.primary_assembly.annotation.gtf
# refGene=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/annovar/GRCh38/GRCh38_refGene.txt
# GenePred ${gtf} ${refGene}
function retrieve(){
    genome=$1
    refGene=$2
    refGeneMrna=$3
    perl ${retrieve} \
    --format refGene \
    --seqfile ${genome} ${refGene} \
    --out ${refGeneMrna}
}
# genome=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa
# refGene=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/GRCm39_refGene.txt
# refGeneMrna=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/GRCm39_refGeneMrna.fa
# retrieve ${genome} ${refGene} ${refGeneMrna}
# genome=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa
# refGene=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/annovar/GRCh38/GRCh38_refGene.txt
# refGeneMrna=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/annovar/GRCh38/GRCh38_refGeneMrna.fa
# retrieve ${genome} ${refGene} ${refGeneMrna}
function convert(){
    vcf=$1
    out=$2
    perl ${convert} \
    -format vcf4 \
    -withfreq ${vcf} > ${out}
}
vcf=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/vcf/SRR17762738.bed
out=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/vcf/SRR17762738.avinput
# convert ${vcf} ${out}
function table(){
    avinput=$1
    db=$2
    out=$3
    buildver=$4
    perl ${table} \
    ${avinput} ${db} \
    -buildver ${buildver} \
    -out ${out} \
    -remove -protocol refGene \
    -operation g \
    -nastring . \
    -csvout
}
avinput=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/annovar/SRR17762738.avinput
db=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/
buildver=GRCm39
out=/ChIP_seq_2/StemCells/RNASNP202503/waitingflow/output/annovar/SRR17762738_table
# table ${avinput} ${db} ${out} ${buildver} 