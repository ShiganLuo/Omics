#!/bin/bash
function genome_download(){
    ID=$1
    outdir=$2
    datasets download genome accession ${ID} --include gff3,genome,gtf
    mv ncbi_dataset.zip ${outdir}
    unzip ${outdir}/ncbi_dataset.zip
    md5sum -c ${outdir}/md5sum.txt
    echo "文件现在并校验成功"
}
genome_dir=/ChIP_seq_2/StemCells/data/genome
##1)
# genome_download GCF_000001405.40 ${genome_dir}
##2)
# wget -P ${genome_dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz
# wget -P ${genome_dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gff3.gz
# wget -P ${genome_dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
# wget -P ${genome_dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.annotation.gtf.gz
# wget -P ${genome_dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.annotation.gff3.gz
# wget -P ${genome_dir} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz

function star_index(){
    genome_index=$1
    genome_fa=$2
    genome_gtf=$3
    read_lenth=$4
    index_core=$5
    ${STAR} --runMode genomeGenerate \
    --runThreadN ${index_core} \
    --genomeDir ${genome_index} \
    --genomeFastaFiles ${genome_fa} \
    --sjdbGTFfile ${genome_gtf} 
}
#    --sjdbOverhang ${read_lenth} 理论上应该是max(read_lenth)-1,实际影响效果微乎其微
export -f star_index
STAR=/opt/STAR-2.7.11b/bin/Linux_x86_64/STAR
genome_index=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/star
genome_fa=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa
genome_gtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/gencode.v47.primary_assembly.annotation.gtf
read_lenth=85
index_core=25
# star_index ${genome_index} ${genome_fa} ${genome_gtf} ${read_lenth} ${index_core}
genome_index=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/star
genome_fa=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa
genome_gtf=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/gencode.vM36.primary_assembly.annotation.gtf
# star_index ${genome_index} ${genome_fa} ${genome_gtf} ${read_lenth} ${index_core}
# 对照：GSM5837959	Wt_Unstarved1  GSM5837962	Wt_Unstarved2 GSM5837965	Wt_Unstarved3
# 样品：GSM5837968	Gcn_Unstarved1 GSM5837971	Gcn_Unstarved2 GSM5837974	Gcn_Unstarved3
source /ChIP_seq_2/StemCells/RNASNP202503/workflow/scripts/ascp.sh
sralist=/ChIP_seq_2/StemCells/RNASNP202503/data/sra.txt
log=/ChIP_seq_2/StemCells/RNASNP202503/output/log/GEOdownload.log
outdir=/ChIP_seq_2/StemCells/RNASNP202503/data/fq
# : > ${log}
# awk -F"\t" '{print $3,$4}' ${sralist} | \
#     parallel -j 6 --colsep ' ' \
#     GEOdownload {1} {2} ${log} ${outdir}
# sralist=/ChIP_seq_2/StemCells/RNASNP202503/output/log/GEOdownload.log
# log=/ChIP_seq_2/StemCells/RNASNP202503/output/log/GEOdownload1.log
# : > ${log}
# awk '$3 == "下载失败" {print $1,$2}' ${sralist} | \
#     parallel -j 2 --colsep ' ' \
#     GEOdownload {1} {2} ${log} ${outdir}
# log=/ChIP_seq_2/StemCells/RNASNP202503/output/log/GEOdownload2.log
# : > ${log}
# echo -e "SRR17762739 PAIRED" | \
#     parallel -j 2 --colsep ' ' \
#     GEOdownload {1} {2} ${log} ${outdir}
function TEtranscripts(){
    bam=$1
    gtf=$2
    TEgtf=$3
    TEcount --sortByPos --format BAM --mode multi -b ${bam} --GTF ${gtf} --TE ${TEgtf} --project RNA
}
export -f TEtranscripts
bam=/ChIP_seq_2/StemCells/RNASNP202503/output/counts/SRR17762739/mouse/SRR17762739Aligned.sortedByCoord.out.bam
gtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/gene.gtf
TEgtf=/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf
TEtranscripts ${bam} ${gtf} ${TEgtf}