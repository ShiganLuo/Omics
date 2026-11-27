#!/bin/bash
set -e
m=$1
n=$2
sradir=/home/lsg/Data/glioblastoma/data/sra
fqdir=/home/lsg/Data/glioblastoma/data/fq
sra=('SRR10353962' 'SRR10353961' 'SRR10353960')
rename=('GBM29' 'GBM28' 'GBM27')
${array[@]:m:n}
for i in {0..2}
do
    sraId=${sra[${i}]}
    renameId=${rename[${i}]}
    echo "${sraId} : ${renameId}"
    echo ${i}
    echo "下载批次：$i"
    time prefetch --output-directory ${sradir} ${sraId}
    mv ${sradir}/${sraId}.sra ${sradir}/${renameId}.sra
    echo "${i}开始提取fsataq文件"
    time parallel-fastq-dump --sra-id ${sradir}/${renameId}.sra --threads 24 --outdir ${fqdir} --split-files --gzip
    echo "${i}的fastaq文件提取完毕"
    mv ${fqdir}/${renameId}_1.fastq.gz ${fqdir}/${renameId}_S1_L001_R1_001.fastq.gz 
    mv ${fqdir}/${renameId}_2.fastq.gz ${fqdir}/${renameId}_S1_L001_R2_001.fastq.gz 
    
done
echo "执行结束"


