#!/bin/bash
set -e
m=$1
n=$2
sradir=/home/lsg/Data/glioblastoma/data/sra
fqdir=/home/lsg/Data/glioblastoma/data/fq
#数组索引是从零开始的
#start代表开始索引0~27
#n代表从start开始取几个元素
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
    #prefetch
    time prefetch --output-directory ${sradir} ${sraId}
    mv ${sradir}/${sraId}.sra ${sradir}/${renameId}.sra
    #parallel-fastq-dump
    echo "${i}开始提取fsataq文件"
    # time fastq-dump --outdir ${fqdir} --split-files --gzip ${fqdir}/${renameId}.sra
    time parallel-fastq-dump --sra-id ${sradir}/${renameId}.sra --threads 24 --outdir ${fqdir} --split-files --gzip
    echo "${i}的fastaq文件提取完毕"
    mv ${fqdir}/${renameId}_1.fastq.gz ${fqdir}/${renameId}_S1_L001_R1_001.fastq.gz 
    mv ${fqdir}/${renameId}_2.fastq.gz ${fqdir}/${renameId}_S1_L001_R2_001.fastq.gz 
    
done
echo "执行结束"


