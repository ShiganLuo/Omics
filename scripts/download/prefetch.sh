#!/bin/bash
sralist=/ChIP_seq_2/StemCells/log/download/ascp_download2.log
function NCBIdownload(){
    ID=$1
    outdir=$2
    log=$3
    pre_donwload(){
        prefetch $1 -O ${outdir}
        mv ${1}/${1}.sra ./
        fastq-dump --split-3 --gzip ${1}.sra
        #--split-files

    }
    if pre_donwload $ID -O /ChIP_seq_2/StemCells/data/fq;then
        echo "${ID} 下载成功" >> ${log}
    else
        echo "${ID} 下载失败" >> ${log}
    

}
export -f NCBIdownload
# outdir=/ChIP_seq_2/StemCells/data/fq
# nohup awk -F " " '$3 == "下载失败" {print $1,$2}' ${sralist} | parallel -j 3 --joblog /ChIP_seq_2/StemCells/log/download/parallel.log --colsep ' ' download {1} ${outdir} ${log} > /ChIP_seq_2/StemCells/log/download/ascp.log 2>&1 &