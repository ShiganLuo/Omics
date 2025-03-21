#!/bin/bash
### download function parameter introduction
## id=$1 should be run number,such as SRR16119550
## type=$2 should be LibraryLayout,such as PAIRED
## log=$3 should log path
## dest=$4 should be outdir path
function GEOdownload(){
    #parllel子进程可能访问不到外部变量
    key=~/miniconda3/envs/RNA-SNP/etc/asperaweb_id_dsa.openssh
    id=$1
    type=$2
    log=$3
    dest=$4
    echo "${id} ${type}开始处理"
    x6=${id:0:6}
    x2=${id: -2}
    x2="0${x2}"
    echo ${x6}
    echo ${x2}
    file_path=(
    "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${x6}/${x2}/${id}/${id}_1.fastq.gz"
    "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${x6}/${x2}/${id}/${id}_2.fastq.gz"
    "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${x6}/${x2}/${id}/${id}.fastq.gz"
    )
    echo ${file_path[@]}
    function ascp_download(){
        ascp -k 1 -T -l 200m -P 33001 --file-checksum=md5 --overwrite=always -i ${key} $1 ${dest}
    }
    if [[ ${type} == "PAIRED" ]];then
        echo "双端测序"
        if ascp_download ${file_path[0]} && ascp_download ${file_path[1]};then
            echo -e "${id} ${type} 下载成功\n" >> ${log}
        elif ascp_download ${file_path[2]};then
            echo -e "${id} ${type}虽然为双端文件，但是可能只代表其中一个。更换下载地址后下载成功\n" >> ${log}
        else
            echo -e "${id} ${type} 下载失败\n" >> ${log}
        fi
    else
        echo "单端测序"
        if ascp_download ${file_path[2]};then
            echo -e "${id} ${type} 下载成功\n" >> ${log}
        else
            echo -e "${id} ${type} 下载失败\n" >> ${log}
        fi
    fi
}
# 导出函数
export -f GEOdownload
# sralist=/ChIP_seq_2/StemCells/log/download/check.log
# log=/ChIP_seq_2/StemCells/log/callSNP/ascp_download1.log
# : > ${log}
# nohup awk -F " " '$3 == "fq文件不存在" {print $1,$2}' ${sralist} | parallel -j 5 --joblog /ChIP_seq_2/StemCells/log/download/parallel.log --colsep ' ' GEOdownload {1} {2} ${log} > /ChIP_seq_2/StemCells/log/callSNP/ascp.log 2>&1 &