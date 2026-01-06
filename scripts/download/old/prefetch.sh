#!/bin/bash
### NCBIdownload function introduction
##ID=$1 the sra run number
##outdir=$2 output dir
##log=$3 log file
function NCBIdownload() {
    local ID="$1"
    local outdir="$2"
    local log="$3"

    mkdir -p "${outdir}"

    if prefetch "${ID}" -O "${outdir}" && \
       mv "${outdir}/${ID}/${ID}.sra" "${outdir}" 2>/dev/null && \
       rmdir "${outdir}/${ID}" 2>/dev/null && \
       (cd "${outdir}" && fastq-dump --split-3 --gzip "${ID}.sra"); then
        rm -f "${outdir}/${ID}.sra"
        echo "$(date '+%F %T') ${ID} 下载并转换成功" >> "${log}"
    else
        echo "$(date '+%F %T') ${ID} 下载或转换失败" >> "${log}"
    fi
}

export -f NCBIdownload

# outdir=/ChIP_seq_2/StemCells/data/fq
# sralist=/ChIP_seq_2/StemCells/log/download/ascp_download2.log
# nohup awk -F " " '$3 == "下载失败" {print $1,$2}' ${sralist} | parallel -j 3 --joblog /ChIP_seq_2/StemCells/log/download/parallel.log --colsep ' ' download {1} ${outdir} ${log} > /ChIP_seq_2/StemCells/log/download/ascp.log 2>&1 &