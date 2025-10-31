#!/bin/bash
## check whether fq has been download already or not

function isFqExist(){
    local sraNumber=$1
    local fqDir=$2
    local log=$3
    cd ${fqDir}
    if [[ (-f ${sraNumber}_1.fastq.gz  && -f ${sraNumber}_2.fastq.gz)  ||  -f ${sraNumber}.fastq.gz ]]; then
        echo -e "${sraNumber}\tfq文件存在" >> ${log}
    else
        echo -e "${sraNumber}\tfq文件不存在" >> ${log}
    fi
}
export -f isFqExist

function gzipTest(){
    local dir=$1
    local log=$2

    for f in ${1}*.gz; do
        gzip -t "$f" >/dev/null 2>&1 \
        && echo "OK:    $f" \
        || echo "ERROR: $f"
    done > ${log} 2>&1
}
export -f gzipTest

function gzipTestForSra(){
    local fqDir=$1
    local sraNumber=$2
    local log=$3
    for f in ${fqDir}/${sraNumber}*.gz; do
        gzip -t "$f" >/dev/null 2>&1 \
        && echo "OK:    $f" \
        || echo "ERROR: $f"
    done >> ${log} 2>&1
}

export -f gzipTestForSra