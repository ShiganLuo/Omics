#!/bin/bash
# conda:RAN-SNP
source workflow/scripts/download/metadata.sh
function GSMInfo(){
    sampleInfo="data/GSM.csv"
    awk -F'\t' '{print $1}' ${sampleInfo} | xargs -I {} bash -c ' 
        GSMRunInfo "{}" "data/GSMRunInfo.csv"
    '
}
export -f GSMInfo
# GSMInfo
source workflow/scripts/download/ascp.sh
function fqDownload(){
    sampleInfo="data/sampleInfo.csv"
    log="log/fqDownload.log"
    dest="data/fq"
    # awk -F'\t' '{print $3,$6}' ${sampleInfo} | while read -r SRR layout;do
    #     echo "ENAdownload ${SRR} ${layout} ${log} ${dest}"
    # done | parallel -j 4
    # ENAdownload SRR25808408 PAIRED log/fqDownload.log data/fq
    # Ensure sufficient network bandwidth,otherwise the download fastq.gz may be damaged.
    ENAdownload SRR25808406 PAIRED log/fqDownload.log data/fq
    ENAdownload SRR25808407 PAIRED log/fqDownload.log data/fq
    ENAdownload SRR25808409 PAIRED log/fqDownload.log data/fq
}
export -f fqDownload
fqDownload 
