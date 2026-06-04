#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# source /disk5/luosg/RNAseq_DicerDHEL120260105/workflow/RNA-SNP/scripts/download/metadata.sh
# meta=/disk5/luosg/RNAseq_DicerDHEL120260105/data/README.md
# html_outdir=/disk5/luosg/RNAseq_DicerDHEL120260105/data/html
# html_log=/disk5/luosg/RNAseq_DicerDHEL120260105/log/download/GetGSMHtml.log
# awk -F'\t' '{print $1}' ${meta} | while read -r GSM;do
#     GetGSMHtml ${GSM} ${html_outdir} ${html_log}
# done
GSM_parser=${SCRIPT_DIR}/GSM_metadata.py
ASCP_downloader=${SCRIPT_DIR}/ascp_download.py
meta_input_generator=${SCRIPT_DIR}/generate_meta_input.py
ENA_ascp_key=${SCRIPT_DIR}/assests/asperaweb_id_dsa.openssh
CNGB_ascp_key=${SCRIPT_DIR}/assests/aspera01.openssh
function download_pipeline(){
    meta=$1
    outdir=$2
    log=$3
    python  ${GSM_parser}\
        --mode both \
        --gsm-file ${meta} \
        --outdir ${outdir} \
        --log ${log}
    
    # python ${ASCP_downloader} \
    #     --meta ${outdir}/sra_metadata.csv \
    #     --srr-col-name Data_id \
    #     --outdir ${outdir}/fastq \
    #     --key ${ascp_key} \
    #     --log ${log}\
    python ${meta_input_generator} \
        -i ${outdir}/sra_metadata.csv \
        -d ${outdir}/fastq \
        -o ${outdir}/meta_input.tsv
}
meat=/data/pub/zhousha/20260417_RNAseq/data/meta/meta.csv
outdir=/data/pub/zhousha/20260417_RNAseq/data/meta
log=/data/pub/zhousha/20260417_RNAseq/log/download/GSM_metadata.log
# download_pipeline ${meat} ${outdir} ${log}
function cngb_download(){
    ip=$1
    outdir=$2
    echo "ascp \
        -P 33001 \
        -v \
        -i ${CNGB_ascp_key} \
        -T \
        -D \
        -l 100m \
        -k1 \
        -d ${ip} \
        ${outdir}"
}
ip=aspera01@download.cncb.ac.cn:gsa6/CRA024880
outdir=/data/pub/zhousha/20260207_Exome/data/tRNA/fastq
cngb_download ${ip} ${outdir}
