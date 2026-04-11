#!/bin/bash
# source /disk5/luosg/RNAseq_DicerDHEL120260105/workflow/RNA-SNP/scripts/download/metadata.sh
# meta=/disk5/luosg/RNAseq_DicerDHEL120260105/data/README.md
# html_outdir=/disk5/luosg/RNAseq_DicerDHEL120260105/data/html
# html_log=/disk5/luosg/RNAseq_DicerDHEL120260105/log/download/GetGSMHtml.log
# awk -F'\t' '{print $1}' ${meta} | while read -r GSM;do
#     GetGSMHtml ${GSM} ${html_outdir} ${html_log}
# done
GSM_parser=/data/pub/zhousha/20260411_MERIPseq/workflow/RNA-SNP/src/download/GSM_metadata.py
ASCP_downloader=/data/pub/zhousha/20260411_MERIPseq/workflow/RNA-SNP/src/download/ascp_download.py
ascp_key=/home/zhousha/miniforge3/envs/RNA/etc/asperaweb_id_dsa.openssh
function download_pipeline(){
    meta=$1
    outdir=$2
    log=$3
    # python  ${GSM_parser}\
    #     --mode both \
    #     --gsm-file ${meta} \
    #     --outdir ${outdir} \
    
    python ${ASCP_downloader} \
        --meta ${outdir}/sra_metadata.csv \
        --srr-col-name Data_id \
        --outdir ${outdir}/fastq \
        --key ${ascp_key} \
        --log ${log}
}
meat=/data/pub/zhousha/20260411_MERIPseq/data/meta.tsv
outdir=/data/pub/zhousha/20260411_MERIPseq/data/meta
log=/data/pub/zhousha/20260411_MERIPseq/log/download/ascp_download.log
download_pipeline ${meat} ${outdir} ${log}



