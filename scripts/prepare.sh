#!/bin/bash
# conda:RAN-SNP
#########################prepare for fq###########################
source workflow/scripts/download/metadata.sh

function GSMInfo(){
    sampleInfo=$1
    awk -F'\t' '{print $1}' ${sampleInfo} | xargs -I {} bash -c ' 
        GSMRunInfo "{}" "data/GSMRunInfo.csv"
    '
}
export -f GSMInfo
# GSMInfo
source workflow/scripts/download/ascp.sh
function fqDownload(){
    # sampleInfo="data/sampleInfo.csv"
    # log="log/fqDownload.log"
    # dest="data/fq"
    # awk -F'\t' '{print $3,$6}' ${sampleInfo} | while read -r SRR layout;do
    #     echo "ENAdownload ${SRR} ${layout} ${log} ${dest}"
    # done | parallel -j 4
    # ENAdownload SRR25808408 PAIRED log/fqDownload.log data/fq
    # Ensure sufficient network bandwidth,otherwise the download fastq.gz may be damaged.
    # ENAdownload SRR15731262 PAIRED log/fqDownload.log data/GSE183522/fq
    # ENAdownload SRR16117863 PAIRED log/fqDownload.log data/GSE185005/fq
    # ENAdownload SRR16117860 PAIRED log/fqDownload.log data/GSE185005/fq
    # ENAdownload SRR16117857 PAIRED log/fqDownload.log data/GSE185005/fq
    # ENAdownload SRR23635563 PAIRED log/fqDownload.log data/GSE224794/fq
    # ENAdownload SRR18181017 PAIRED log/fqDownload.log data/GSE204801/fq
    # ENAdownload SRR18181018 PAIRED log/fqDownload.log data/GSE204801/fq
    # ENAdownload SRR18181020 PAIRED log/fqDownload.log data/GSE204801/fq
    # ENAdownload SRR13633379 PAIRED log/fqDownload.log data/GSE166216/fq
    # ENAdownload SRR13633380 PAIRED log/fqDownload.log data/GSE166216/fq
    # ENAdownload SRR13633381 PAIRED log/fqDownload.log data/GSE166216/fq
    # ENAdownload SRR13633383 PAIRED log/fqDownload.log data/GSE166216/fq
    ENAdownload SRR25690002 PAIRED log/fqDownload.log data/GSE224794/fq
    ENAdownload SRR25690001 PAIRED log/fqDownload.log data/GSE224794/fq
    

}
export -f fqDownload
# fqDownload
#########################prepare for annovar###########################
######human
# # snakemake -s workflow/snakemake/RNASNP-human.smk \
# #     --cores 30 --directory workflow/snakemake \
# #     --config indir="../../data/GSE224794" outdir="../../output/GSE224794" \
# #     --use-conda --dry-run
# nohup time snakemake -s workflow/snakemake/RNASNP-human.smk \
#     --cores 30 --directory workflow/snakemake \
#     --config indir="../../data/GSE224794" outdir="../../output/GSE224794" \
#     --use-conda >> log/GSE224794.log 2>&1 &

# # snakemake -s workflow/snakemake/RNASNP-human.smk \
# #     --cores 30 --directory workflow/snakemake \
# #     --config indir="../../data/GSE204801" outdir="../../output/GSE204801" \
# #     --use-conda --dry-run
# nohup time snakemake -s workflow/snakemake/RNASNP-human.smk \
#     --cores 30 --directory workflow/snakemake \
#     --config indir="../../data/GSE204801" outdir="../../output/GSE204801" \
#     --use-conda >> log/GSE204801.log 2>&1 &
######mouse
# snakemake -s workflow/snakemake/RNASNP-noRemove.smk \
#     --cores 30 --directory workflow/snakemake \
#     --config indir="../../data/GSE166216" outdir="../../output/GSE166216" \
#     --use-conda --dry-run
# nohup time snakemake -s workflow/snakemake/RNASNP-noRemove.smk \
#     --cores 30 --directory workflow/snakemake \
#     --config indir="../../data/GSE166216" outdir="../../output/GSE166216" \
#     --use-conda >> log/GSE166216.log 2>&1 &

# snakemake -s workflow/snakemake/RNASNP-noRemove.smk \
#     --cores 30 --directory workflow/snakemake \
#     --config indir="../../data/GSE183522" outdir="../../output/GSE183522" \
#     --use-conda --dry-run
nohup time snakemake -s workflow/snakemake/RNASNP-noRemove.smk \
    --cores 30 --directory workflow/snakemake \
    --config indir="../../data/GSE183522" outdir="../../output/GSE183522" \
    --use-conda >> log/GSE183522.log 2>&1 &

# snakemake -s workflow/snakemake/RNASNP-noRemove.smk \
#     --cores 30 --directory workflow/snakemake \
#     --config indir="../../data/GSE185005" outdir="../../output/GSE185005" \
#     --use-conda --dry-run
nohup time snakemake -s workflow/snakemake/RNASNP-noRemove.smk \
    --cores 30 --directory workflow/snakemake \
    --config indir="../../data/GSE185005" outdir="../../output/GSE185005" \
    --use-conda >> log/GSE185005.log 2>&1 &