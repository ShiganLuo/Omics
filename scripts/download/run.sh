#!/bin/bash
# source /disk5/luosg/RNAseq_DicerDHEL120260105/workflow/RNA-SNP/scripts/download/metadata.sh
# meta=/disk5/luosg/RNAseq_DicerDHEL120260105/data/README.md
# html_outdir=/disk5/luosg/RNAseq_DicerDHEL120260105/data/html
# html_log=/disk5/luosg/RNAseq_DicerDHEL120260105/log/download/GetGSMHtml.log
# awk -F'\t' '{print $1}' ${meta} | while read -r GSM;do
#     GetGSMHtml ${GSM} ${html_outdir} ${html_log}
# done
# python /disk5/luosg/RNAseq_DicerDHEL120260105/workflow/RNA-SNP/scripts/download/GSM_metadata.py \
#   --mode both \
#   --gsm-file /disk5/luosg/RNAseq_DicerDHEL120260105/data/meta.tsv \
#   --outdir /disk5/luosg/RNAseq_DicerDHEL120260105/data

python /disk5/luosg/RNAseq_DicerDHEL120260105/workflow/RNA-SNP/scripts/download/ascp_download.py \
  --meta /disk5/luosg/RNAseq_DicerDHEL120260105/data/sra_metadata.csv \
  --outdir /disk5/luosg/RNAseq_DicerDHEL120260105/data/fastq \
  --key /home/luosg/miniconda3/envs/RNA_SNP/etc/asperaweb_id_dsa.openssh \
  --log /disk5/luosg/RNAseq_DicerDHEL120260105/log/download/ascp_download.log


