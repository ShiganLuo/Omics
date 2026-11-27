#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)





















python ${SCRIPT_DIR}/utils/gene_model.py \
    -t /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf \
    -g Alk -g Ciita -g Fgfr2 -g Fhit -g Lrp1b \
    -m /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20_only.vcf \
    -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/gene \
    -f pdf \
    -j 5