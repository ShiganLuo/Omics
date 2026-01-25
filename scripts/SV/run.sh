#!/bin/bash
python /data/pub/zhousha/Totipotent20251031/workflow/RNA-SNP/scripts/SV/PlaB_only.py \
    -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/DMSO_P6.sv.vcf.gz \
    -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/PlaB_P6.sv.vcf.gz \
    -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06

python /data/pub/zhousha/Totipotent20251031/workflow/RNA-SNP/scripts/SV/PlaB_only.py \
    -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/DMSO_P20.sv.vcf \
    -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/PlaB_P20.sv.vcf \
    -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20

python /data/pub/zhousha/Totipotent20251031/workflow/RNA-SNP/scripts/SV/PlaB_only.py \
    -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/DMSO_P6.sv.vcf.gz \
    -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/DMSO_P20.sv.vcf \
    -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06

python /data/pub/zhousha/Totipotent20251031/workflow/RNA-SNP/scripts/SV/PlaB_only.py \
    -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/PlaB_P6.sv.vcf.gz  \
    -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/PlaB_P20.sv.vcf \
    -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06