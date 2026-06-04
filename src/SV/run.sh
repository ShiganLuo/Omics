#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# python ${SCRIPT_DIR}/exp_specific.py \
#     -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/DMSO_P6.sv.vcf.gz \
#     -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/PlaB_P6.sv.vcf.gz \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB06_vs_DMSO06 \
#     --image_format pdf

# python ${SCRIPT_DIR}/exp_specific.py \
#     -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/DMSO_P20.sv.vcf \
#     -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/PlaB_P20.sv.vcf \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20 \
#     --image_format pdf

# python ${SCRIPT_DIR}/exp_specific.py \
#     -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/DMSO_P6.sv.vcf.gz \
#     -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/DMSO_P20.sv.vcf \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/DMSO20_vs_DMSO06 \
#     --image_format pdf

# python ${SCRIPT_DIR}/exp_specific.py \
#     -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/PlaB_P6.sv.vcf.gz  \
#     -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/PlaB_P20.sv.vcf \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/PlaB20_vs_PlaB06 \
#     --image_format pdf

### circos

python ${SCRIPT_DIR}/run_circos.py \
    --vcf /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/PlaB20_vs_PlaB06_only.vcf  \
    --fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm38/GRCm38.primary_assembly.genome.fa \
    --outdir /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/circos/PlaB20_vs_PlaB06 \
    --image_formats pdf

python ${SCRIPT_DIR}/run_circos.py \
    --vcf /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/DMSO20_vs_DMSO06_only.vcf  \
    --fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm38/GRCm38.primary_assembly.genome.fa \
    --outdir /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/circos/DMSO20_vs_DMSO06 \
    --image_formats pdf

python ${SCRIPT_DIR}/run_circos.py \
    --vcf /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20_only.vcf  \
    --fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm38/GRCm38.primary_assembly.genome.fa \
    --outdir /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/circos/PlaB20_vs_DMSO20 \
    --image_formats pdf

python ${SCRIPT_DIR}/run_circos.py \
    --vcf /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB06_vs_DMSO06_only.vcf  \
    --fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm38/GRCm38.primary_assembly.genome.fa \
    --outdir /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/circos/PlaB06_vs_DMSO06 \
    --image_formats pdf

### enrichment

