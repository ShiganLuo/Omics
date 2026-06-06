#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

### exp specific SV analysis

# python ${SCRIPT_DIR}/exp_specific.py \
#     -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/DMSO_P6.sv.vcf.gz \
#     -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/PlaB_P6.sv.vcf.gz \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB06_vs_DMSO06 \
#     -f pdf

# python ${SCRIPT_DIR}/exp_specific.py \
#     -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/DMSO_P20.sv.vcf \
#     -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/PlaB_P20.sv.vcf \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20 \
#     -f pdf

# python ${SCRIPT_DIR}/exp_specific.py \
#     -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/DMSO_P6.sv.vcf.gz \
#     -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/DMSO_P20.sv.vcf \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/DMSO20_vs_DMSO06 \
#     -f pdf

# python ${SCRIPT_DIR}/exp_specific.py \
#     -c /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/PlaB_P6.sv.vcf.gz  \
#     -e /data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/PlaB_P20.sv.vcf \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/PlaB20_vs_PlaB06 \
#     -f pdf

### pbsv diff analysis
# python ${SCRIPT_DIR}/pbsv_sv_diff_analysis.py \
#     -g PlaB06_vs_DMSO06:/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB06_vs_DMSO06_only.vcf \
#     -g PlaB20_vs_PlaB06:/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/PlaB20_vs_PlaB06_only.vcf \
#     -g DMSO20_vs_DMSO06:/data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/DMSO20_vs_DMSO06_only.vcf \
#     -g PlaB20_vs_DMSO20:/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20_only.vcf \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/pbsv_specific_svtype_diff \
#     -f pdf

### circos

# python ${SCRIPT_DIR}/run_circos.py \
#     --vcf /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/PlaB20_vs_PlaB06_only.vcf  \
#     --fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm38/GRCm38.primary_assembly.genome.fa \
#     --outdir /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/circos/PlaB20_vs_PlaB06 \
#     -f pdf

# python ${SCRIPT_DIR}/run_circos.py \
#     --vcf /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/DMSO20_vs_DMSO06_only.vcf  \
#     --fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm38/GRCm38.primary_assembly.genome.fa \
#     --outdir /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/circos/DMSO20_vs_DMSO06 \
#     -f pdf

# python ${SCRIPT_DIR}/run_circos.py \
#     --vcf /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20_only.vcf  \
#     --fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm38/GRCm38.primary_assembly.genome.fa \
#     --outdir /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/circos/PlaB20_vs_DMSO20 \
#     -f pdf

# python ${SCRIPT_DIR}/run_circos.py \
#     --vcf /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB06_vs_DMSO06_only.vcf  \
#     --fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm38/GRCm38.primary_assembly.genome.fa \
#     --outdir /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/circos/PlaB06_vs_DMSO06 \
#     -f pdf

# ### enrichment

# python ${SCRIPT_DIR}/run_enrichment.py \
#     -a /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB06_vs_DMSO06_annotated.tab \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/enrichment \
#     -g /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf \
#     -fi /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa.fai \
#     -f pdf

# python ${SCRIPT_DIR}/run_enrichment.py \
#     -a /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20_annotated.tab \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/enrichment \
#     -g /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf \
#     -fi /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa.fai \
#     -f pdf

# python ${SCRIPT_DIR}/run_enrichment.py \
#     -a /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/DMSO20_vs_DMSO06_annotated.tab \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/enrichment \
#     -g /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf \
#     -fi /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa.fai \
#     -f pdf

# python ${SCRIPT_DIR}/run_enrichment.py \
#     -a /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/PlaB20_vs_PlaB06_annotated.tab \
#     -o /data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/enrichment \
#     -g /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf \
#     -fi /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa.fai \
#     -f pdf

### OncoPrint

python ${SCRIPT_DIR}/run_OncoPrint.py \
    -g PlaB06_vs_DMSO06:/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB06_vs_DMSO06_annotated.tab \
    -g PlaB20_vs_PlaB06:/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/PlaB20_vs_PlaB06_annotated.tab \
    -g DMSO20_vs_DMSO06:/data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/DMSO20_vs_DMSO06_annotated.tab \
    -g PlaB20_vs_DMSO20:/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20_annotated.tab \
    -d PlaB06:/data/pub/zhousha/Totipotent20251031/data/Pacbio/RNAseq/P6.tsv \
    -d PlaB20:/data/pub/zhousha/Totipotent20251031/data/Pacbio/RNAseq/P20.tsv \
    --cosmic_file /data/pub/zhousha/Reference/human/Cosmic_CancerGeneCensusHallmarksOfCancer_v103_GRCh38.tsv \
    -o /data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint \
    -f pdf
