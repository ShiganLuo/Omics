#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# python ${SCRIPT_DIR}/run.py \
#     -m /data/pub/zhousha/20260207_Exome/data/Exome/samplesheet.csv \
#     -w Mutation \
#     -o /data/pub/zhousha/20260207_Exome/output \
#     -t 48 \
#     --log /data/pub/zhousha/20260207_Exome/log/Mutation.log \
#     --conda-prefix /data/pub/zhousha/env/mutation_0.1/ \
#     --genome.fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa

# python ${SCRIPT_DIR}/run.py \
#     -m /data/pub/zhousha/20260207_Exome/data/RNAseq/Rawdata \
#     -w RNAseq \
#     -o /data/pub/zhousha/20260207_Exome/output \
#     -t 48 \
#     --log /data/pub/zhousha/20260207_Exome/log/RNAseq.log \
#     --conda-prefix /data/pub/zhousha/env/mutation_0.1 \
#     --genome.fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa \
#     --rerun-triggers input \
#     --dry-run

python ${SCRIPT_DIR}/run.py \
    -m /data/pub/zhousha/20260207_Exome/data/PacBio/samplesheet.csv \
    -w PacVar \
    -o /data/pub/zhousha/20260207_Exome/output \
    -t 48 \
    --log /data/pub/zhousha/20260207_Exome/log/PacVar.log \
    --conda-prefix /data/pub/zhousha/env/mutation_0.1 \
    --genome.fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa \
    --rerun-triggers input \
    --dry-run