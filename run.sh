#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# python ${SCRIPT_DIR}/run.py \
#     -m /rna_seq_1/luoshg/Chipseq_20260709/data/meta.tsv \
#     -w PeakCalling \
#     -o /rna_seq_1/luoshg/Chipseq_20260709/output \
#     -t 48 \
#     --log /rna_seq_1/luoshg/Chipseq_20260709/logs/PeakCalling.log \
#     --conda-prefix /rna_seq_1/luoshg/env \
#     --Params.macs3.cutoff_analysis \
#     --rerun-triggers mtime

python ${SCRIPT_DIR}/run.py \
    -m /rna_seq_1/luoshg/Chipseq_20260709/data/Srp54/meta_input.tsv \
    -w ncRNAseq \
    -o /rna_seq_1/luoshg/Chipseq_20260709/output \
    -t 48 \
    --log /rna_seq_1/luoshg/Chipseq_20260709/logs/ncRNAseq.log \
    --conda-prefix /rna_seq_1/luoshg/env \
    --rerun-triggers input \
    --Procedure.aligner star_3pass \
    --dry-run

# python ${SCRIPT_DIR}/run.py \
#     -m /data/pub/zhousha/20260207_Exome/data/tRNA/meta.tsv \
#     -w tRNAseq \
#     -o /data/pub/zhousha/20260207_Exome/output \
#     -t 48 \
#     --log /data/pub/zhousha/20260207_Exome/log/tRNAseq.log \
#     --conda-prefix /data/pub/zhousha/env/mutation_0.1 \
#     --genome.fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa \
#     --rerun-triggers mtime \
#     --Params.cutadapt.match_read_wildcards True \
#     --Params.cutadapt.cut 16 \
#     --Params.cutadapt.trimmed_only True \
#     --Params.cutadapt.adapter_r1 AAAAAAAACAAAAAAAAAA AAAA$ AAA$ AA$ A$ \
#     --Params.cutadapt.minimum_length 50 \
#     --Params.cutadapt.maximum_length 110 \


# python ${SCRIPT_DIR}/run.py \
#     -m /data/pub/zhousha/20260207_Exome/data/PacBio/samplesheet.csv\
#     -w PacVar \
#     -o /data/pub/zhousha/20260207_Exome/output \
#     -t 48 \
#     --log /data/pub/zhousha/20260207_Exome/log/PacVar.log \
#     --conda-prefix /data/pub/zhousha/env/mutation_0.1/ \
#     --genome.fasta /data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa \
#     --snakemake-args \
#     --sdm apptainer \
#     --singularity-args '--bind /data/pub/zhousha/Reference'