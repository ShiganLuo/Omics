#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python ${SCRIPT_DIR}/run.py \
    -m /home/luosg/Data/genomeStability/data/Rn7sk/meta_input.tsv \
    -w RNAseq \
    -o /home/luosg/Data/genomeStability/output \
    -t 48 \
    --log /home/luosg/Data/genomeStability/log/Rn7sk_RNAseq.log \
    --conda-prefix /home/luosg/Database/env \
    --rerun-triggers mtime \
    --Params.function.gmt /home/luosg/Data/genomeStability/workflow/Omics/assests/geneset/2C_mouse.gmt

# python ${SCRIPT_DIR}/run.py \
#     -m /home/luosg/Data/genomeStability/data/Srp54/meta_input.tsv \
#     -w ncRNAseq \
#     -o /home/luosg/Data/genomeStability/output \
#     -t 48 \
#     --log /home/luosg/Data/genomeStability/log/ncRNAseq.log \
#     --conda-prefix  ~/env \
#     --rerun-triggers mtime \
#     --Procedure.aligner star_3pass_gene
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