#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# --conda-prefix  /home/luosg/Database/env \
# --use-conda
# --snakemake-args \
# --sdm apptainer \
# --singularity-args '--bind /home/luosg/Database,/home/luosg/Data/genomeStability,/tmp'

# python ${SCRIPT_DIR}/run.py \
#     -m /home/luosg/Data/genomeStability/data/EED_public/meta_input.tsv \
#     -w RNAseq \
#     -o /home/luosg/Data/genomeStability/output/EED \
#     -t 48 \
#     --log /home/luosg/Data/genomeStability/log/EED_RNAseq.log \
#     --sdm \
#     --rerun-triggers mtime \
#     --forcerun RNAseq_generate_report \
#     --dry-run

# python ${SCRIPT_DIR}/run.py \
#     -m /home/luosg/Data/genomeStability/data/Rn7sk/meta_input.tsv \
#     -w ncRNAseq \
#     -o /home/luosg/Data/genomeStability/output/Rn7sk \
#     -t 48 \
#     --log /home/luosg/Data/genomeStability/log/ncRNAseq_Rn7sk.log \
#     --Params.workflow.aligner star \
#     --sdm

# python ${SCRIPT_DIR}/run.py \
#     -m /home/luosg/Data/genomeStability/data/20260820_RNAseq/meta_input.tsv \
#     -w ncRNAseq \
#     -o /home/luosg/Data/genomeStability/output/Hsd17b10 \
#     -t 48 \
#     --log /home/luosg/Data/genomeStability/log/ncRNAseq_Hsd17b10.log \
#     --Params.workflow.aligner star \
#     --sdm

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

# python ${SCRIPT_DIR}/run.py \
#     -m /home/luosg/Data/genomeStability/data/MS/sample_input.tsv \
#     -w QuantMS \
#     -o /home/luosg/Data/genomeStability/output/MS \
#     -t 48 \
#     --log /home/luosg/Data/genomeStability/log/MS.log \
#     --sdm

# python ${SCRIPT_DIR}/run.py \
#     -m /home/luosg/Data/genomeStability/data/Rnp/meta_input.tsv \
#     -w PeakCalling \
#     -o /home/luosg/Data/genomeStability/output/Rnp \
#     -t 48 \
#     --log /home/luosg/Data/genomeStability/log/PeakCalling.log \
#     --sdm \
#     --dry-run

python ${SCRIPT_DIR}/run.py \
    -m /home/luosg/Data/genomeStability/data/20260820_scRNAseq/meta_input.tsv \
    -w scRNAseq \
    -o /home/luosg/Data/genomeStability/output/luancao \
    -t 48 \
    --log /home/luosg/Data/genomeStability/log/scRNAseq.log \
    --sdm \
    --counters scTE cellranger \
    --aligner cellranger \
    --rerun-triggers mtime \
    --snakemake-args
    # --forcerun \
    #     scanpy_auto:tissue=Uterus,counter=scTE \
    #     scanpy_auto:tissue=Uterus,counter=cellranger
