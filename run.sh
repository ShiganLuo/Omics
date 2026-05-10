#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python ${SCRIPT_DIR}/run.py \
    -m /data/pub/zhousha/20260422_ClIPseq/data/meta/fastq \
    -w CLIP \
    -o /data/pub/zhousha/20260422_ClIPseq/output \
    -t 48 \
    --log /data/pub/zhousha/20260422_ClIPseq/log/CLIP.log \
    --conda-prefix /data/pub/zhousha/env/mutation_0.1/ \
    --Params.trim_galore.quality 10 \
    --Params.umi_tools.bc_pattern "NNNNNNNNNNNNNNN" \
    --Params.umi_tools.bc_pattern2 "NNNNNNNNNNNNNNN" \
    --rerun-trigger input code mtime params software-env