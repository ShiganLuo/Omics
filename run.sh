#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python ${SCRIPT_DIR}/run.py \
    -m /data/pub/zhousha/20260422_ClIPseq/data/meta/fastq \
    -w CLIP \
    -o /data/pub/zhousha/20260422_ClIPseq/output \
    -t 48 \
    --log /data/pub/zhousha/20260422_ClIPseq/log/CLIP.log \
    --conda-prefix /data/pub/zhousha/env/mutation_0.1/ \
    --Params.trim_galore.quality 10