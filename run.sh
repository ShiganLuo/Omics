#!/bin/bash
python /data/pub/zhousha/20260422_RNAseq/workflow/RNA-SNP/run.py \
    -m /data/pub/zhousha/20260422_RNAseq/data/meta.tsv \
    -w RNAseq -o /data/pub/zhousha/20260422_RNAseq/output \
    -t 48 \
    --log /data/pub/zhousha/20260422_RNAseq/log/RNAseq.log \
    --conda-prefix /data/pub/zhousha/env/mutation_0.1/