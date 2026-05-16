#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python ${SCRIPT_DIR}/run.py \
    -m /data/pub/zhousha/professional_practice/data/meta.tsv \
    -w Mutation \
    -o /data/pub/zhousha/professional_practice/output \
    -t 48 \
    --log /data/pub/zhousha/professional_practice/log/Mutation.log \
    --conda-prefix /data/pub/zhousha/env/mutation_0.1/

