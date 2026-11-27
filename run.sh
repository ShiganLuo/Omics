#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)






python ${SCRIPT_DIR}/run.py \
    -m /data/16T/luosg/data/meta_DR_Totipotency.tsv \
    -w DR \
    -o /data/16T/luosg/output/Totipotency \
    --log /data/16T/luosg/log/DR_Totipotency.log \
    --sdm \
    --rerun-triggers mtime









