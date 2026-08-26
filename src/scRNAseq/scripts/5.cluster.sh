#!/bin/bash
nohup python scripts/python/cellranger/1.cluster.py --options 2 --input /home/lsg/Data/glioblastoma/output/h5ad/GBM.h5ad --out GBM > log/py.1.clusteGBM.log 2>&1 &