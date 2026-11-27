#!/usr/bin/env bash
set -euo pipefail

cd /home/luosg/Data/genomeStability
export CONDA_NO_PLUGINS=true



python3 workflow/Omics/src/scRNAseq/scripts/python/utils/enrichment.py go-kegg \
  --deg-dir /home/luosg/Data/genomeStability/output/luancao/scRNAseq/analysis/results/DEG/deg_pseudobulk_gene \
  --skip-existing

python3 workflow/Omics/src/scRNAseq/scripts/python/utils/enrichment.py te-subfamily \
  --te-dir /home/luosg/Data/genomeStability/output/luancao/scRNAseq/analysis/results/DEG/deg_pseudobulk_te \
  --te-bed /home/luosg/Database/Reference/mulatta/ENSEMBL/Mmul_10/rheMac10_rmsk_TE.bed \
  --skip-existing