#!/usr/bin/env bash
set -euo pipefail

cd /home/luosg/Data/genomeStability
export CONDA_NO_PLUGINS=true

# ── 完整调用命令（供审核） ──
# --deg-mode gene te all: 三种模式分别分析
#   gene → 排除 TE，只用蛋白编码基因等
#   te   → 只用 TE（gene_type == "TE"）
#   all  → 基因 + TE 全部（原始默认行为）
#
# 输出目录：
#   pseudobulk_counts / deg_pseudobulk          (all mode)
#   pseudobulk_counts_gene / deg_pseudobulk_gene (gene mode)
#   pseudobulk_counts_te / deg_pseudobulk_te     (te mode)
#
# --skip-plot: 火山图由 generate_report.py 单独生成

# python workflow/Omics/src/scRNAseq/scripts/python/utils/pseudobulk.py \
#   --out-dir /home/luosg/Data/genomeStability/output/luancao/scRNAseq/analysis/results/DEG \
#   --h5ad ovaries=output/luancao/scRNAseq/common/5_combine_h5ad/ovaries/ovaries_scTE_auto.h5ad \
#   --h5ad uterus=output/luancao/scRNAseq/common/5_combine_h5ad/Uterus/Uterus_scTE_auto.h5ad \
#   --compare ovaries:luanchao-21310-10XSC3:luanchao-11238-10XSC3=Youth_vs_Aged \
#   --compare uterus:zigong-21310-10XSC3:ZIGONG-21224-10XSC3=Sham_vs_OV-POI \
#   --compare uterus:zigong-11238-10XSC3:zigong-1411200-10XSC3=Baseline_vs_OV \
#   --compare uterus:zigong-1411200-10XSC3:zigong-2200671-10XSC3=OV_vs_Transplant \
#   --deg-mode gene te all

# ── GO/KEGG enrichment for gene DEGs ──
python3 workflow/Omics/src/scRNAseq/scripts/python/utils/enrichment.py go-kegg \
  --deg-dir /home/luosg/Data/genomeStability/output/luancao/scRNAseq/analysis/results/DEG/deg_pseudobulk_gene \
  --skip-existing

# ── TE subfamily analysis for TE DEGs ──
python3 workflow/Omics/src/scRNAseq/scripts/python/utils/enrichment.py te-subfamily \
  --te-dir /home/luosg/Data/genomeStability/output/luancao/scRNAseq/analysis/results/DEG/deg_pseudobulk_te \
  --te-bed /home/luosg/Database/Reference/mulatta/ENSEMBL/Mmul_10/rheMac10_rmsk_TE.bed \
  --skip-existing