#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

BASE = "/data/pub/zhousha/20260207_Exome/output/RNAseq/results/DEG"
OUT_DIR = "/data/pub/zhousha/20260207_Exome/output/RNAseq/results/TE_activation"
os.makedirs(OUT_DIR, exist_ok=True)

comparisons = [
    "KSR_vs_E14",
    "TLSCS_vs_E14",
    "TLSCS_vs_KSR",
]

def load_deseq(path):
    df = pd.read_csv(path, sep="\t")
    # keep standard columns
    return df

def count_sig(df, padj_col="padj", lfc_col="log2FoldChange", padj=0.05, lfc=1.0):
    df = df.copy()
    df = df[(~df[padj_col].isna()) & (~df[lfc_col].isna())]
    up = df[(df[padj_col] < padj) & (df[lfc_col] > lfc)].shape[0]
    dn = df[(df[padj_col] < padj) & (df[lfc_col] < -lfc)].shape[0]
    return up, dn

def te_activation_score(te_df, padj=0.05, lfc=0.58):
    # score = (TE_up - TE_down) / (TE_up + TE_down + 1)
    up, dn = count_sig(te_df, padj=padj, lfc=lfc)
    return (up - dn) / (up + dn + 1), up, dn

def gene_activation_score(gene_df, padj=0.05, lfc=0.58):
    # gene perturbation magnitude
    up, dn = count_sig(gene_df, padj=padj, lfc=lfc)
    return (up + dn), up, dn

rows = []
for comp in comparisons:
    te_path = os.path.join(BASE, comp, "DESeq2", "TEcount_TE.tsv")
    gene_path = os.path.join(BASE, comp, "DESeq2", "TEcount_Gene.tsv")

    te_df = load_deseq(te_path)
    gene_df = load_deseq(gene_path)

    te_score, te_up, te_dn = te_activation_score(te_df)
    gene_perturb, g_up, g_dn = gene_activation_score(gene_df)

    rows.append({
        "comparison": comp,
        "TE_up": te_up,
        "TE_down": te_dn,
        "TE_activation_score": te_score,
        "Gene_up": g_up,
        "Gene_down": g_dn,
        "Gene_perturbation": gene_perturb
    })

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT_DIR, "te_activation_summary.csv"), index=False)

# 1) Dual-axis barplot: TE up/down and gene perturbation
plt.figure(figsize=(8, 4.8))
ax = plt.gca()
x = np.arange(len(summary))
ax.bar(x - 0.15, summary["TE_up"], width=0.3, label="TE_up")
ax.bar(x + 0.15, summary["TE_down"], width=0.3, label="TE_down")
ax.set_xticks(x)
ax.set_xticklabels(summary["comparison"], rotation=15)
ax.set_ylabel("TE counts (padj<0.05, |LFC|>0.58)")
ax.legend(loc="upper right")

ax2 = ax.twinx()
ax2.plot(x, summary["Gene_perturbation"], color="black", marker="o", label="Gene_perturbation")
ax2.set_ylabel("Gene perturbation (up+down)")
plt.title("TE activation vs gene perturbation")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "te_vs_gene_perturbation.png"), dpi=300)
plt.close()

# 2) Scatter: TE activation score vs gene perturbation
plt.figure(figsize=(5.2, 4.8))
sns.scatterplot(
    data=summary,
    x="TE_activation_score",
    y="Gene_perturbation",
    s=80
)
for _, r in summary.iterrows():
    plt.text(r["TE_activation_score"], r["Gene_perturbation"], r["comparison"], fontsize=9, ha="left", va="bottom")
plt.xlabel("TE activation score (up-down)/(up+down+1)")
plt.ylabel("Gene perturbation (up+down)")
plt.title("Is gene perturbation driven by TE activation?")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "te_activation_scatter.png"), dpi=300)
plt.close()

# 3) Stacked proportions: TE_up vs TE_down proportion
prop = summary[["comparison", "TE_up", "TE_down"]].copy()
prop["total"] = prop["TE_up"] + prop["TE_down"]
prop["TE_up_frac"] = prop["TE_up"] / prop["total"].replace(0, np.nan)
prop["TE_down_frac"] = prop["TE_down"] / prop["total"].replace(0, np.nan)

plt.figure(figsize=(6.5, 4.2))
x = np.arange(len(prop))
plt.bar(x, prop["TE_up_frac"], label="TE_up_fraction")
plt.bar(x, prop["TE_down_frac"], bottom=prop["TE_up_frac"], label="TE_down_fraction")
plt.xticks(x, prop["comparison"], rotation=15)
plt.ylabel("Fraction of significant TE")
plt.title("TE activation direction by comparison")
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "te_updown_fraction.png"), dpi=300)
plt.close()

print("Done. Outputs in:", OUT_DIR)