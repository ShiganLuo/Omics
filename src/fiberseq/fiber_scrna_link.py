#!/usr/bin/env python3
"""Link Fiber-seq chromatin state with single-cell RNA-seq (scTE) profiles.

Not a single-cell pairing (different modalities, different cells): a
co-location association across features.

Layer 1 - TE families (main deliverable):
  Fiber side: per TE family/subfamily -> median MSP accessibility, mCpG and
  m6A density (from te_landscape.py's te_copy_metrics.tsv or computed here
  from element BEDs + TE bed).
  RNA side: per TE family/subfamily -> mean scTE expression per cell type
  from an scTE h5ad (raw matrix, var gene_type == 'TE' with family_id/class_id).
  Output: joined table + per-cell-type scatter (x = chromatin metric,
  y = expression) with Spearman rho.

Layer 2 - genes (optional):
  Fiber side: promoter FIRE signal per gene (MSP coverage over TSS+/-1kb from
  msp.bed.gz, gene models from GTF via --gtf).
  RNA side: mean expression per cell type (gene-level part of the same h5ad).
  Same scatter/join.

Run inside an env with anndata/scanpy (e.g. scanpy.sif).

Usage:
  python fiber_scrna_link.py --h5ad ovaries_scTE_auto.h5ad \
      --te-metrics te_copy_metrics.tsv -o link_out/ [--gtf genes.gtf \
      --msp sample.msp.bed.gz --sample-label luanchao-21310-10XSC3]
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

ANNDATA_MISSING = False
try:
    import anndata as ad
    import scipy.sparse as sp
except ImportError:  # pragma: no cover
    ANNDATA_MISSING = True


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("FiberScRnaLink", log_file=str(log_file))
    else:
        lg = setup_logger("FiberScRnaLink")
    for h in lg.handlers:
        if type(h) is logging.StreamHandler:
            h.setStream(sys.stdout)
    return lg


logger = _get_logger()


# --------------------------------------------------------------------------- #
# RNA side
# --------------------------------------------------------------------------- #
def te_expression_per_celltype(h5ad_path: str, level: str = "family",
                               layer: str = "counts") -> pd.DataFrame:
    """Mean TE expression (raw counts CPM-ish, log1p) per cell type x TE key."""
    a = ad.read_h5ad(h5ad_path)
    var = a.raw.var if a.raw is not None else a.var
    X = a.raw.X if a.raw is not None else a.X
    mask = var["gene_type"].astype(str) == "TE"
    key_col = {"family": "family_id", "class": "class_id",
               "subfamily": None}.get(level)
    if key_col is None:  # subfamily = var index
        keys = pd.Index([str(v) for v in var.index])
    else:
        keys = var[key_col].astype(str)
    keys = pd.Series(keys, index=var.index)[mask]
    if hasattr(X, "toarray"):
        X = X.toarray()
    df = pd.DataFrame(np.asarray(X), index=[str(c) for c in a.obs_names],
                      columns=[str(v) for v in var.index])
    cell_types = a.obs["cell_type"].astype(str).values
    out = {}
    grouped = df.groupby(cell_types, observed=True)
    for ct, block in grouped:
        sub = block.loc[:, keys.index]
        # aggregate features -> key (sum counts), then log1p mean per cell
        agg = sub.T.groupby(keys.values).sum().T
        out[ct] = np.log1p(agg.mean(axis=0))
    table = pd.DataFrame(out)
    return table


def gene_expression_per_celltype(h5ad_path: str, top_n: int = 5000) -> pd.DataFrame:
    """Mean log1p gene expression per cell type (gene-level features)."""
    a = ad.read_h5ad(h5ad_path)
    var = a.raw.var if a.raw is not None else a.var
    X = a.raw.X if a.raw is not None else a.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    genes = [str(v) for v in var.index]
    df = pd.DataFrame(np.asarray(X), index=[str(c) for c in a.obs_names],
                      columns=genes)
    # cap to genes with any signal
    keep = df.columns[(df > 0).sum(axis=0) > df.shape[0] * 0.01]
    df = df[keep]
    cell_types = a.obs["cell_type"].astype(str).values
    return df.groupby(cell_types, observed=True).mean()


# --------------------------------------------------------------------------- #
# Fiber side
# --------------------------------------------------------------------------- #
def fiber_te_metrics(te_metrics_tsv: str, sample_label: Optional[str]) -> pd.DataFrame:
    """Median per-TE-key chromatin metrics from te_landscape output."""
    df = pd.read_csv(te_metrics_tsv, sep="\t")
    if sample_label and "sample" in df.columns:
        df = df[df["sample"] == sample_label]
    return df


def promoter_msp_coverage(msp_path: str, gtf_path: str, flank: int = 1000,
                          max_genes: int = 20000) -> pd.DataFrame:
    """Mean MSP coverage over TSS +/- flank per gene (streamed join)."""
    tss: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    import gzip
    op = gzip.open if str(msp_path).endswith(".gz") else open
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            attrs = dict(kv.strip().split(" ", 1)
                         for kv in f[8].split(";") if " " in kv)
            gid = attrs.get("gene_id", "").strip('"')
            gname = attrs.get("gene_name", gid).strip('"')
            start = int(f[3]) - 1
            tss[f[0]].append((start if f[6] == "+" else int(f[4]), gname))
    for c in tss:
        tss[c].sort()
    cov: Dict[str, List[float]] = defaultdict(list)
    with op(msp_path, "rt") as fh:
        pending = {c: 0 for c in tss}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] not in tss:
                continue
            s, e = int(f[1]), int(f[2])
            if len(f) > 11:  # BED12: use blocks
                sizes = [int(x) for x in f[10].split(",") if x]
                starts = [int(x) for x in f[11].split(",") if x]
                blocks = [(s + b, s + b + sz) for b, sz in zip(starts, sizes)]
            else:
                blocks = [(s, e)]
            i = pending[f[0]]
            anchors = tss[f[0]]
            while i < len(anchors) and anchors[i][0] < e:
                apos, gname = anchors[i]
                for bs, be in blocks:
                    s0, e0 = max(bs, apos - flank), min(be, apos + flank)
                    if e0 > s0:
                        cov[gname].append((e0 - s0) / (2 * flank))
                i += 1
            pending[f[0]] = max(pending[f[0]], i - 1 if i < len(anchors) else 0)
    return pd.DataFrame({"gene": list(cov.keys()),
                         "msp_promoter_cov": [float(np.mean(v)) for v in cov.values()]})


# --------------------------------------------------------------------------- #
# association
# --------------------------------------------------------------------------- #
def join_and_test(fiber: pd.DataFrame, rna: pd.DataFrame,
                  fiber_key: str, metrics: Sequence[str]) -> pd.DataFrame:
    """Join on key and Spearman-correlate each fiber metric with each cell type."""
    from scipy.stats import spearmanr

    rows = []
    for metric in metrics:
        if metric not in fiber.columns:
            continue
        f = fiber[[fiber_key, metric]].groupby(fiber_key)[metric].median()
        shared = f.index.intersection(rna.index)
        if len(shared) < 5:
            continue
        f = f.loc[shared]
        for ct in rna.columns:
            r = rna.loc[shared, ct]
            rho, p = spearmanr(f.values, r.values)
            rows.append({"metric": metric, "cell_type": ct,
                         "n_families": len(shared), "spearman_rho": rho,
                         "spearman_p": p})
    return pd.DataFrame(rows)


def plot_scatter(fiber: pd.DataFrame, rna: pd.DataFrame, fiber_key: str,
                 metric: str, outdir: Path, formats: Sequence[str]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    f = fiber[[fiber_key, metric]].groupby(fiber_key)[metric].median()
    shared = f.index.intersection(rna.index)
    if len(shared) < 5:
        return
    f = f.loc[shared]
    n = len(rna.columns)
    fig, axes = plt.subplots(1, n, figsize=(3 * n, 3), squeeze=False)
    for ax, ct in zip(axes[0], rna.columns):
        ax.scatter(f.values, rna.loc[shared, ct].values, s=12, alpha=0.7)
        ax.set_xlabel(metric, fontsize=8)
        ax.set_ylabel(f"log1p expr", fontsize=8)
        ax.set_title(ct, fontsize=8)
    fig.suptitle(f"TE {fiber_key}: fiber chromatin vs scRNA expression")
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(outdir / f"te_{fiber_key}_{metric}_scatter.{fmt}", dpi=250)
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="fiber_scrna_link.py",
        description="Co-location association of Fiber-seq chromatin and scRNA/scTE.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--h5ad", required=True, help="scTE auto h5ad (cell_type obs)")
    p.add_argument("--te-metrics", required=True,
                   help="te_landscape te_copy_metrics.tsv (fiber side)")
    p.add_argument("--te-level", default="family", choices=["class", "family", "subfamily"])
    p.add_argument("--fiber-sample", default=None,
                   help="restrict fiber side to this sample label")
    p.add_argument("--gtf", default=None, help="gene GTF for promoter layer")
    p.add_argument("--msp", default=None, help="msp.bed.gz for promoter layer")
    p.add_argument("--promoter-flank", type=int, default=1000)
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("-f", "--formats", nargs="+", default=["png"],
                   choices=["png", "pdf", "svg", "eps"])
    p.add_argument("-l", "--log", type=Path, default=None)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    if ANNDATA_MISSING:
        logger.error("anndata not importable -- run inside scanpy.sif or an env with anndata")
        return 2
    args.outdir.mkdir(parents=True, exist_ok=True)

    # ---- TE layer
    fiber = fiber_te_metrics(args.te_metrics, args.fiber_sample)
    key = {"class": "class", "family": "family", "subfamily": "subfamily"}[args.te_level]
    rna = te_expression_per_celltype(args.h5ad, level=args.te_level)
    rna.to_csv(args.outdir / "te_expression_per_celltype.tsv", sep="\t")
    assoc = join_and_test(fiber, rna, key,
                          ["msp_frac", "nuc_frac", "mcpg_per_kb_overlap",
                           "m6a_per_kb_overlap"])
    assoc.to_csv(args.outdir / "te_chromatin_vs_expression.tsv", sep="\t", index=False)
    if not args.no_plot and "msp_frac" in fiber.columns:
        plot_scatter(fiber, rna, key, "msp_frac", args.outdir, args.formats)
    logger.info(f"TE layer: {len(assoc)} associations")

    # ---- gene promoter layer (optional)
    if args.gtf and args.msp:
        prom = promoter_msp_coverage(args.msp, args.gtf, args.promoter_flank)
        gex = gene_expression_per_celltype(args.h5ad)
        merged = prom.merge(gex.T.reset_index().rename(columns={"index": "gene"}),
                            on="gene", how="inner")
        merged.to_csv(args.outdir / "gene_promoter_link.tsv", sep="\t", index=False)
        logger.info(f"gene layer: {len(merged)} genes joined")
    logger.info(f"done -> {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
