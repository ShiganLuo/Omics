#!/usr/bin/env python3
"""Differential accessibility between sample groups across regions.

Consumes a merged quant matrix (regions x samples) and a group design, and
tests per-region differences for every contrast (Mann-Whitney U + effect
size log2FC, BH-FDR). When an annotation table is given, results are also
summarized per region class (e.g. promoter/enhancer, or CDR core vs flanking
as in Nat Genet 2025 [3]).

Design input (TSV), one row per sample:
    sample    group    contrast
    s1        ctrl     ctrl_vs_treat
    s2        treat    ctrl_vs_treat

Usage:
    python diff_accessibility.py -m percent_accessible_matrix.tsv \
        -d design.tsv -o results/diff --metric percent_accessible
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
_ROOT_DIR = _BIN_DIR.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

from src.common.util.FiberlibUtil import bh_fdr, write_tsv  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("DiffAccessibility", log_file=str(log_file))
    else:
        lg = setup_logger("DiffAccessibility")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def read_design(path: str) -> Dict[str, dict]:
    """Parse design TSV: contrast -> {group -> [samples]}."""
    contrasts: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            contrasts[f[col["contrast"]]][f[col["group"]]].append(f[col["sample"]])
    return {c: dict(g) for c, g in contrasts.items()}


def read_matrix(path: str):
    """Read regions x samples matrix; returns (header, rows as dicts)."""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, f)))
    return header, rows


def run(args: argparse.Namespace) -> int:
    from scipy.stats import mannwhitneyu

    header, rows = read_matrix(str(args.matrix))
    samples_in_matrix = header[4:]
    contrasts = read_design(str(args.design))
    logger.info(f"Matrix: {len(rows)} regions x {len(samples_in_matrix)} samples; "
                f"contrasts: {list(contrasts)}")

    class_by_key: Dict[str, str] = {}
    gene_by_key: Dict[str, str] = {}
    if args.annot is not None:
        with open(args.annot) as fh:
            ah = fh.readline().rstrip("\n").split("\t")
            acol = {name: i for i, name in enumerate(ah)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                key = f"{f[acol['chrom']]}:{f[acol['start']]}-{f[acol['end']]}|{f[acol['name']]}"
                class_by_key[key] = f[acol.get("region_class", 5)]
                if "nearest_gene" in acol:
                    gene_by_key[key] = f[acol["nearest_gene"]]

    out_tables = {}
    for contrast, groups in contrasts.items():
        if args.contrast and contrast != args.contrast:
            continue
        if len(groups) != 2:
            logger.warning(f"contrast {contrast} has {len(groups)} groups; skipping (need 2)")
            continue
        (g1, s1), (g2, s2) = sorted(groups.items())
        res_rows = []
        pvals = []
        for r in rows:
            key = f"{r['chrom']}:{r['start']}-{r['end']}|{r['name']}"
            v1 = [float(r[s]) for s in s1 if r.get(s) not in (None, "", "NA")]
            v2 = [float(r[s]) for s in s2 if r.get(s) not in (None, "", "NA")]
            if len(v1) < 2 or len(v2) < 2:
                continue
            try:
                _, p = mannwhitneyu(v1, v2, alternative="two-sided")
            except ValueError:
                p = 1.0
            mean1, mean2 = float(np.mean(v1)), float(np.mean(v2))
            l2fc = float(np.log2((mean1 + 1e-9) / (mean2 + 1e-9)))
            pvals.append(p)
            res_rows.append([
                r["chrom"], r["start"], r["end"], r["name"],
                class_by_key.get(key, "NA"),
                g1, g2, round(mean1, 6), round(mean2, 6), round(l2fc, 4),
                p,
            ])
        qvals = bh_fdr(pvals)
        for row, q in zip(res_rows, qvals):
            row.append(round(q, 6))
            row.append("yes" if q < args.alpha else "no")
        out_tables[contrast] = res_rows
        logger.info(f"{contrast}: {g1}({len(s1)}) vs {g2}({len(s2)}) -> {len(res_rows)} regions tested")

    if not out_tables:
        logger.error("no valid contrast produced results")
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    for contrast, res_rows in out_tables.items():
        out = args.outdir / f"diff_{contrast}.tsv"
        write_tsv(str(out),
                  ["chrom", "start", "end", "name", "region_class",
                   "group1", "group2", "mean_group1", "mean_group2",
                   "log2FC_group1_vs_group2", "p_value", "q_value", "significant"],
                  res_rows)
        logger.info(f"Wrote {out}")
        _gene_table(args, contrast, res_rows, gene_by_key)
        _volcano(args, contrast, res_rows)
        if args.annot is not None:
            _class_summary(args, contrast, res_rows)
    return 0


def _gene_table(args, contrast: str, res_rows, gene_by_key: Dict[str, str]) -> None:
    """Emit a gene-level table compatible with modules/function go-kegg.r.

    Output: <outdir>/<contrast>.TEcount_Gene.name.tsv with columns
    GeneName / log2FoldChange / padj (one row per gene; the most significant
    peak per gene wins). Consumable by rule function_go_kegg unchanged
    (function expects <indir>/<contrast>/<contrast>.TEcount_Gene.name.tsv, so
    the rule passes <outdir> = .../<contrast>).
    """
    best: Dict[str, List] = {}
    for r in res_rows:
        key = f"{r[0]}:{r[1]}-{r[2]}|{r[3]}"
        gene = gene_by_key.get(key) or r[3]
        if not gene or gene == ".":
            continue
        p = float(r[10])
        if gene not in best or p < float(best[gene][10]):
            best[gene] = r
    rows = [[gene, r[9], r[11]] for gene, r in best.items()]
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"{contrast}.TEcount_Gene.name.tsv"
    write_tsv(str(out), ["GeneName", "log2FoldChange", "padj"], rows)
    logger.info(f"Wrote {out} ({len(rows)} genes)")


def _volcano(args, contrast: str, res_rows) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = [r[9] for r in res_rows]
    y = [-np.log10(max(r[10], 1e-300)) for r in res_rows]
    sig = [r[12] == "yes" for r in res_rows]
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter([xi for xi, s in zip(x, sig) if not s],
               [yi for yi, s in zip(y, sig) if not s], s=6, c="#bdc3c7", label="ns")
    ax.scatter([xi for xi, s in zip(x, sig) if s],
               [yi for yi, s in zip(y, sig) if s], s=8, c="#c0392b", label="q<alpha")
    ax.set_xlabel("log2FC (group1 vs group2)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(contrast)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.outdir / f"volcano_{contrast}.pdf")
    plt.close(fig)


def _class_summary(args, contrast: str, res_rows) -> None:
    by_class = defaultdict(lambda: [0, 0])
    for r in res_rows:
        by_class[r[4]][0] += 1
        if r[12] == "yes":
            by_class[r[4]][1] += 1
    rows = [[contrast, cls, n, sig] for cls, (n, sig) in sorted(by_class.items())]
    out = args.outdir / f"class_summary_{contrast}.tsv"
    write_tsv(str(out), ["contrast", "region_class", "n_regions", "n_significant"], rows)
    logger.info(f"Wrote {out}")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="diff_accessibility.py",
        description="Differential accessibility between sample groups.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-m", "--matrix", type=Path, required=True,
                   help="merged quant matrix TSV (regions x samples)")
    p.add_argument("-d", "--design", type=Path, required=True,
                   help="design TSV with sample/group/contrast columns")
    p.add_argument("-a", "--annot", type=Path, default=None,
                   help="annotated peaks TSV (adds region_class stratification)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("-C", "--contrast", type=str, default=None,
                   help="run only this contrast (default: all contrasts in design)")
    p.add_argument("--metric", type=str, default="percent_accessible",
                   help="metric name (for logging)")
    p.add_argument("-p", "--alpha", type=float, default=0.05, help="FDR threshold")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in filter(None, [args.matrix, args.design, args.annot]):
        if not Path(f).exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would test {args.matrix} by {args.design} -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
