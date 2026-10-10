#!/usr/bin/env python3
"""Summarize a DIA-NN report into per-run QC and a protein matrix.

Reads the main DIA-NN report (report.tsv) and produces:
  - diann_qc.tsv: per-run precursor/protein identifications at q-value
    threshold, median q-value, and PG.MaxLFQ intensity coverage
  - diann_protein_matrix.tsv: Protein.Group x Run matrix of PG.MaxLFQ
    intensities (q-filtered precursors only)

Usage:
    python diann_summary.py -i diann.report.tsv -o outdir
    python diann_summary.py -i diann.report.tsv -o outdir --qvalue 0.01 -l log.txt
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

_SRC_DIR = Path(__file__).resolve().parent.parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("DiaNNSummary", log_file=str(log_file))
    else:
        lg = setup_logger("DiaNNSummary")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def run(args: argparse.Namespace) -> int:
    logger.info(f"Loading report {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    logger.info(f"Report rows: {len(df)}, columns: {len(df.columns)}")

    # DIA-NN column names differ slightly across versions; resolve defensively.
    def col(*candidates: str) -> Optional[str]:
        for c in candidates:
            if c in df.columns:
                return c
        return None

    run_col = col("Run", "File.Name", "Run_ID")
    q_col = col("Global.Q.Value", "Q.Value", "Lib.Q.Value")
    pg_col = col("Protein.Group", "Protein.Ids")
    maxlfq_col = col("PG.MaxLFQ", "PG.Quantity", "Genes.Quantity")
    if run_col is None or pg_col is None:
        logger.error(f"Required columns missing (need Run and Protein.Group); "
                     f"got: {list(df.columns)[:20]}...")
        return 1

    qval = pd.to_numeric(df[q_col], errors="coerce") if q_col else pd.Series(1.0, index=df.index)
    keep = qval <= args.qvalue
    dfq = df[keep].copy()
    logger.info(f"q<={args.qvalue}: {len(dfq)} precursor rows")

    # Per-run QC
    qc_rows: List[List] = []
    for run_name, sub in dfq.groupby(run_col):
        q_sub = qval.loc[sub.index]
        intensity = pd.to_numeric(sub[maxlfq_col], errors="coerce") if maxlfq_col else pd.Series(dtype=float)
        qc_rows.append([
            run_name,
            sub["Precursor.Id"].nunique() if "Precursor.Id" in sub else len(sub),
            sub[pg_col].nunique(),
            round(float(q_sub.median()), 6) if len(q_sub) else "NA",
            int((intensity > 0).sum()) if len(intensity) else 0,
        ])
    qc = pd.DataFrame(qc_rows, columns=[
        "run", "precursors_q", "protein_groups_q", "median_q_value", "pg_maxlfq_nonzero"])
    qc_path = args.outdir / "diann_qc.tsv"
    qc.to_csv(qc_path, sep="\t", index=False)
    logger.info(f"Wrote {qc_path}")

    # Protein matrix (PG.MaxLFQ, q-filtered)
    if maxlfq_col:
        mat = (dfq.assign(_int=pd.to_numeric(dfq[maxlfq_col], errors="coerce"))
               .pivot_table(index=pg_col, columns=run_col, values="_int",
                            aggfunc="max", fill_value=0))
        mat_path = args.outdir / "diann_protein_matrix.tsv"
        mat.to_csv(mat_path, sep="\t")
        logger.info(f"Wrote {mat_path} ({mat.shape[0]} proteins x {mat.shape[1]} runs)")
    else:
        logger.warning("No PG.MaxLFQ column found; protein matrix skipped")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="diann_summary.py",
        description="Per-run QC and protein matrix from a DIA-NN report.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", type=Path, required=True,
                   help="DIA-NN report.tsv")
    p.add_argument("-o", "--outdir", type=Path, required=True,
                   help="output directory")
    p.add_argument("-q", "--qvalue", type=float, default=0.01,
                   help="q-value threshold for identifications")
    p.add_argument("-l", "--log", type=Path, default=None,
                   help="log file path (default: log to stdout)")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    if not args.input.exists():
        logger.error(f"input not found: {args.input}")
        return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would summarize {args.input} -> {args.outdir}")
        return 0
    args.outdir.mkdir(parents=True, exist_ok=True)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
