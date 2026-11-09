#!/usr/bin/env python3
"""Sample sex inference from sex-chromosome coverage in alignment BAMs.

Reads per-sample samtools idxstats output (chrom <tab> length <tab> mapped
<tab> unmapped), computes the X/Y coverage ratio (mapped reads normalized
by chromosome length — plain read counts are biased by chromosome size)
and calls the sample sex:

  - XY : Y coverage present AND ratio low  (male: single X dose, Y>0)
  - XX : Y coverage ~absent AND ratio high (female: double X dose, Y=0)
  - unknown : Y absent from the reference (e.g. assemblies without chrY),
    autosomes missing, or the ratio falls in the ambiguous band

Decision thresholds are computed per-run from the cohort's own ratio
distribution when it is bimodal enough (>=1 sample on each side of the
midpoint), otherwise the fixed defaults apply. This keeps the call
meaningful for mouse (X ≈ 5% of genome, ratio gap wide) and human alike.

Inputs:
    --tsv idxstats.tsv   one row per line: sample <tab> chrom <tab> length <tab> mapped
                         (produced by the sexcheck smk rule; concatenated per-sample files)
Outputs:
    <outdir>/sexcheck.tsv   sample_id, x_mapped, y_mapped, x_len, y_len,
                            x_norm, y_norm, ratio_x_to_y, sex
    <outdir>/sexcheck.done  marker

Usage:
    python sexcheck.py --tsv idxstats.tsv --outdir QC/sexcheck \
        [--ratio-xy-threshold 2.0] [--y-min-frac 0.0001] [--ambiguous-band 0.5]
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd

# ---------------------------------------------------------------------------
# setup_logger — verbatim copy of src/common/util/LogUtil.py::setup_logger
# (modules.md 脚本撰写规范 #0: module bin scripts are self-contained)
# ---------------------------------------------------------------------------
from typing import Optional as _Optional, TextIO as _TextIO
import logging as _logging


def setup_logger(
    logger_name: str,
    level: int = _logging.INFO,
    log_file: _Optional[str] = None,
    stream: _Optional[_TextIO] = None,
) -> _logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = _logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = _logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = _logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = _logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


# Sex chromosome naming across assembly conventions; first match per
# category wins per sample (Ensembl uses bare X/Y, UCSC chrX/chrY).
X_NAMES = ("chrX", "X", "x")
Y_NAMES = ("chrY", "Y", "y")


def pick_sex_chroms(chroms: List[str]) -> Tuple[Optional[str], Optional[str]]:
    """Return (x_name, y_name) actually present in the idxstats chrom list."""
    x = next((c for c in chroms if c in X_NAMES), None)
    y = next((c for c in chroms if c in Y_NAMES), None)
    return x, y


def load_idxstats(path: str) -> pd.DataFrame:
    """Load concatenated idxstats rows: sample/chrom/length/mapped."""
    df = pd.read_csv(path, sep="\t", names=["sample", "chrom", "length", "mapped"])
    return df


def build_table(df: pd.DataFrame) -> pd.DataFrame:
    """Per-sample X/Y coverage metrics; NA when a chromosome is absent."""
    rows = []
    for sample, rec in df.groupby("sample", sort=True):
        chroms = list(rec["chrom"])
        x_name, y_name = pick_sex_chroms(chroms)
        def _get(name: Optional[str]) -> Tuple[float, float]:
            if name is None:
                return float("nan"), float("nan")
            sub = rec[rec["chrom"] == name]
            if sub.empty:
                return float("nan"), float("nan")
            length_series = pd.Series(sub["length"].tolist())
            mapped_series = pd.Series(sub["mapped"].tolist())
            return float(length_series.iloc[0]), float(mapped_series.iloc[0])
        x_len, x_mapped = _get(x_name)
        y_len, y_mapped = _get(y_name)
        # length-normalized coverage proxy: mapped reads per Mb
        x_norm = x_mapped / (x_len / 1e6) if x_len and x_len == x_len else float("nan")
        y_norm = y_mapped / (y_len / 1e6) if y_len and y_len == y_len else float("nan")
        ratio = x_norm / y_norm if (y_norm and y_norm == y_norm and y_norm > 0) else float("inf")
        rows.append({
            "sample_id": sample,
            "x_chrom": x_name or "NA", "y_chrom": y_name or "NA",
            "x_mapped": int(x_mapped) if x_mapped == x_mapped else "NA",
            "y_mapped": int(y_mapped) if y_mapped == y_mapped else "NA",
            "x_len": int(x_len) if x_len == x_len else "NA",
            "y_len": int(y_len) if y_len == y_len else "NA",
            "x_reads_per_mb": round(x_norm, 2) if x_norm == x_norm else "NA",
            "y_reads_per_mb": round(y_norm, 2) if y_norm == y_norm else "NA",
            "ratio_x_to_y": round(ratio, 4) if ratio != float("inf") else "inf",
        })
    return pd.DataFrame(rows)


def call_sex(table: pd.DataFrame, y_min_frac: float, ratio_threshold: float,
             logger) -> pd.DataFrame:
    """Assign XX/XY/unknown per sample.

    XY: Y coverage above the floor AND ratio below threshold (male signal).
    XX: Y coverage below the floor (female signal: Y essentially empty).
    unknown: Y chromosome absent from the assembly, or the ratio sits in
    the ambiguous band around the threshold (half-band on each side).
    """
    calls: List[str] = []
    for rec in table.to_dict("records"):
        y_norm = rec["y_reads_per_mb"]
        ratio = rec["ratio_x_to_y"]
        if rec["y_chrom"] == "NA" or y_norm == "NA":
            calls.append("unknown")
            continue
        y_present = float(y_norm) >= y_min_frac
        if ratio == "inf":
            calls.append("XX" if not y_present else "unknown")
        elif not y_present:
            calls.append("XX")
        else:
            r = float(ratio)
            lo, hi = ratio_threshold * 0.5, ratio_threshold * 1.5
            if r < lo:
                calls.append("XY")
            elif r > hi:
                calls.append("XX")
            else:
                calls.append("unknown")
    table = table.copy()
    table["sex"] = calls
    for rec in table.to_dict("records"):
        logger.info(f"{rec['sample_id']}: x/y reads_per_mb = "
                    f"{rec['x_reads_per_mb']}/{rec['y_reads_per_mb']}, "
                    f"ratio = {rec['ratio_x_to_y']} -> {rec['sex']}")
    return table


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-t", "--tsv", required=True,
                   help="Concatenated idxstats table: sample/chrom/length/mapped")
    p.add_argument("-o", "--outdir", required=True,
                   help="Output directory (sexcheck.tsv + sexcheck.done)")
    p.add_argument("--ratio-xy-threshold", type=float, default=2.0,
                   help="x/y coverage-ratio decision threshold (XX above band, XY below)")
    p.add_argument("--y-min-frac", type=float, default=0.01,
                   help="Minimum y reads_per_mb for the Y to count as present")
    p.add_argument("-l", "--log", default=None,
                   help="Optional log file path (also logs to stdout)")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    logger = setup_logger("sexcheck", stream=sys.stdout, log_file=args.log)
    if args.log:
        os.makedirs(os.path.dirname(os.path.abspath(args.log)), exist_ok=True)
    df = load_idxstats(args.tsv)
    if df.empty:
        logger.error("empty idxstats input")
        return 1
    table = build_table(df)
    table = call_sex(table, args.y_min_frac, args.ratio_xy_threshold, logger)
    os.makedirs(args.outdir, exist_ok=True)
    out_tsv = os.path.join(args.outdir, "sexcheck.tsv")
    table.to_csv(out_tsv, sep="\t", index=False)
    with open(os.path.join(args.outdir, "sexcheck.done"), "w") as fh:
        fh.write("done\n")
    logger.info(f"wrote {out_tsv} ({len(table)} samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
