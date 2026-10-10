#!/usr/bin/env python3
"""Quantify per-region Fiber-seq metrics and merge them into a matrix.

Subcommand `run`   : per-sample region quantification (percent accessibility,
                     nucleosome occupancy, m6A density, FIRE fraction).
Subcommand `merge` : combine per-sample region_quant.tsv files into one
                     matrix (regions x samples) per metric.

Usage:
    python quant_matrix.py run -p sample.fire_peaks.bed \
        --msp sample.msp.bed.gz --nuc sample.nuc.bed.gz --m6a sample.m6a.bed.gz \
        --fire sample.fire.bed.gz -o results/sample/quant --sample sample
    python quant_matrix.py merge -q s1.region_quant.tsv,s2.region_quant.tsv \
        -o results/matrix
"""
from __future__ import annotations

import argparse
import logging
import sys
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

from src.common.util.FiberlibUtil import RegionIndex, read_bed, safe_div, write_tsv  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("QuantMatrix", log_file=str(log_file))
    else:
        lg = setup_logger("QuantMatrix")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()

REGION_KEY = ("chrom", "start", "end", "name")


def region_key(row: Sequence) -> str:
    return f"{row[0]}:{row[1]}-{row[2]}|{row[3]}"


def fibers_overlapping(index: RegionIndex, records, chrom, start, end) -> set:
    """Set of fiber names with any element overlapping the region."""
    return {records[i].name for i in index.query(chrom, start, end)}


def run_sample(args: argparse.Namespace) -> int:
    peaks = read_bed(str(args.peaks))
    msp = read_bed(str(args.msp))
    nuc = read_bed(str(args.nuc))
    m6a = read_bed(str(args.m6a))
    fire = read_bed(str(args.fire)) if args.fire else []
    msp_i, nuc_i, m6a_i = RegionIndex(msp), RegionIndex(nuc), RegionIndex(m6a)
    fire_i = RegionIndex(fire) if fire else None
    logger.info(f"peaks={len(peaks)} msp={len(msp)} nuc={len(nuc)} m6a={len(m6a)} fire={len(fire)}")

    rows = []
    for p in peaks:
        span = max(1, p.end - p.start)
        covering = fibers_overlapping(nuc_i, nuc, p.chrom, p.start, p.end) | \
            fibers_overlapping(msp_i, msp, p.chrom, p.start, p.end)
        n_cover = len(covering)
        acc = len(fibers_overlapping(msp_i, msp, p.chrom, p.start, p.end))
        occ = len(fibers_overlapping(nuc_i, nuc, p.chrom, p.start, p.end))
        m6a_n = sum(
            1 for i in m6a_i.query(p.chrom, p.start, p.end)
            for blk in m6a[i].blocks
            if blk[0] < p.end and blk[1] > p.start
        )
        fire_n = len(fibers_overlapping(fire_i, fire, p.chrom, p.start, p.end)) if fire_i else 0
        rows.append([
            p.chrom, p.start, p.end, p.name,
            n_cover, acc, occ, fire_n,
            round(safe_div(acc, n_cover), 6),
            round(safe_div(occ, n_cover), 6),
            round(safe_div(m6a_n, n_cover), 6),
            round(safe_div(fire_n, n_cover), 6) if fire_i else "NA",
        ])

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"{args.sample}.region_quant.tsv"
    write_tsv(str(out),
              ["chrom", "start", "end", "name", "n_fibers", "n_accessible",
               "n_nucleosome", "n_fire", "percent_accessible", "nuc_occupancy",
               "m6a_per_fiber", "fire_frac"], rows)
    logger.info(f"Wrote {out} ({len(rows)} regions)")
    return 0


def merge_samples(args: argparse.Namespace) -> int:
    """Merge per-sample region_quant.tsv into regions x samples matrices."""
    tables: Dict[str, Dict[str, Dict[str, str]]] = {}
    for spec in args.quant.split(","):
        spec = spec.strip()
        if not spec:
            continue
        path = Path(spec)
        sample = path.name.split(".region_quant.tsv")[0]
        tables[sample] = {}
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                key = f"{f[0]}:{f[1]}-{f[2]}|{f[3]}"
                tables[sample][key] = dict(zip(header, f))
        logger.info(f"Loaded {len(tables[sample])} regions from {path} (sample={sample})")

    if not tables:
        logger.error("no quant tables provided")
        return 1
    # preserve first-table region order
    first = next(iter(tables))
    order = list(tables[first].keys())
    samples = list(tables.keys())

    metrics = ["percent_accessible", "nuc_occupancy", "m6a_per_fiber", "fire_frac"]
    args.outdir.mkdir(parents=True, exist_ok=True)
    for metric in metrics:
        rows = []
        for key in order:
            chrom_pos, name = key.split("|")
            chrom, span = chrom_pos.split(":")
            start, end = span.split("-")
            row = [chrom, start, end, name]
            for s in samples:
                row.append(tables[s].get(key, {}).get(metric, "NA"))
            rows.append(row)
        out = args.outdir / f"{metric}_matrix.tsv"
        write_tsv(str(out), ["chrom", "start", "end", "name"] + samples, rows)
        logger.info(f"Wrote {out}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="quant_matrix.py",
        description="Per-region Fiber-seq quantification and matrix merge.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("run", help="quantify one sample")
    pr.add_argument("-p", "--peaks", type=Path, required=True, help="regions/peaks BED")
    pr.add_argument("-s", "--msp", type=Path, required=True, help="MSP BED12 (.bed.gz)")
    pr.add_argument("-u", "--nuc", type=Path, required=True, help="nucleosome BED12 (.bed.gz)")
    pr.add_argument("-a", "--m6a", type=Path, required=True, help="m6A BED12 (.bed.gz)")
    pr.add_argument("-f", "--fire", type=Path, default=None,
                    help="per-fiber FIRE BED (from `ft fire --extract`)")
    pr.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    pr.add_argument("--sample", type=str, default="sample", help="sample name")
    pr.add_argument("-l", "--log", type=Path, default=None, help="log file")
    pr.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")

    pm = sub.add_parser("merge", help="merge per-sample quant tables")
    pm.add_argument("-q", "--quant", type=str, required=True,
                    help="comma-separated <sample>.region_quant.tsv paths")
    pm.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    pm.add_argument("-l", "--log", type=Path, default=None, help="log file")
    pm.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    if args.command == "run":
        for f in filter(None, [args.peaks, args.msp, args.nuc, args.m6a, args.fire]):
            if not Path(f).exists():
                logger.error(f"input not found: {f}")
                return 1
        if args.dry_run:
            logger.info(f"[DRY-RUN] would quantify {args.peaks} -> {args.outdir}")
            return 0
        return run_sample(args)
    if args.dry_run:
        logger.info(f"[DRY-RUN] would merge {args.quant} -> {args.outdir}")
        return 0
    return merge_samples(args)


if __name__ == "__main__":
    sys.exit(main())
