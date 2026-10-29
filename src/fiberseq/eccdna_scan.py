#!/usr/bin/env python3
"""Detect circular-DNA-informative reads in a Fiber-seq fire BAM and profile them.

Concatemeric reads -- a circular template sequenced through multiple times --
appear in the BAM as one alignment chain with periodic large internal
deletions (period ~ circle length, high identity) or as supplementary
alignments tiling the same locus repeatedly. For each candidate:

  - candidate circle locus + estimated circle size (median period)
  - supporting reads, periods, mappability proxy (MAPQ)
  - TE composition of the circle interval (fraction bp in each TE class /
    family from a TE BED) -> "what sequences the circles are made of"
  - methylation on the supporting reads' intervals: m6A / mCpG call density
    per kb from the ft element BEDs (answers "is eccDNA methylated?")

Also writes a per-sample summary comparable across samples (DMSO vs PlaB).

Usage:
  python eccdna_scan.py --bam sample.fire.sorted.bam --fai genome.fa.fai \
      --te-bed te.bed -s SAMPLE=MSP,NUC,CPG,M6A -o ecc_out/ --sample DMSO
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pysam

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402
from src.common.util.FiberlibUtil import write_tsv, safe_div  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("EccDnaScan", log_file=str(log_file))
    else:
        lg = setup_logger("EccDnaScan")
    for h in lg.handlers:
        if type(h) is logging.StreamHandler:
            h.setStream(sys.stdout)
    return lg


logger = _get_logger()


def find_concatemer_candidates(bam_path: str, min_mapq: int, min_period: int,
                               max_period: int, min_identity: float,
                               min_periods: int, max_reads: int) -> List[dict]:
    """Reads whose primary alignment shows >= min_periods large internal dels.

    A periodic internal deletion chain (equal-sized D blocks within
    [min_period, max_period], similar sizes) indicates the read passed over
    the same circle multiple times.
    """
    candidates = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for n, read in enumerate(bam.fetch(until_eof=True)):
            if n >= max_reads:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapq < min_mapq or read.query_length < 2 * min_period:
                continue
            cig = read.cigartuples or []
            dels = [size for op, size in cig if op == 2 and size >= min_period]
            if len(dels) < min_periods:
                continue
            # periodicity: dels sizes similar AND the walking position advances
            # by ~circle length between consecutive dels
            pos = read.reference_start
            periods = []
            for op, size in cig:
                if op in (0, 2, 3, 7, 8):
                    if op == 2 and size >= min_period:
                        periods.append(size)
                    if op in (0, 7, 8):
                        pos += size
                elif op == 1:
                    pos += 0
            if len(periods) < min_periods:
                continue
            med = float(np.median(periods))
            if not (min_period <= med <= max_period):
                continue
            ok = sum(1 for p in periods if abs(p - med) <= 0.25 * med)
            if ok < min_periods:
                continue
            candidates.append({
                "read": read.query_name,
                "chrom": read.reference_name,
                "start": read.reference_start,
                "end": read.reference_end,
                "n_periods": len(periods),
                "period_median": med,
                "mapq": read.mapq,
            })
    return candidates


def cluster_candidates(cands: List[dict], window: int) -> List[dict]:
    """Group candidate reads by (chrom, period) proximity into circle loci."""
    by_chrom = defaultdict(list)
    for c in cands:
        by_chrom[c["chrom"]].append(c)
    loci = []
    for chrom, items in sorted(by_chrom.items()):
        items.sort(key=lambda c: c["start"])
        cluster = []
        for c in items:
            if cluster and c["start"] - cluster[-1]["start"] > window:
                loci.append(_summarize_cluster(cluster))
                cluster = []
            cluster.append(c)
        if cluster:
            loci.append(_summarize_cluster(cluster))
    return loci


def _summarize_cluster(cluster: List[dict]) -> dict:
    periods = [c["period_median"] for c in cluster]
    ends = [c["end"] for c in cluster if c["end"] is not None]
    return {
        "chrom": cluster[0]["chrom"],
        "start": min(c["start"] for c in cluster),
        "end": max(ends) if ends else max(c["start"] for c in cluster),
        "n_reads": len({c["read"] for c in cluster}),
        "circle_size_median": float(np.median(periods)),
        "period_cv": float(np.std(periods) / max(np.mean(periods), 1e-9)),
        "reads": cluster,
    }


def te_composition(loci: List[dict], te_bed: str) -> List[dict]:
    """Fraction of each circle interval covered by each TE class/family."""
    per_chrom = defaultdict(list)
    with open(te_bed) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                per_chrom[f[0]].append((int(f[1]), int(f[2]),
                                        f[3] if len(f) > 3 else ".",
                                        f[6] if len(f) > 6 else ".",
                                        f[7] if len(f) > 7 else "."))
            except ValueError:
                continue
    for c in per_chrom:
        per_chrom[c].sort()
    out = []
    for locus in loci:
        s, e = locus["start"], locus["end"]
        cov: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        total_te = []
        for ts, te, sub, fam, cls in per_chrom.get(locus["chrom"], []):
            if ts >= e:
                break
            if te <= s:
                continue
            s0, e0 = max(ts, s), min(te, e)
            cov[f"{cls}:{fam}"].append((s0, e0))
            total_te.append((s0, e0))
        span = e - s

        def union_bp(spans):
            spans = sorted(spans)
            bp, cs, ce = 0, None, None
            for a, b in spans:
                if cs is None:
                    cs, ce = a, b
                elif a <= ce:
                    ce = max(ce, b)
                else:
                    bp += ce - cs
                    cs, ce = a, b
            if cs is not None:
                bp += ce - cs
            return bp

        comp = {k: round(union_bp(v) / span, 4) for k, v in sorted(cov.items())
                if union_bp(v) / span >= 0.01}
        locus["te_fraction_total"] = round(union_bp(total_te) / span, 4)
        locus["te_composition"] = comp
        out.append(locus)
    return out


def methylation_on_loci(loci: List[dict], paths: Dict[str, str]) -> List[dict]:
    """m6A / mCpG call density per kb over each circle locus (all samples)."""
    import subprocess

    for key in ("m6a", "cpg"):
        for locus in loci:
            region = f"{locus['chrom']}:{locus['start']}-{locus['end']}"
            try:
                res = subprocess.run(
                    ["tabix", paths[key], region],
                    capture_output=True, text=True, timeout=120)
                lines = res.stdout.splitlines() if res.returncode == 0 else []
            except (OSError, subprocess.TimeoutExpired):
                lines = []
            n_calls = 0
            read_names = set()
            for ln in lines:
                f = ln.split("\t")
                if len(f) < 12:
                    continue
                sizes = [int(x) for x in f[10].split(",") if x]
                starts = [int(x) for x in f[11].split(",") if x]
                blocks = list(zip(starts, [s + b for s, b in zip(starts, sizes)]))
                if len(blocks) > 2:
                    n_calls += len(blocks) - 2  # drop sentinels
                read_names.add(f[3])
            locus[f"{key}_calls"] = n_calls
            locus[f"{key}_reads"] = len(read_names)
    return loci


def parse_sample(spec: str) -> Tuple[str, Dict[str, str]]:
    label, _, rest = spec.partition("=")
    parts = [p for p in rest.split(",") if p]
    d = dict(zip(["msp", "nuc", "cpg", "m6a"], parts))
    d["label"] = label
    return label, d


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="eccdna_scan.py",
        description="Concatemer-read eccDNA detection + TE composition + methylation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam", required=True, help="fire.sorted.bam (indexed)")
    p.add_argument("--te-bed", required=True)
    p.add_argument("-s", "--sample", action="append", required=True,
                   metavar="LABEL=MSP,NUC,CPG,M6A",
                   help="element BEDs used for methylation extraction")
    p.add_argument("--sample-name", required=True, help="sample label for outputs")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-period", type=int, default=200,
                   help="min circle size (bp) = min internal deletion size")
    p.add_argument("--max-period", type=int, default=100000)
    p.add_argument("--min-periods", type=int, default=3,
                   help="min repeats of the circle within one read")
    p.add_argument("--cluster-window", type=int, default=5000)
    p.add_argument("--max-reads", type=int, default=20_000_000)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("-l", "--log", type=Path, default=None)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    logger.info(f"scanning {args.bam} for concatemer reads ...")
    cands = find_concatemer_candidates(args.bam, args.min_mapq, args.min_period,
                                       args.max_period, 0, args.min_periods,
                                       args.max_reads)
    logger.info(f"candidate reads: {len(cands)}")
    loci = cluster_candidates(cands, args.cluster_window)
    logger.info(f"circle loci: {len(loci)}")
    loci = te_composition(loci, args.te_bed)

    samples = dict(parse_sample(s) for s in args.sample)
    for label, paths in samples.items():
        loci = methylation_on_loci(loci, paths)

    rows = []
    for locus in loci:
        rows.append({
            "chrom": locus["chrom"], "start": locus["start"], "end": locus["end"],
            "n_reads": locus["n_reads"],
            "circle_size_median": locus["circle_size_median"],
            "period_cv": round(locus["period_cv"], 4),
            "te_fraction_total": locus["te_fraction_total"],
            "te_composition": ";".join(f"{k}={v}" for k, v in locus["te_composition"].items()),
            **{f"{label}_{k}": locus.get(k, 0) for label in samples for k in ("m6a_calls", "m6a_reads", "cpg_calls", "cpg_reads")},
        })
    header = ["chrom", "start", "end", "n_reads", "circle_size_median", "period_cv",
              "te_fraction_total", "te_composition"] + \
             [f"{label}_{k}" for label in samples
              for k in ("m6a_calls", "m6a_reads", "cpg_calls", "cpg_reads")]
    write_tsv(str(args.outdir / f"{args.sample_name}.eccdna_loci.tsv"), header,
              [[r.get(c, "") for c in header] for r in rows])

    size_hist: Dict[str, int] = defaultdict(int)
    for locus in loci:
        bucket = int(min(locus["circle_size_median"] // 1000, 100))
        size_hist[f"{bucket}kb"] += 1
    write_tsv(str(args.outdir / f"{args.sample_name}.eccdna_size_hist.tsv"),
              ["size_bucket", "n_loci"], sorted(size_hist.items()))
    logger.info(f"done -> {args.outdir} ({len(loci)} loci)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
