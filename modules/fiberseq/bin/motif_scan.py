#!/usr/bin/env python3
"""Scan motif occurrences (MEME motif format) in region sequences.

Parses a MEME-format motif file (e.g. JASPAR export), extracts region
sequences from a genome FASTA (via pysam FastaFile), and scans both strands
with a log-odds PWM. Hits with score >= --pwm-threshold * max_score are
reported as a motif site BED (used by motif_profile / footprint_profile).

Usage:
    python motif_scan.py --motif-file JASPAR.meme --regions sample.fire_peaks.bed \
        --fasta genome.fa -o results/sample/motif --sample sample
"""
from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
_ROOT_DIR = _BIN_DIR.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

from src.common.util.FiberlibUtil import read_bed, write_tsv  # noqa: E402

BASES = ("A", "C", "G", "T")
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("MotifScan", log_file=str(log_file))
    else:
        lg = setup_logger("MotifScan")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def parse_meme(path: str) -> Dict[str, np.ndarray]:
    """Parse MEME motif file into name -> position-probability matrix."""
    motifs: Dict[str, List[List[float]]] = {}
    current: Optional[str] = None
    in_matrix = False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("MOTIF"):
                parts = s.split()
                current = parts[1] if len(parts) > 1 else f"motif{len(motifs) + 1}"
                motifs[current] = []
                in_matrix = False
                continue
            if s.startswith("letter-probability"):
                in_matrix = True
                continue
            if s.startswith("END MATRIX"):
                in_matrix = False
                continue
            if in_matrix and current is not None and not s.startswith(("URL", "MOTIF")):
                parts = s.split()
                if len(parts) == 4:
                    try:
                        motifs[current].append([float(x) for x in parts])
                    except ValueError:
                        in_matrix = False
    return {name: np.array(pwm) for name, pwm in motifs.items() if pwm}


def pwm_log_odds(ppm: np.ndarray, pseudocount: float = 0.001) -> np.ndarray:
    """Convert a probability matrix to log-odds against uniform background."""
    p = np.clip(ppm, pseudocount, 1.0)
    return np.log2(p / 0.25)


def revcomp_scores(scores: np.ndarray) -> np.ndarray:
    """Reverse-complement score matrix (positions x 4 in ACGT order)."""
    return scores[::-1][:, [3, 2, 1, 0]]


def scan_sequence(seq: str, scores: np.ndarray, threshold: float) -> List[Tuple[int, float, str]]:
    """Return (offset, score, strand) hits above threshold * max_score."""
    L = scores.shape[0]
    max_score = float(np.sum(np.max(scores, axis=1)))
    min_score = threshold * max_score
    hits = []
    seq = seq.upper()
    for strand, mat in (("+", scores), ("-", revcomp_scores(scores))):
        for i in range(len(seq) - L + 1):
            window = seq[i: i + L]
            score = 0.0
            for j, base in enumerate(window):
                idx = BASES.index(base) if base in BASES else -1
                score += mat[j][idx] if idx >= 0 else -2.0
            if score >= min_score:
                hits.append((i, float(score), strand))
    return hits


def run(args: argparse.Namespace) -> int:
    try:
        import pysam
    except ImportError:
        logger.error("pysam is required for FASTA access")
        return 1

    motifs = parse_meme(str(args.motif_file))
    if not motifs:
        logger.error(f"no motifs parsed from {args.motif_file} (MEME format expected)")
        return 1
    logger.info(f"Loaded {len(motifs)} motifs: {list(motifs)}")

    regions = read_bed(str(args.regions))
    fasta = pysam.FastaFile(str(args.fasta))
    logger.info(f"Scanning {len(regions)} regions on {args.fasta}")

    rows = []
    for r in regions:
        seq = fasta.fetch(r.chrom, r.start, r.end)
        for name, ppm in motifs.items():
            scores = pwm_log_odds(ppm)
            for offset, score, strand in scan_sequence(seq, scores, args.pwm_threshold):
                site_start = r.start + offset
                rows.append([
                    r.chrom, site_start, site_start + scores.shape[0],
                    f"{name}@{r.name}", round(score, 3), strand,
                ])
    if not rows:
        logger.warning("no motif hit above threshold; consider lowering --pwm-threshold")
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"{args.sample}.motif_sites.bed"
    write_tsv(str(out), ["chrom", "start", "end", "name", "score", "strand"], rows)
    logger.info(f"Wrote {len(rows)} motif sites to {out}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="motif_scan.py",
        description="Scan MEME motifs in region sequences from a genome FASTA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-m", "--motif-file", type=Path, required=True, help="MEME motif file")
    p.add_argument("-r", "--regions", type=Path, required=True, help="regions BED")
    p.add_argument("-f", "--fasta", type=Path, required=True, help="genome FASTA (indexed)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-t", "--pwm-threshold", type=float, default=0.8,
                   help="hit threshold as fraction of max PWM score")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in (args.motif_file, args.regions, args.fasta):
        if not f.exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would scan {args.motif_file} in {args.regions} -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
