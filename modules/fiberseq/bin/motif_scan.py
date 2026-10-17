#!/usr/bin/env python3
"""Scan motif occurrences (MEME motif format) in region sequences.

Parses a MEME-format motif file (e.g. JASPAR export), extracts region
sequences from a genome FASTA (via pysam FastaFile), and scans both strands
with a log-odds PWM. Hits with score >= --pwm-threshold * max_score are
reported as a motif site BED (used by motif_profile / footprint_profile).

Implementation note (vectorized):
    Regions are concatenated into chunked pseudo-sequences separated by
    (Lmax-1) 'N' bases, so no sliding window can span two regions (any
    window touching an 'N' scores at most 4 - 2*L, far below
    threshold * max_score). Window scores are accumulated per motif
    position with numpy slice additions in j = 0..L-1 order, which
    reproduces the legacy scalar loop bit-exactly (same float64 additions,
    same order), so output rows are identical to the legacy implementation
    -- same hits, same scores, same order. Progress is logged per chunk;
    --threads parallelizes over chunks.

Usage:
    python motif_scan.py --motif-file JASPAR.meme --regions sample.fire_peaks.bed \
        --fasta genome.fa -o results/sample/motif --sample sample [--threads 4]
"""
from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import sys
import time
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

# Encode A/C/G/T -> 0..3; anything else (N and other IUPAC codes) -> 4,
# which is scored as -2.0 exactly like the legacy scalar loop.
_BASE_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(BASES):
    _BASE_LUT[ord(_b)] = _i

_CHUNK_BP = 2_000_000  # pseudo-sequence base budget per worker chunk

# Globals inherited by forked workers (avoids re-pickling per task).
_FASTA_PATH: Optional[str] = None
_MOTIFS: List[Tuple[str, np.ndarray]] = []
_THRESHOLD: float = 0.8


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
    """Legacy scalar scanner (kept as reference / for equivalence tests).

    Return (offset, score, strand) hits above threshold * max_score.
    """
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


def _score_windows(idx: np.ndarray, scores: np.ndarray,
                   threshold: float) -> Tuple[Optional[np.ndarray], float]:
    """Sliding-window log-odds scores over an encoded pseudo-sequence.

    Accumulates motif position j contributions in j ascending order into the
    window score vector -- the exact float64 operations of the legacy scalar
    loop, so scores are bit-identical. Returns (window_scores, min_score);
    window_scores is None when no full window fits.
    """
    L = scores.shape[0]
    n = idx.shape[0]
    W = n - L + 1
    max_score = float(np.sum(np.max(scores, axis=1)))
    min_score = threshold * max_score
    if W <= 0:
        return None, min_score
    scores5 = np.concatenate([scores, np.full((L, 1), -2.0, dtype=np.float64)], axis=1)
    ws = np.zeros(W, dtype=np.float64)
    for j in range(L):
        ws += scores5[j][idx[j: j + W]]
    return ws, min_score


def _scan_chunk(region_chunk: List[Tuple[str, int, int, str]]) -> List[list]:
    """Scan one chunk of regions; return rows in legacy output order."""
    import pysam

    fasta = pysam.FastaFile(_FASTA_PATH)
    lmax = max(sc.shape[0] for _, sc in _MOTIFS)
    sep = max(0, lmax - 1)

    parts: List[bytes] = []
    seg_starts: List[int] = []
    seg_lens: List[int] = []
    pos = 0
    for i, (chrom, start, end, _name) in enumerate(region_chunk):
        seq = fasta.fetch(chrom, start, end).upper()
        if i:
            parts.append(b"N" * sep)
            pos += sep
        seg_starts.append(pos)
        seg_lens.append(len(seq))
        parts.append(seq.encode("ascii"))
        pos += len(seq)
    fasta.close()

    idx = _BASE_LUT[np.frombuffer(b"".join(parts), dtype=np.uint8)]
    seg_starts_a = np.asarray(seg_starts, dtype=np.int64)
    seg_ends_a = seg_starts_a + np.asarray(seg_lens, dtype=np.int64)

    buckets: List[List[list]] = [[] for _ in region_chunk]
    for name, ppm in _MOTIFS:  # MEME file order (matches legacy row order)
        scores = pwm_log_odds(ppm)
        L = scores.shape[0]
        for strand, sc in (("+", scores), ("-", revcomp_scores(scores))):
            ws, min_score = _score_windows(idx, sc, _THRESHOLD)
            if ws is None:
                continue
            mask = ws >= min_score
            if not mask.any():
                continue
            hit_pos = np.flatnonzero(mask)
            ridx = np.searchsorted(seg_ends_a, hit_pos, side="right")
            offsets = hit_pos - seg_starts_a[ridx]
            svals = ws[hit_pos]
            for k in range(hit_pos.shape[0]):
                r = int(ridx[k])
                chrom, start, _end, rname = region_chunk[r]
                site_start = start + int(offsets[k])
                buckets[r].append([
                    chrom, site_start, site_start + L,
                    f"{name}@{rname}", round(float(svals[k]), 3), strand,
                ])
    rows: List[list] = []
    for b in buckets:
        rows.extend(b)
    return rows


def run(args: argparse.Namespace) -> int:
    global _FASTA_PATH, _MOTIFS, _THRESHOLD
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
    fasta = pysam.FastaFile(str(args.fasta))  # validates FASTA + index early
    fasta.close()
    _FASTA_PATH = str(args.fasta)
    _MOTIFS = list(motifs.items())
    _THRESHOLD = float(args.pwm_threshold)

    total_bp = sum(r.end - r.start for r in regions)
    chunks: List[List[Tuple[str, int, int, str]]] = []
    cur: List[Tuple[str, int, int, str]] = []
    cur_bp = 0
    for r in regions:
        cur.append((r.chrom, r.start, r.end, r.name))
        cur_bp += r.end - r.start
        if cur_bp >= _CHUNK_BP:
            chunks.append(cur)
            cur, cur_bp = [], 0
    if cur:
        chunks.append(cur)

    n_workers = max(1, int(args.threads))
    logger.info(
        f"Scanning {len(regions)} regions ({total_bp / 1e6:.1f} Mbp) on {args.fasta} "
        f"in {len(chunks)} chunk(s), {n_workers} worker(s)"
    )

    rows: List[list] = []
    t0 = time.time()
    done_bp = 0
    if n_workers == 1 or len(chunks) <= 1:
        results = map(_scan_chunk, chunks)
        pool = None
    else:
        ctx = mp.get_context("fork")
        pool = ctx.Pool(processes=min(n_workers, len(chunks)))
        results = pool.imap(_scan_chunk, chunks)
    try:
        for ci, chunk_rows in enumerate(results):
            rows.extend(chunk_rows)
            chunk_bp = sum(e - s for _c, s, e, _n in chunks[ci])
            done_bp += chunk_bp
            dt = time.time() - t0
            pct = done_bp / total_bp * 100 if total_bp else 100.0
            logger.info(
                f"chunk {ci + 1}/{len(chunks)} done: {len(chunks[ci])} regions, "
                f"{chunk_bp / 1e6:.1f} Mbp, {dt:.1f}s elapsed ({pct:.0f}% of bases)"
            )
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    if not rows:
        logger.warning("no motif hit above threshold; consider lowering --pwm-threshold")
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"{args.sample}.motif_sites.bed"
    write_tsv(str(out), ["chrom", "start", "end", "name", "score", "strand"], rows)
    logger.info(f"Wrote {len(rows)} motif sites to {out} "
                f"(total {time.time() - t0:.1f}s)")
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
    p.add_argument("--threads", type=int, default=4,
                   help="worker processes for chunk-parallel scanning")
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
