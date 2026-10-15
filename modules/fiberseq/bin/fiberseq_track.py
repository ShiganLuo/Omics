#!/usr/bin/env python3
"""Generate genome-browser tracks from Fiber-seq element BEDs.

Produces binned bedGraph tracks (accessibility = MSP fiber fraction,
m6A density = m6A marks per fiber per kb) and converts them to bigWig with
bedGraphToBigWig (skipped with a warning when the binary is unavailable).
Also writes a UCSC trackDb snippet listing the generated tracks.

Input BEDs must be grouped by chromosome (sort -k1,1). Records are streamed
one chromosome at a time, so peak memory is bounded by the largest
chromosome instead of the whole file.

Usage:
    python fiberseq_track.py --msp sample.msp.bed.gz --m6a sample.m6a.bed.gz \
        --nuc sample.nuc.bed.gz --chrom-sizes genome.chrom.sizes \
        -o results/sample/track --sample sample
"""
from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
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

from src.common.util.FiberlibUtil import (  # noqa: E402
    atomic_write,
    block_signal,
    iter_bed_chrom,
)


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("FiberseqTrack", log_file=str(log_file))
    else:
        lg = setup_logger("FiberseqTrack")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def read_chrom_sizes(path: str) -> Dict[str, int]:
    """Read .chrom.sizes or .fai (first two columns)."""
    sizes: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    return sizes


def binned_rows_for_chrom(records, chrom: str, size: int, bin_size: int,
                          uniq_fibers: bool) -> List[Tuple[str, int, int, float]]:
    """Binned mean signal rows for one chromosome's records."""
    rows: List[Tuple[str, int, int, float]] = []
    n_bins = size // bin_size + 1
    signal = block_signal(records, chrom, 0, n_bins * bin_size,
                          uniq_fibers=uniq_fibers)
    for b in range(n_bins):
        s0, e0 = b * bin_size, min((b + 1) * bin_size, size)
        if e0 <= s0:
            continue
        val = float(np.mean(signal[s0:e0]))
        if val > 0:
            rows.append((chrom, s0, e0, round(val, 5)))
    return rows


def write_bedgraph(path: Path, rows) -> None:
    def _w(fh):
        for chrom, s, e, v in rows:
            fh.write(f"{chrom}\t{s}\t{e}\t{v}\n")
    atomic_write(str(path), _w)


def to_bigwig(bedgraph: Path, chrom_sizes: Path, bigwig: Path) -> bool:
    """Convert bedGraph to bigWig when bedGraphToBigWig is on PATH."""
    exe = shutil.which("bedGraphToBigWig")
    if exe is None:
        logger.warning("bedGraphToBigWig not found on PATH; keeping bedGraph only")
        return False
    cmd = [exe, str(bedgraph), str(chrom_sizes), str(bigwig)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        logger.error(f"bedGraphToBigWig failed: {res.stderr.strip()}")
        return False
    return True


class _Stream:
    """One pass over chromosome-grouped BEDs, yielding per-chrom record lists.

    All sources are advanced in lockstep, so at any moment only one
    chromosome's records (per source) are resident. Files must share the
    same chromosome order (sort -k1,1 on each).
    """

    def __init__(self, paths: Sequence[str]) -> None:
        self._iters = [iter_bed_chrom(p) for p in paths]
        self._states: List[Optional[Tuple[str, list]]] = [
            next(it, None) for it in self._iters
        ]

    def advance(self, chrom: str) -> List[list]:
        """Return per-source records for *chrom*; exhausted sources give []."""
        out: List[list] = []
        for i, s in enumerate(self._states):
            if s is not None and s[0] == chrom:
                out.append(s[1])
                self._states[i] = next(self._iters[i], None)
            else:
                out.append([])
        return out

    def remaining_chroms(self) -> List[str]:
        return [s[0] for s in self._states if s is not None]


def run(args: argparse.Namespace) -> int:
    sizes = read_chrom_sizes(str(args.chrom_sizes))
    paths = [str(args.msp), str(args.m6a)]
    if args.nuc:
        paths.append(str(args.nuc))
    n_total = 0

    args.outdir.mkdir(parents=True, exist_ok=True)
    acc_bg = args.outdir / f"{args.sample}.accessibility.bedGraph"
    dens_bg = args.outdir / f"{args.sample}.m6a_density.bedGraph"
    acc_fh = open(acc_bg, "w")
    dens_fh = open(dens_bg, "w")
    try:
        cov_map: Dict[Tuple[str, int], float] = {}
        stream = _Stream(paths)
        for chrom in sizes:
            msp_c, m6a_c, nuc_c = stream.advance(chrom)[:3]
            n_total += len(msp_c) + len(m6a_c) + len(nuc_c)
            acc_rows = binned_rows_for_chrom(msp_c, chrom, sizes[chrom],
                                             args.bin, uniq_fibers=True)
            cov_rows = binned_rows_for_chrom(
                msp_c + nuc_c, chrom, sizes[chrom], args.bin, uniq_fibers=True)
            m6a_rows = binned_rows_for_chrom(m6a_c, chrom, sizes[chrom],
                                             args.bin, uniq_fibers=False)
            cov_map.update({(c, s): v for c, s, e, v in cov_rows})
            for c, s, e, v in acc_rows:
                acc_fh.write(f"{c}\t{s}\t{e}\t{v}\n")
            for c, s, e, v in m6a_rows:
                dens = round(v / max(cov_map.get((c, s), 1), 1), 5)
                dens_fh.write(f"{c}\t{s}\t{e}\t{dens}\n")
            del msp_c, m6a_c, nuc_c, acc_rows, cov_rows, m6a_rows
        leftover = stream.remaining_chroms()
        if leftover:
            raise ValueError(
                f"input BEDs are not grouped/sorted by chromosome "
                f"(sort -k1,1 each); chromosome(s) {sorted(set(leftover))} "
                f"appear outside the .fai order")
    finally:
        acc_fh.close()
        dens_fh.close()
    logger.info(f"records processed={n_total} chroms={len(sizes)}")

    tracks = []
    acc_bw = args.outdir / f"{args.sample}.accessibility.bw"
    if to_bigwig(acc_bg, args.chrom_sizes, acc_bw):
        tracks.append((f"{args.sample}_accessibility", "FIRE-seq MSP accessibility", acc_bw))
    else:
        tracks.append((f"{args.sample}_accessibility", "FIRE-seq MSP accessibility", acc_bg))

    dens_bw = args.outdir / f"{args.sample}.m6a_density.bw"
    if to_bigwig(dens_bg, args.chrom_sizes, dens_bw):
        tracks.append((f"{args.sample}_m6a", "FIRE-seq m6A density", dens_bw))
    else:
        tracks.append((f"{args.sample}_m6a", "FIRE-seq m6A density", dens_bg))

    def _w(fh):
        for name, desc, path in tracks:
            fh.write(f"track type=bigWig name=\"{name}\" description=\"{desc}\"\n")
            fh.write(f"bigDataUrl {path}\n\n")
    atomic_write(str(args.outdir / "trackDb.txt"), _w)
    logger.info(f"Wrote {len(tracks)} tracks and trackDb.txt to {args.outdir}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="fiberseq_track.py",
        description="Generate bedGraph/bigWig browser tracks from Fiber-seq BEDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-s", "--msp", type=Path, required=True, help="MSP BED12 (.bed.gz)")
    p.add_argument("-a", "--m6a", type=Path, required=True, help="m6A BED12 (.bed.gz)")
    p.add_argument("-u", "--nuc", type=Path, default=None, help="nucleosome BED12 (.bed.gz)")
    p.add_argument("-c", "--chrom-sizes", type=Path, required=True,
                   help="chrom sizes file (.chrom.sizes or .fai)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-b", "--bin", type=int, default=10, help="bin size (bp)")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in filter(None, [args.msp, args.m6a, args.nuc, args.chrom_sizes]):
        if not Path(f).exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would build tracks -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
