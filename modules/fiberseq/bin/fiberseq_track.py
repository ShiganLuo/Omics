#!/usr/bin/env python3
"""Generate genome-browser tracks from Fiber-seq element BEDs.

Produces binned bedGraph tracks (accessibility = MSP fiber fraction,
m6A density = m6A marks per fiber per kb) and converts them to bigWig with
bedGraphToBigWig (skipped with a warning when the binary is unavailable).
Also writes a UCSC trackDb snippet listing the generated tracks.

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
from typing import Dict, List, Optional, Sequence

import numpy as np

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
_ROOT_DIR = _BIN_DIR.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

from src.common.util.FiberlibUtil import atomic_write, block_signal, read_bed  # noqa: E402


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


def binned_bedgraph(records, sizes: Dict[str, int], bin_size: int,
                    uniq_fibers: bool) -> List[tuple]:
    """Compute binned mean signal per chromosome as bedGraph rows."""
    rows = []
    for chrom, size in sorted(sizes.items()):
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


def run(args: argparse.Namespace) -> int:
    msp = read_bed(str(args.msp))
    m6a = read_bed(str(args.m6a))
    nuc = read_bed(str(args.nuc)) if args.nuc else []
    sizes = read_chrom_sizes(str(args.chrom_sizes))
    logger.info(f"msp={len(msp)} m6a={len(m6a)} nuc={len(nuc)} chroms={len(sizes)}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    tracks = []

    acc_rows = binned_bedgraph(msp, sizes, args.bin, uniq_fibers=True)
    acc_bg = args.outdir / f"{args.sample}.accessibility.bedGraph"
    write_bedgraph(acc_bg, acc_rows)
    acc_bw = args.outdir / f"{args.sample}.accessibility.bw"
    if to_bigwig(acc_bg, args.chrom_sizes, acc_bw):
        tracks.append((f"{args.sample}_accessibility", "FIRE-seq MSP accessibility", acc_bw))
    else:
        tracks.append((f"{args.sample}_accessibility", "FIRE-seq MSP accessibility", acc_bg))

    cov_rows = binned_bedgraph(list(msp) + list(nuc), sizes, args.bin, uniq_fibers=True)
    m6a_rows = binned_bedgraph(m6a, sizes, args.bin, uniq_fibers=False)
    cov_map = {(c, s): v for c, s, e, v in cov_rows}
    dens_rows = [
        (c, s, e, round(v / max(cov_map.get((c, s), 1), 1), 5))
        for c, s, e, v in m6a_rows
    ]
    dens_bg = args.outdir / f"{args.sample}.m6a_density.bedGraph"
    write_bedgraph(dens_bg, dens_rows)
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
