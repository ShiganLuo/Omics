#!/usr/bin/env python3
"""Build telomere / centromere-satellite / mito region BEDs from a RepeatMasker .out.

Three region sets (all with a region-class name in column 4, consumed by
region_profile.py):

  telomere : merged terminal (CCCTAA)n/(TTAGGG)n Simple_repeat arrays within
             --terminal-dist of a chromosome end (per fai), plus 100 kb
             subtelomeric flanks. Classes: tel_array_5p/3p, subtel_flank_5p/3p.
  censat   : Satellite-class hits clustered into arrays (gap <= --sat-gap),
             classified by consensus-name patterns (--major-pattern GSAT,
             --minor-pattern MMSAT for mouse), plus 10 kb flanks.
             Classes: majorSat_array / minorSat_core / otherSat / flank_10kb.
  mito     : the mitochondrial contig (auto-guessed as the shortest fai
             sequence shorter than 30 kb, override with --mito-name).

Usage:
  python build_region_bed.py --rm-out asm.bp.p_ctg.fa.out --fai genome.fa.fai \
      -o regions/ --mode both
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

TELO_PATTERNS = ("acccta", "ttaggg", "ccctaa", "ctaggg")


def _get_logger(log_file: Optional[Path] = None):
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("BuildRegionBed", log_file=str(log_file))
    else:
        lg = setup_logger("BuildRegionBed")
    for h in lg.handlers:
        if type(h) is logging.StreamHandler:
            h.setStream(sys.stdout)
    return lg


logger = _get_logger()


def read_fai(path: Path) -> Dict[str, int]:
    """Parse .fai into {contig: length}."""
    lengths: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2:
                lengths[f[0]] = int(f[1])
    return lengths


def iter_rm_hits(rm_out: Path):
    """Yield (chrom, start0, end, name, cls) from a RepeatMasker .out."""
    in_table = False
    with open(rm_out) as fh:
        for line in fh:
            if not in_table:
                if line.lstrip().startswith("SW") and "perc" in line:
                    in_table = True
                continue
            if not line.strip() or line.startswith("----"):
                continue
            f = line.split()
            if len(f) < 11:
                continue
            try:
                beg, end = int(f[5]), int(f[6])
            except ValueError:
                continue
            yield f[4], beg - 1, end, f[9], f[10]


def merge_spans(spans: List[Tuple[int, int]], gap: int) -> List[Tuple[int, int]]:
    """Merge sorted (start, end) spans with gap tolerance."""
    merged: List[Tuple[int, int]] = []
    for s, e in sorted(spans):
        if merged and s - merged[-1][1] <= gap:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def build_telomere(hits, lengths: Dict[str, int], out: Path, terminal_dist: int,
                   tel_gap: int, subtel: int) -> None:
    """Terminal telomeric-repeat arrays + subtelomeric flanks."""
    per_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, s, e, name, cls in hits:
        if "simple_repeat" not in cls.lower():
            continue
        low = name.lower()
        if not any(p in low for p in TELO_PATTERNS):
            continue
        clen = lengths.get(chrom)
        if clen is None:
            continue
        if s <= terminal_dist or e >= clen - terminal_dist:
            per_chrom.setdefault(chrom, []).append((s, e))
    rows = []
    for chrom, spans in sorted(per_chrom.items()):
        clen = lengths[chrom]
        for s, e in merge_spans(spans, tel_gap):
            mid = (s + e) / 2
            side = "5p" if mid < clen / 2 else "3p"
            rows.append((chrom, s, e, f"tel_array_{side}"))
            fs, fe = (max(0, s - subtel), s) if side == "5p" else (e, min(clen, e + subtel))
            if fe > fs:
                rows.append((chrom, fs, fe, f"subtel_flank_{side}"))
    rows.sort()
    with open(out, "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\t.\t.\n")
    n_arr = sum(1 for r in rows if r[3].startswith("tel_array"))
    logger.info(f"telomere: {n_arr} terminal arrays (+flanks) -> {out}")


def build_censat(hits, lengths: Dict[str, int], out: Path, sat_gap: int,
                 flank: int, major_pat: str, minor_pat: str,
                 class_pat: str) -> None:
    """Clustered satellite arrays, classified, plus flanks."""
    per_key: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for chrom, s, e, name, cls in hits:
        if class_pat and class_pat.lower() not in cls.lower():
            continue
        low = name.lower()
        if minor_pat and minor_pat.lower() in low:
            kind = "minorSat_core"
        elif major_pat and major_pat.lower() in low:
            kind = "majorSat_array"
        else:
            kind = "otherSat"
        per_key.setdefault((chrom, kind), []).append((s, e))
    rows = []
    for (chrom, kind), spans in sorted(per_key.items()):
        clen = lengths.get(chrom)
        for s, e in merge_spans(spans, sat_gap):
            rows.append((chrom, s, e, kind))
            if clen is not None:
                for fs, fe in ((max(0, s - flank), s), (e, min(clen, e + flank))):
                    if fe > fs:
                        rows.append((chrom, fs, fe, "flank_10kb"))
    rows.sort()
    with open(out, "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\t.\t.\n")
    from collections import Counter
    cnt = Counter(r[3] for r in rows)
    logger.info(f"censat: {dict(cnt)} -> {out}")


def build_mito(lengths: Dict[str, int], out: Path, mito_name: Optional[str]) -> None:
    """Single whole-mitochondrion region (region_profile bins it)."""
    name = mito_name
    if name is None:
        cands = [(ln, c) for c, ln in lengths.items() if ln < 30000]
        if not cands:
            logger.warning("no mito contig found (shortest contig >= 30kb); skipping mito bed")
            return
        cands.sort()
        name = cands[0][1]
    if name not in lengths:
        logger.error(f"--mito-name {name} not in fai; skipping")
        return
    with open(out, "w") as fh:
        fh.write(f"{name}\t0\t{lengths[name]}\tMT_genome\t.\t.\n")
    logger.info(f"mito: {name} ({lengths[name]} bp) -> {out}")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="build_region_bed.py",
        description="Build telomere/censat/mito region BEDs from RepeatMasker .out + fai.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--rm-out", type=Path, required=True, help="RepeatMasker .out file")
    p.add_argument("--fai", type=Path, required=True, help="genome .fai")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--mode", choices=["telomere", "censat", "mito", "all"], default="all")
    p.add_argument("--prefix", default=None, help="output prefix (default: rm-out stem)")
    p.add_argument("--terminal-dist", type=int, default=1000000,
                   help="distance from chromosome end for terminal arrays")
    p.add_argument("--tel-gap", type=int, default=50000, help="merge gap for telomere arrays")
    p.add_argument("--subtel", type=int, default=100000, help="subtelomeric flank size")
    p.add_argument("--sat-gap", type=int, default=10000, help="merge gap for satellite arrays")
    p.add_argument("--flank", type=int, default=10000, help="censat flank size")
    p.add_argument("--major-pattern", default="GSAT", help="major satellite name pattern")
    p.add_argument("--minor-pattern", default="MMSAT", help="minor satellite name pattern")
    p.add_argument("--class-pattern", default="Satellite", help="keep hits whose class contains this")
    p.add_argument("--mito-name", default=None, help="override mito contig name")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file (default stdout)")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    lengths = read_fai(args.fai)
    logger.info(f"fai: {len(lengths)} contigs")
    prefix = args.prefix or args.rm_out.name.split(".")[0]
    hits = list(iter_rm_hits(args.rm_out))
    logger.info(f"RepeatMasker hits: {len(hits)}")
    if args.mode in ("telomere", "all"):
        build_telomere(hits, lengths, args.outdir / f"{prefix}.telomere.bed",
                       args.terminal_dist, args.tel_gap, args.subtel)
    if args.mode in ("censat", "all"):
        build_censat(hits, lengths, args.outdir / f"{prefix}.censat.bed",
                     args.sat_gap, args.flank, args.major_pattern,
                     args.minor_pattern, args.class_pattern)
    if args.mode in ("mito", "all"):
        build_mito(lengths, args.outdir / f"{prefix}.mito.bed", args.mito_name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
