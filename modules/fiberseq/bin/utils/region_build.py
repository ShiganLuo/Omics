#!/usr/bin/env python3
"""Build telomere / centromere-satellite region BEDs for Fiber-seq region analysis.

Subcommand `telomere`: convert the contig-end scan TSV from
`modules/telomere/bin/scan_assembly_telomere.py` into a region BED with
classes `telomere_core` (the telomeric repeat tract) and `subtelomere`
(a fixed window proximal to the tract).

Subcommand `sat`: convert a RepeatMasker .out into clustered satellite
array BED (one array per collinear hit cluster, max gap `--gap`), with the
region class taken from a consensus-name -> class map TSV (2 columns:
name pattern, class; first substring match wins; fallback `--other-class`).
The clustering mirrors modules/censatbuild/bin/censatbuild_rm2bed.py but
supports >2 satellite classes.

Usage:
    python region_build.py telomere --scan-tsv Mmul_T2T_assembly_telomere.tsv \
        --subtelomere 100000 -o Mmul_T2T.telomere.bed
    python region_build.py sat --rm-out genome.rm.out \
        --class-map sat_class_map.tsv --gap 10000 -o Mmul_T2T.sat_arrays.bed
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

_BIN_DIR = Path(__file__).resolve().parent
_ROOT_DIR = _BIN_DIR.parent.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("RegionBuild", log_file=str(log_file))
    else:
        lg = setup_logger("RegionBuild")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def write_bed(path: Path, rows: Sequence[Tuple[str, int, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3]))
    with open(path, "w") as out:
        for chrom, start, end, name in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    logger.info(f"Wrote {path} ({len(rows)} regions)")


def build_telomere(args: argparse.Namespace) -> int:
    rows: List[Tuple[str, int, int, str]] = []
    with open(args.scan_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom = f[idx["contig"]]
            clen = int(f[idx["contig_length"]])
            tel_start = int(f[idx["tel_start"]])
            tel_end = int(f[idx["tel_end"]])
            if tel_end - tel_start < args.min_len:
                logger.warning(f"{chrom} telomere tract {tel_end - tel_start}bp "
                               f"< --min-len {args.min_len}; skipped")
                continue
            rows.append((chrom, tel_start, tel_end, "telomere_core"))
            if args.subtelomere > 0:
                if f[idx["end"]] == "5prime":
                    sub_start, sub_end = tel_end, min(clen, tel_end + args.subtelomere)
                else:
                    sub_start, sub_end = max(0, tel_start - args.subtelomere), tel_start
                if sub_end > sub_start:
                    rows.append((chrom, sub_start, sub_end, "subtelomere"))
    write_bed(args.out, rows)
    return 0


def load_class_map(path: Path) -> List[Tuple[str, str]]:
    """First-match (pattern, class) list; pattern matches the consensus name."""
    mapping: List[Tuple[str, str]] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            mapping.append((f[0].strip().lower(), f[1].strip()))
    return mapping


def classify(name: str, mapping: Sequence[Tuple[str, str]],
             other_class: str) -> str:
    low = name.lower()
    for pat, cls in mapping:
        if pat in low:
            return cls
    return other_class


def build_cen(args: argparse.Namespace) -> int:
    """Pick the dominant array per chromosome as cen_core + flank windows."""
    best: Dict[str, Tuple[int, int]] = {}
    with open(args.arrays) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if f[3] != args.array_class:
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if chrom not in best or end - start > best[chrom][1] - best[chrom][0]:
                best[chrom] = (start, end)
    rows: List[Tuple[str, int, int, str]] = []
    for chrom, (start, end) in sorted(best.items()):
        if end - start < args.min_len:
            logger.warning(f"{chrom} dominant {args.array_class} array "
                           f"{end - start}bp < --min-len {args.min_len}; skipped")
            continue
        rows.append((chrom, start, end, "cen_core"))
        if args.flank > 0:
            if start - args.flank >= 0:
                rows.append((chrom, start - args.flank, start, "cen_flank"))
            rows.append((chrom, end, end + args.flank, "cen_flank"))
    logger.info(f"cen_core arrays: {sum(1 for r in rows if r[3] == 'cen_core')}")
    write_bed(args.out, rows)
    return 0


def build_sat(args: argparse.Namespace) -> int:
    mapping = load_class_map(args.class_map)
    hits: List[Tuple[str, int, int, str]] = []
    with open(args.rm_out) as fh:
        in_table = False
        for line in fh:
            if line.lstrip().startswith("SW") and "perc" in line:
                in_table = True
                continue
            if not in_table or line.startswith("----") or not line.strip():
                continue
            f = line.split()
            if len(f) < 11:
                continue
            try:
                start, end = int(f[5]), int(f[6])
            except ValueError:
                continue
            if start > end:  # minus-strand hits are still query coords
                start, end = end, start
            hits.append((f[4], start - 1, end,
                         classify(f[9], mapping, args.other_class)))
    logger.info(f"satellite hits kept: {len(hits)}")

    # cluster per (chrom, class), gap <= args.gap (as censatbuild_rm2bed)
    by_key: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    for chrom, start, end, cls in hits:
        by_key[(chrom, cls)].append((start, end))

    rows: List[Tuple[str, int, int, str]] = []
    for (chrom, cls), spans in sorted(by_key.items()):
        spans.sort()
        cs, ce = spans[0]
        for s, e in spans[1:]:
            if s - ce <= args.gap:
                ce = max(ce, e)
            else:
                rows.append((chrom, cs, ce, cls))
                cs, ce = s, e
        rows.append((chrom, cs, ce, cls))

    if args.min_len > 1:
        kept = [r for r in rows if r[2] - r[1] >= args.min_len]
        logger.info(f"arrays: {len(rows)} total, {len(kept)} >= {args.min_len}bp")
        rows = kept
    write_bed(args.out, rows)
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="region_build.py",
        description="Build telomere / satellite-array region BEDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)

    pt = sub.add_parser("telomere", help="telomere tracts + subtelomere windows")
    pt.add_argument("--scan-tsv", type=Path, required=True,
                    help="scan_assembly_telomere.py detail TSV")
    pt.add_argument("--subtelomere", type=int, default=100000,
                    help="subtelomere window proximal to the tract (bp, 0=off)")
    pt.add_argument("--min-len", type=int, default=1000,
                    help="min telomere tract length to keep (bp)")
    pt.add_argument("-o", "--out", type=Path, required=True)
    pt.add_argument("-l", "--log", type=Path, default=None)

    ps = sub.add_parser("sat", help="cluster RepeatMasker satellite hits")
    ps.add_argument("--rm-out", type=Path, required=True,
                    help="RepeatMasker .out file")
    ps.add_argument("--class-map", type=Path, required=True,
                    help="TSV: name pattern -> region class (first match wins)")
    ps.add_argument("--other-class", default="otherSat",
                    help="class for names not matching the map")
    ps.add_argument("--gap", type=int, default=10000,
                    help="max gap between hits within one array (bp)")
    ps.add_argument("--min-len", type=int, default=1,
                    help="min array length to keep (bp)")
    ps.add_argument("-o", "--out", type=Path, required=True)
    ps.add_argument("-l", "--log", type=Path, default=None)

    pc = sub.add_parser("cen", help="dominant satellite array per chromosome")
    pc.add_argument("--arrays", type=Path, required=True,
                    help="clustered arrays BED from `sat`")
    pc.add_argument("--array-class", default="alphaSat",
                    help="array class considered centromeric")
    pc.add_argument("--flank", type=int, default=100000,
                    help="pericentromere flank window each side (bp, 0=off)")
    pc.add_argument("--min-len", type=int, default=100000,
                    help="min dominant array length to accept (bp)")
    pc.add_argument("-o", "--out", type=Path, required=True)
    pc.add_argument("-l", "--log", type=Path, default=None)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    if args.command == "telomere":
        return build_telomere(args)
    if args.command == "cen":
        return build_cen(args)
    return build_sat(args)


if __name__ == "__main__":
    sys.exit(main())
