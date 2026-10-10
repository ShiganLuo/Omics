#!/usr/bin/env python3
"""Annotate FIRE peaks with gene context and repeat overlap.

Classifies each peak as promoter / exon / intron / intergenic against a GTF
annotation and flags repeat overlap against a repeat BED (e.g. RepeatMasker).
Region class priority: promoter > exon > intron > intergenic.

Usage:
    python fire_annotate.py -p sample.fire_peaks.bed -g gencode.gtf \
        -r repeat.bed -o results/sample/annotate --sample sample
    python fire_annotate.py -p peaks.bed -g genes.gtf -o out/ -n  # dry-run
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
_ROOT_DIR = _BIN_DIR.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

from src.common.util.FiberlibUtil import BedRecord, RegionIndex, read_bed, safe_div, write_tsv  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("FireAnnotate", log_file=str(log_file))
    else:
        lg = setup_logger("FireAnnotate")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def parse_gtf(gtf: str) -> Tuple[Dict[str, dict], List[BedRecord], List[BedRecord]]:
    """Parse GTF into per-gene info, gene-body records, and exon records.

    Returns
    -------
    genes : dict
        gene_id -> {"symbol", "chrom", "strand", "tss", "start", "end"}
    bodies : list of BedRecord
        Gene-body intervals (name = gene_id) for intron classification.
    exons : list of BedRecord
        Exon intervals (name = gene_id) for exon classification.
    """
    genes: Dict[str, dict] = {}
    exon_spans: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            attrs = {}
            for kv in f[8].split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                parts = kv.split(None, 1)
                if len(parts) == 2:
                    attrs[parts[0]] = parts[1].strip('"')
            gid = attrs.get("gene_id")
            if not gid:
                continue
            chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            if f[2] == "gene":
                genes[gid] = {
                    "symbol": attrs.get("gene_name", gid),
                    "chrom": chrom,
                    "strand": strand,
                    "tss": start if strand != "-" else end - 1,
                    "start": start,
                    "end": end,
                }
            else:
                exon_spans[gid].append((start, end))
    bodies = [
        BedRecord(g["chrom"], g["start"], g["end"], gid)
        for gid, g in genes.items()
    ]
    exons = [
        BedRecord(genes[gid]["chrom"], s, e, gid)
        for gid, spans in exon_spans.items()
        if gid in genes
        for (s, e) in spans
    ]
    return genes, bodies, exons


def annotate(
    peaks: Sequence[BedRecord],
    genes: Dict[str, dict],
    bodies: Sequence[BedRecord],
    exons: Sequence[BedRecord],
    repeat_bed: Optional[str],
    tss_flank: int,
    repeat_frac: float,
) -> List[List]:
    """Build annotation rows for all peaks."""
    tss_records = [
        BedRecord(g["chrom"], max(0, g["tss"] - tss_flank), g["tss"] + tss_flank + 1, gid)
        for gid, g in genes.items()
    ]
    tss_index = RegionIndex(tss_records)
    tss_name = {i: r.name for i, r in enumerate(tss_records)}
    exon_index = RegionIndex(exons)
    gene_index = RegionIndex(bodies)

    repeats: List[BedRecord] = read_bed(repeat_bed) if repeat_bed else []
    rep_index = RegionIndex(repeats)

    rows = []
    for p in peaks:
        mid = (p.start + p.end) // 2
        cls, gene_sym = "intergenic", ""
        best = (1 << 60, "")
        for gi in tss_index.query(p.chrom, mid - tss_flank, mid + tss_flank + 1):
            gid = tss_name[gi]
            g = genes[gid]
            d = abs(g["tss"] - mid)
            if d < best[0]:
                best = (d, g["symbol"])
        if best[1]:
            cls, gene_sym = "promoter", best[1]
        else:
            for gi in exon_index.query(p.chrom, p.start, p.end):
                cls, gene_sym = "exon", genes[exons[gi].name]["symbol"]
                break
            else:
                for gi in gene_index.query(p.chrom, p.start, p.end):
                    cls, gene_sym = "intron", genes[bodies[gi].name]["symbol"]
                    break
        if not gene_sym:
            # intergenic (or exon/intron without symbol): nearest gene by TSS
            nearest = (1 << 60, "")
            for gid, g in genes.items():
                if g["chrom"] != p.chrom:
                    continue
                d = abs(g["tss"] - mid)
                if d < nearest[0]:
                    nearest = (d, g["symbol"])
            gene_sym = nearest[1]

        rep_bases = 0
        rep_names = set()
        plen = max(1, p.end - p.start)
        for ri in rep_index.query(p.chrom, p.start, p.end):
            r = repeats[ri]
            ov = min(p.end, r.end) - max(p.start, r.start)
            if ov > 0:
                rep_bases += ov
                rep_names.add(r.name if r.name != "." else "repeat")
        rep_frac = safe_div(rep_bases, plen)
        rows.append([
            p.chrom, p.start, p.end, p.name, p.score,
            cls, gene_sym,
            "yes" if rep_frac >= repeat_frac else "no",
            round(rep_frac, 4),
            ",".join(sorted(rep_names)) if rep_names else ".",
        ])
    return rows


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="fire_annotate.py",
        description="Annotate FIRE peaks with gene context and repeat overlap.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-p", "--peaks", type=Path, required=True, help="FIRE peaks BED")
    p.add_argument("-g", "--gtf", type=Path, required=True, help="gene annotation GTF/GFF")
    p.add_argument("-r", "--repeat-bed", type=Path, default=None,
                   help="repeat annotation BED (e.g. RepeatMasker)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("-s", "--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-t", "--tss-flank", type=int, default=1000,
                   help="promoter window around TSS (bp)")
    p.add_argument("-f", "--repeat-frac", type=float, default=0.25,
                   help="min fraction of peak overlapping repeat to flag it")
    p.add_argument("-l", "--log", type=Path, default=None,
                   help="log file path (default: log to stdout)")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")

    if not args.gtf.exists():
        logger.error(f"GTF not found: {args.gtf} — set genome.references.<asm>.gtf in config JSON")
        return 1
    if args.repeat_bed is not None and not args.repeat_bed.exists():
        logger.error(f"repeat BED not found: {args.repeat_bed}")
        return 1

    if args.dry_run:
        logger.info(f"[DRY-RUN] would annotate {args.peaks} with {args.gtf} "
                    f"(repeat={args.repeat_bed}) -> {args.outdir}")
        return 0

    logger.info("Parsing GTF")
    genes, bodies, exons = parse_gtf(str(args.gtf))
    logger.info(f"Loaded {len(genes)} genes, {len(exons)} exon intervals")

    peaks = read_bed(str(args.peaks))
    logger.info(f"Loaded {len(peaks)} peaks from {args.peaks}")
    rows = annotate(peaks, genes, bodies, exons,
                    str(args.repeat_bed) if args.repeat_bed else None,
                    args.tss_flank, args.repeat_frac)

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_table = args.outdir / f"{args.sample}.peaks_annotated.tsv"
    write_tsv(str(out_table),
              ["chrom", "start", "end", "name", "score", "region_class",
               "nearest_gene", "repeat_flag", "repeat_frac", "repeat_names"],
              rows)

    class_counts: Dict[str, int] = defaultdict(int)
    rep_counts = {"yes": 0, "no": 0}
    for r in rows:
        class_counts[r[5]] += 1
        rep_counts[r[7]] += 1
    summary_rows = (
        [["region_class", k, v] for k, v in sorted(class_counts.items())]
        + [["repeat_flag", k, v] for k, v in sorted(rep_counts.items())]
        + [["total", "peaks", len(rows)]]
    )
    out_summary = args.outdir / f"{args.sample}.peaks_summary.tsv"
    write_tsv(str(out_summary), ["category", "key", "count"], summary_rows)
    logger.info(f"Wrote {out_table} and {out_summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
