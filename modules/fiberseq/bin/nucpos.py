#!/usr/bin/env python3
"""Nucleosome positioning analysis around anchors (TSS or peak centers).

Collects nucleosome centers per fiber in a window around each anchor,
assigns each nucleosome to a relative slot (i-th nearest upstream/downstream
of the anchor), and computes per-slot mean offset (phasing) and dispersion
decay, following the analysis in Cell Rep Methods 2024 [1].

Usage:
    python nucpos.py -a anchors.bed --nuc sample.nuc.bed.gz \
        -o results/sample/nucpos --sample sample
    python nucpos.py --gtf genes.gtf --anchor tss --nuc sample.nuc.bed.gz -o out/
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
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

from src.common.util.FiberlibUtil import RegionIndex, read_bed, write_tsv  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("NucPos", log_file=str(log_file))
    else:
        lg = setup_logger("NucPos")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def tss_anchors_from_gtf(gtf: str) -> List[Tuple[str, int, str, str]]:
    """Return (chrom, tss, strand, gene_symbol) for every gene in the GTF."""
    out = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = {}
            for kv in f[8].split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                parts = kv.split(None, 1)
                if len(parts) == 2:
                    attrs[parts[0]] = parts[1].strip('"')
            strand = f[6]
            tss = int(f[3]) - 1 if strand != "-" else int(f[4])
            out.append((f[0], tss, strand, attrs.get("gene_name", attrs.get("gene_id", "."))))
    return out


def slot_offsets(
    nuc_index: RegionIndex,
    nuc_records,
    chrom: str,
    center: int,
    max_slots: int,
    max_reach: int,
) -> Tuple[Dict[int, List[int]], int]:
    """Per-slot nucleosome-center offsets for one anchor.

    Nucleosomes are taken per fiber and assigned to the i-th nearest slot
    upstream (negative) or downstream (positive) of the anchor. Returns
    slot -> list of (center - anchor) offsets and the fiber count.
    """
    per_fiber: Dict[str, List[int]] = defaultdict(list)
    for ri in nuc_index.query(chrom, center - max_reach, center + max_reach):
        r = nuc_records[ri]
        for bs, be in r.blocks:
            mid = (bs + be) // 2
            if center - max_reach <= mid <= center + max_reach:
                per_fiber[r.name].append(mid)
    slots: Dict[int, List[int]] = defaultdict(list)
    for centers in per_fiber.values():
        centers.sort()
        up = [c for c in centers if c < center]
        dn = [c for c in centers if c >= center]
        up.reverse()
        for i, c in enumerate(up[:max_slots], start=1):
            slots[-i].append(c - center)
        for i, c in enumerate(dn[:max_slots], start=1):
            slots[i].append(c - center)
    return slots, len(per_fiber)


def run(args: argparse.Namespace) -> int:
    if args.gtf is not None:
        anchors = tss_anchors_from_gtf(str(args.gtf))
        anchor_rows = [
            {"chrom": c, "center": t, "name": sym, "strand": st}
            for c, t, st, sym in anchors
        ]
    else:
        anchor_rows = []
        for r in read_bed(str(args.anchors)):
            anchor_rows.append({
                "chrom": r.chrom,
                "center": (r.start + r.end) // 2,
                "name": r.name,
                "strand": r.strand,
            })
    logger.info(f"Loaded {len(anchor_rows)} anchors")

    nuc = read_bed(str(args.nuc))
    nuc_index = RegionIndex(nuc)
    logger.info(f"Loaded {len(nuc)} fiber nucleosome rows")

    max_slots = max(args.upstream_nucs, args.downstream_nucs)
    max_reach = (max_slots + 1) * 300
    all_slots: Dict[int, List[int]] = defaultdict(list)
    anchor_counts = []
    rows = []
    for a in anchor_rows:
        slots, nfibers = slot_offsets(nuc_index, nuc, a["chrom"], a["center"],
                                      max_slots, max_reach)
        if nfibers < args.min_fibers:
            continue
        anchor_counts.append(nfibers)
        for slot, offsets in sorted(slots.items()):
            if slot < -args.upstream_nucs or slot > args.downstream_nucs:
                continue
            all_slots[slot].extend(offsets)
            rows.append([
                a["chrom"], a["center"], a["name"], a["strand"], slot,
                nfibers, len(offsets),
                round(float(np.mean(offsets)), 2),
                round(float(np.std(offsets)), 2),
            ])

    if not all_slots:
        logger.error(f"No anchor passed min_fibers={args.min_fibers}; no outputs")
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    per_anchor_tsv = args.outdir / f"{args.sample}.nuc_offset_per_anchor.tsv"
    write_tsv(str(per_anchor_tsv),
              ["chrom", "center", "anchor", "strand", "slot", "n_fibers",
               "n_obs", "mean_offset", "sd_offset"], rows)

    summary_rows = []
    for slot in sorted(all_slots):
        offs = np.array(all_slots[slot], dtype=float)
        summary_rows.append([
            slot, len(offs), round(float(np.mean(offs)), 2),
            round(float(np.median(offs)), 2), round(float(np.std(offs)), 2),
        ])
    summary_tsv = args.outdir / f"{args.sample}.nuc_offset_summary.tsv"
    write_tsv(str(summary_tsv),
              ["slot", "n_obs", "mean_offset", "median_offset", "sd_offset"],
              summary_rows)

    _plot(args, summary_rows)
    logger.info(f"Anchors used: {len(anchor_counts)}; wrote {per_anchor_tsv}, {summary_tsv}")
    return 0


def _plot(args, summary_rows) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    slots = [r[0] for r in summary_rows]
    meds = [r[3] for r in summary_rows]
    sds = [r[4] for r in summary_rows]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].boxplot(
        [[r[3]] for r in summary_rows], positions=slots, widths=0.5,
        medianprops={"color": "#c0392b"},
    )
    axes[0].plot(slots, meds, "o-", color="#2c3e50", ms=4)
    axes[0].axhline(0, ls="--", lw=0.8, color="grey")
    axes[0].set_xlabel("nucleosome slot relative to anchor (0 = nearest downstream)")
    axes[0].set_ylabel("median offset from anchor (bp)")
    axes[0].set_title(f"{args.sample} nucleosome positioning")

    axes[1].plot([abs(s) for s in slots], sds, "o-", color="#8e44ad", ms=4)
    axes[1].set_xlabel("|slot| distance to anchor")
    axes[1].set_ylabel("offset SD (bp)")
    axes[1].set_title("phasing decay")
    fig.tight_layout()
    fig.savefig(args.outdir / f"{args.sample}.nucpos.pdf")
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="nucpos.py",
        description="Nucleosome positioning/phasing analysis around anchors.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-a", "--anchors", type=Path,
                     help="anchor BED (center of each interval used as anchor)")
    src.add_argument("-g", "--gtf", type=Path,
                     help="GTF to derive TSS anchors (use with --anchor tss)")
    p.add_argument("--anchor", type=str, default="tss", choices=["tss", "bed"],
                   help="anchor source semantics")
    p.add_argument("-u", "--nuc", type=Path, required=True, help="nucleosome BED12 (.bed.gz)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-5", "--upstream-nucs", type=int, default=5, help="slots upstream")
    p.add_argument("-3", "--downstream-nucs", type=int, default=5, help="slots downstream")
    p.add_argument("-i", "--min-fibers", type=int, default=10,
                   help="minimum fibers per anchor to include")
    p.add_argument("-l", "--log", type=Path, default=None,
                   help="log file path (default: log to stdout)")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    if args.gtf is None and args.anchors is None:
        logger.error("provide --anchors or --gtf")
        return 1
    for f in filter(None, [args.gtf, args.anchors, args.nuc]):
        if not Path(f).exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would analyze nucleosome positioning -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
