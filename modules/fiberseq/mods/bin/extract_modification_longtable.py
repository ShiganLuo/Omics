#!/usr/bin/env python3
"""Aggregate per-fiber modification calls (m6A / CpG) over FIRE peaks into
the long table sigtest consumes (feature x group x value).

For every FIRE peak region and every sample, counts m6A/CpG call sites
falling inside the region and normalizes by the number of fibers covering
it (calls per fiber, same denominator convention as quant_matrix's
m6a_per_fiber). Output is one row per region x sample:

    feature        <chrom>:<start>-<end>|<peak name>   (sigtest feature id)
    sample         sample id
    group          the sample's group from the design TSV
    value          modification calls per covering fiber (0 when no fiber
                   covers the region — the region stays testable)

Modification BEDs are BED12 (calls as blockStarts/blockSizes on the fiber
read); only call blocks overlapping the peak interval are counted. Files
are streamed with ChromSweep-style position sweep (multi-GB inputs), so
memory stays proportional to the peak neighborhood.

Usage:
    python extract_modification_longtable.py \
        --mod m6a:sampleA.m6a.bed.gz,sampleB.m6a.bed.gz \
        --peaks sampleA.fire_peaks.bed ... (one per sample, same order) \
        --design design.tsv --metric m6a_per_fiber \
        --out mods/sigtest/m6a_per_fiber_long.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO, Tuple

# Self-contained module bin scripts import the module-private helper
# (verbatim src/common copies) — modules.md 脚本撰写规范 #0.
from _fiberlib import BedRecord, ChromSweep, read_bed, setup_logger

logger = setup_logger("ModsLongTable", stream=sys.stdout)


def parse_group_samples(spec: str) -> Dict[str, str]:
    """Parse 'gA:s1+s2,gB:s3' (node.py group_samples) into sample -> group."""
    sample_group: Dict[str, str] = {}
    for tok in spec.split(","):
        if ":" not in tok:
            raise ValueError(f"malformed group spec token {tok!r} "
                             "(expected group:sample[+sample])")
        group, s_str = tok.split(":", 1)
        for s in s_str.split("+"):
            if s.strip():
                sample_group[s.strip()] = group.strip()
    return sample_group


def count_region_calls(
    peaks: Sequence[BedRecord], mod_bed: str
) -> List[Optional[int]]:
    """Call sites per peak region in one modification BED (None = no fiber).

    Streams the BED12 against the peak anchors; a call block counts when
    it overlaps the peak interval (blockStart < end and blockEnd > start).
    """
    sweep = ChromSweep([mod_bed], [(p.chrom, p.start, p.end) for p in peaks])
    counts: List[Optional[int]] = [None] * len(peaks)
    for k, (_chrom, _s, _e), ov in sweep.advance():
        n = 0
        p = peaks[k]
        # blocks are absolute (chromStart + blockStart, + blockSizes)
        for rec in ov[0]:
            for b_start, b_end in rec.iter_blocks():
                if b_start < p.end and b_end > p.start:
                    n += 1
        counts[k] = n
    return counts


def count_region_fibers(
    peaks: Sequence[BedRecord], msp_bed: str
) -> List[int]:
    """Fibers covering each peak region (MSP = accessibility footprint).

    The per-fiber denominator uses MSP records overlapping the region —
    the same ``n_fibers`` definition quant_matrix uses, keeping the two
    quantifications comparable.
    """
    sweep = ChromSweep([msp_bed], [(p.chrom, p.start, p.end) for p in peaks])
    fibers: List[int] = [0] * len(peaks)
    for k, (_chrom, _s, _e), ov in sweep.advance():
        fibers[k] = len({r.name for r in ov[0]})
    return fibers


def build_long_table(args: argparse.Namespace) -> int:
    """Aggregate every sample's modification BED over its peak set."""
    # sample-parallel inputs: mod beds, peak beds, msp beds share order
    entries = [e for e in args.mods.split(",") if e.strip()]
    sample_ids = [e.split(":")[0] for e in entries]
    mod_beds = [e.split(":", 1)[1] for e in entries]
    peak_beds = [p for p in args.peaks.split(",") if p.strip()]
    msp_beds = [m for m in args.msp.split(",") if m.strip()]
    if not (len(peak_beds) == len(mod_beds) == len(msp_beds)):
        logger.error(
            f"--mods/--peaks/--msp must list one path per sample in the "
            f"same order (mods={len(mod_beds)} peaks={len(peak_beds)} "
            f"msp={len(msp_beds)})"
        )
        return 1
    sample_group = parse_group_samples(args.group_samples)
    missing = [s for s in sample_ids if s not in sample_group]
    if missing:
        logger.error(f"samples absent from design TSV: {missing}")
        return 1
    logger.info(
        f"Aggregating {args.metric} over peaks for {len(sample_ids)} "
        f"samples: {sample_ids}"
    )

    # Region axis = union of every sample's peak names, first-seen order;
    # a sample missing a region (its own peaks differ) contributes no row
    # for that region — sigtest handles unequal n per group natively.
    region_order: List[str] = []
    region_peaks: Dict[str, BedRecord] = {}
    per_sample: Dict[str, Dict[str, float]] = {s: {} for s in sample_ids}
    for sid, mod_bed, peak_bed, msp_bed in zip(
        sample_ids, mod_beds, peak_beds, msp_beds
    ):
        peaks = read_bed(peak_bed)
        if not peaks:
            logger.error(
                f"sample {sid}: 0 peaks read from {peak_bed} — refusing to "
                f"build a long table with a missing sample (an empty peaks "
                f"file is usually call_peaks' zero-read placeholder or a "
                f"clobbered ft_call_peaks output; check its rule log)"
            )
            return 1
        n_calls = count_region_calls(peaks, mod_bed)
        n_fibers = count_region_fibers(peaks, msp_bed)
        for peak, calls, fibers in zip(peaks, n_calls, n_fibers):
            key = f"{peak.chrom}:{peak.start}-{peak.end}|{peak.name}"
            if key not in region_peaks:
                region_peaks[key] = peak
                region_order.append(key)
            # 0 when uncovered: region stays in the test instead of
            # silently dropping (a real biological zero)
            per_sample[sid][key] = calls / fibers if calls is not None and fibers else 0.0
        logger.info(f"sample {sid}: {len(peaks)} peaks aggregated")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("feature\tsample\tgroup\tvalue\n")
        for key in region_order:
            for sid in sample_ids:
                if key in per_sample[sid]:
                    fh.write(
                        f"{key}\t{sid}\t{sample_group[sid]}\t"
                        f"{per_sample[sid][key]:.6f}\n"
                    )
    logger.info(f"Wrote long table: {out} ({len(region_order)} regions)")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="extract_modification_longtable.py",
        description="Aggregate m6A/CpG calls over FIRE peaks into a sigtest long table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--mods", required=True,
                   help="comma list of sample:mod_bed.bed.gz (m6a or cpg)")
    sample_ids_help = "value is calls per covering fiber"
    p.add_argument("--peaks", required=True,
                   help="comma list of per-sample fire_peaks.bed (same order)")
    p.add_argument("--msp", required=True,
                   help="comma list of per-sample msp.bed.gz (fiber denominator)")
    p.add_argument("--group-samples", required=True,
                   help="group spec 'gA:s1+s2,gB:s3' (sample -> group mapping)")
    p.add_argument("--metric", default="m6a_per_fiber",
                   help=f"metric name in outputs; {sample_ids_help}")
    p.add_argument("--out", required=True, help="output long TSV path")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    if args.log:
        logger = setup_logger("ModsLongTable", log_file=str(args.log),
                              stream=sys.stdout)
    logger.info(f"Command: {' '.join(sys.argv)}")
    if args.dry_run:
        logger.info(f"[DRY-RUN] would aggregate {args.mods} -> {args.out}")
        return 0
    return build_long_table(args)


if __name__ == "__main__":
    sys.exit(main())
