#!/usr/bin/env python3
"""Derive a hybrid-capture target BED from aligned BAM coverage.

Vendor target BEDs are often unavailable or built on an older genome build.
This tool reconstructs the effective capture region from the actual aligned
reads instead: per-base depths from ``samtools depth`` are thresholded, called
per-sample, intersected across samples (majority mode optional), and merged
with a minimum gap so that tiny interspersed gaps inside a captured exon do
not fragment the region.

Usage:
    derive_targets.py -b sample1.bam sample2.bam -o targets.bed \
        [--min-depth 10] [--min-bq 20] [--min-mq 20] [--min-gap 50] \
        [--fraction 1.0] [--min-region 50]

Output: sorted BED (chrom, start, end) of merged captured regions.
"""
import argparse
import subprocess
import sys
from pathlib import Path


def stream_depth(bam: str, min_depth: int, min_bq: int, min_mq: int):
    """Yield (chrom, pos0, covered) per base from samtools depth (streamed)."""
    cmd = [
        "samtools", "depth", "-a",
        "-q", str(min_bq),
        "-Q", str(min_mq),
        "-G", "UNMAP,SECONDARY,QCFAIL,DUP",
        bam,
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        chrom, pos, depth = parts[0], int(parts[1]) - 1, int(parts[2])
        yield chrom, pos, depth >= min_depth
    if proc.wait() != 0:
        raise RuntimeError(f"samtools depth failed for {bam} (exit {proc.returncode})")


def call_sample_regions(bam, min_depth, min_bq, min_mq):
    """Call covered runs (chrom, start, end) for one BAM."""
    regions = []
    cur_chrom, cur_start, cur_end, cur_covered = None, None, None, False
    for chrom, pos, covered in stream_depth(bam, min_depth, min_bq, min_mq):
        if chrom != cur_chrom:
            if cur_covered and cur_chrom is not None:
                regions.append((cur_chrom, cur_start, cur_end))
            cur_chrom, cur_start, cur_end, cur_covered = chrom, pos, pos + 1, covered
            continue
        if covered and not cur_covered:
            cur_start = pos
        elif not covered and cur_covered:
            regions.append((cur_chrom, cur_start, cur_end))
        cur_end = pos + 1
        cur_covered = covered
    if cur_covered and cur_chrom is not None:
        regions.append((cur_chrom, cur_start, cur_end))
    return regions


def merge_regions(regions, min_gap):
    """Merge regions within min_gap of each other (per chromosome)."""
    merged = []
    for chrom, start, end in sorted(regions):
        if merged and merged[-1][0] == chrom and start - merged[-1][1] <= min_gap:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def intersect_or_majority(per_sample, fraction):
    """Combine per-sample region lists into consensus covered intervals."""
    from collections import defaultdict
    idx = defaultdict(list)  # (chrom, bucket) -> per-base vote lists are heavy;
    # instead count regions overlapping via interval sweep on union boundaries
    events = defaultdict(list)  # chrom -> list of (pos, delta)
    for regions in per_sample:
        for chrom, start, end in regions:
            events[chrom].append((start, 1))
            events[chrom].append((end, -1))
    need = max(1, round(fraction * len(per_sample)))
    consensus = []
    for chrom, evs in events.items():
        evs.sort()
        depth = 0
        prev_pos = None
        for pos, delta in evs:
            if prev_pos is not None and depth >= need and pos > prev_pos:
                consensus.append((chrom, prev_pos, pos))
            depth += delta
            prev_pos = pos
    # merge touching intervals from the sweep
    out = []
    for chrom, start, end in sorted(consensus):
        if out and out[-1][0] == chrom and start <= out[-1][2]:
            out[-1] = (chrom, out[-1][1], max(out[-1][2], end))
        else:
            out.append((chrom, start, end))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-b", "--bams", nargs="+", required=True, help="aligned BAM files")
    ap.add_argument("-o", "--out", required=True, help="output BED path")
    ap.add_argument("--min-depth", type=int, default=10,
                    help="per-base depth threshold to call a base captured [10]")
    ap.add_argument("--min-bq", type=int, default=20, help="min base quality [20]")
    ap.add_argument("--min-mq", type=int, default=20, help="min mapping quality [20]")
    ap.add_argument("--min-gap", type=int, default=50,
                    help="merge gaps up to this many bp inside captured regions [50]")
    ap.add_argument("--fraction", type=float, default=1.0,
                    help="fraction of samples that must cover a region [1.0 = intersection]")
    ap.add_argument("--min-region", type=int, default=50,
                    help="drop final regions shorter than this many bp [50]")
    args = ap.parse_args()

    per_sample = []
    for bam in args.bams:
        regions = call_sample_regions(bam, args.min_depth, args.min_bq, args.min_mq)
        if not regions:
            print(f"WARNING: no regions passed the depth threshold in {bam}",
                  file=sys.stderr)
        per_sample.append(regions)

    consensus = intersect_or_majority(per_sample, args.fraction)
    merged = merge_regions(consensus, args.min_gap)
    final = [(c, s, e) for c, s, e in merged if e - s >= args.min_region]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w") as fh:
        for chrom, start, end in final:
            fh.write(f"{chrom}\t{start}\t{end}\n")
    tmp.replace(out)

    total_bp = sum(e - s for _, s, e in final)
    print(f"Derived {len(final)} target regions ({total_bp} bp) from "
          f"{len(args.bams)} BAM(s) -> {out}", file=sys.stderr)
    if not final:
        print("ERROR: empty target BED; lower --min-depth or check BAM alignment",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
