#!/usr/bin/env python3
"""Measure aligned sequencing depth at telomere regions and quantify alignment bias.

Diagnostic (NOT a length estimator): computes mean *aligned* read depth over
telomere regions (telomere_core / subtelomere) and a euchromatic control, then
reports the recovery fraction (telomere depth / control depth). Pure (TTAGGG)n
telomeres align very poorly, so region-based alignment depth drastically
undercounts them — this script quantifies that bias to demonstrate why
alignment/region-depth is NOT a viable telomere-length estimator (use
read_density_telomere.py / telogator2 instead).

Usage:
    python telomere_region_depth.py -b s1.bam -d telomere_regions.bed \
        -c NC_133406.1:120000000-120100000 -s s1 -o results/

    # bed: 3+ columns, optional 4th column = class label (telomere_core/subtelomere/…)
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import pysam

try:
    from .LogUtil import setup_logger
    from .plot_util import normalize_formats, save_fig
except ImportError:
    from LogUtil import setup_logger
    from plot_util import normalize_formats, save_fig

import matplotlib.pyplot as plt

logger = setup_logger("telomere_region_depth", stream=sys.stdout)


DEFAULT_MIN_RECOVERY = 0.5


@dataclass
class RegionDepth:
    """Aligned-depth summary for one region."""
    contig: str
    start: int
    stop: int
    region_class: str
    mean_depth: float
    n_reads: int
    length: int


def mean_aligned_depth(bam: pysam.AlignmentFile, contig: str, start: int,
                       stop: int) -> tuple[float, int]:
    """Exact mean aligned depth over [start, stop): sum(read overlap bp)/length.

    Counts only primary mapped reads (skips unmapped/secondary/supplementary).
    Memory-light: accumulates overlap bases, does not build per-base arrays.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Open indexed BAM.
    contig, start, stop : str, int, int
        Region in 0-based half-open coordinates.

    Returns
    -------
    (mean_depth, n_reads) : tuple[float, int]
    """
    overlap = 0
    n = 0
    for read in bam.fetch(contig, start, stop):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        rs = read.reference_start
        re_ = read.reference_end
        if re_ is None:
            continue
        ov = min(re_, stop) - max(rs, start)
        if ov > 0:
            overlap += ov
            n += 1
    return overlap / (stop - start), n


def parse_region(spec: str) -> tuple[str, int, int]:
    """Parse 'contig:start-end' (1-based inclusive) -> (contig, start, stop) 0-based."""
    contig, rng = spec.rsplit(":", 1)
    s, e = rng.replace(",", "").split("-")
    return contig, int(s) - 1, int(e)


def load_bed(path: Path) -> list[tuple[str, int, int, str]]:
    """Read a BED file into (contig, start, stop, class) tuples (col4 optional)."""
    rows: list[tuple[str, int, int, str]] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        contig, start, stop = f[0], int(f[1]), int(f[2])
        region_class = f[3] if len(f) > 3 else "region"
        rows.append((contig, start, stop, region_class))
    return rows


def plot_depth_diagnostic(depths: list, class_mean: dict, control_depth,
                          tel_class, sample_name: str, out_base: str,
                          formats: Optional[List[str]] = None) -> List[str]:
    """Figure 4: aligned depth by region class + recovery (diagnostic QC).

    Panel 1 bars mean aligned depth per region class (+ euchromatic control);
    panel 2 shows the telomere-vs-control recovery fraction against the viability
    threshold. Demonstrates pure-telomere depth collapse (method QC, not biology).
    """
    names = list(class_mean.keys())
    vals = [class_mean[k] for k in names]
    if control_depth is not None:
        names.append("control")
        vals.append(control_depth)
    def _c(n):
        if n in ("telomere_core", "tel_core", "core"):
            return "#C0504D"
        return "#3C6E9F" if n == "control" else "#4C9F70"

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    axes[0].bar(names, vals, color=[_c(n) for n in names])
    axes[0].set_ylabel("Mean aligned depth (x)")
    axes[0].set_title("Aligned depth by region class")
    axes[0].tick_params(axis="x", rotation=20)

    ax = axes[1]
    if control_depth and tel_class and tel_class in class_mean:
        recovery = class_mean[tel_class] / control_depth
        ax.bar(["recovery"], [recovery], color="#C0504D")
        ax.axhline(0.5, ls="--", color="gray", label="viability threshold 0.5")
        ax.set_ylabel("Telomere depth / control depth")
        ax.set_title(f"Alignment recovery at '{tel_class}'")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "no control / telomere class", ha="center", va="center")
        ax.set_axis_off()

    fig.suptitle(f"Alignment-bias diagnostic — {sample_name}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    written = save_fig(fig, out_base, formats)
    plt.close(fig)
    return written


def main() -> None:
    """Parse CLI, compute per-region depth, and write TSV + a verdict."""
    parser = argparse.ArgumentParser(
        description="Quantify alignment bias at telomere regions (diagnostic).")
    parser.add_argument("-b", "--bam", required=True,
                        help="Coordinate-sorted indexed BAM (reads aligned to reference)")
    parser.add_argument("-d", "--bed", required=True,
                        help="Telomere regions BED (col4 = class label)")
    parser.add_argument("-c", "--control", default="",
                        help="Euchromatic control region 'contig:start-end' for baseline depth")
    parser.add_argument("-s", "--sample_name", required=True,
                        help="Sample identifier")
    parser.add_argument("-o", "--output_dir", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("-m", "--min_recovery", type=float, default=DEFAULT_MIN_RECOVERY,
                        help="Recovery fraction below which region depth is deemed unusable "
                             f"(default: {DEFAULT_MIN_RECOVERY})")
    parser.add_argument("-l", "--log", default="",
                        help="Optional log file path (default: log to stdout only)")
    parser.add_argument("-F", "--format", action="append", dest="formats", metavar="FMT",
                        help="Output figure format(s), repeatable (default: png)")
    args = parser.parse_args()

    log = logger
    if args.log:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        log = setup_logger("telomere_region_depth", log_file=args.log,
                           stream=sys.stdout)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    regions = load_bed(Path(args.bed))
    bam = pysam.AlignmentFile(args.bam, "rb")

    depths: list[RegionDepth] = []
    for contig, start, stop, rclass in regions:
        md, n = mean_aligned_depth(bam, contig, start, stop)
        depths.append(RegionDepth(contig, start, stop, rclass, md, n, stop - start))

    control_depth: Optional[float] = None
    if args.control:
        c_contig, c_start, c_stop = parse_region(args.control)
        control_depth, _ = mean_aligned_depth(bam, c_contig, c_start, c_stop)
    bam.close()

    by_class: dict[str, list[RegionDepth]] = {}
    for d in depths:
        by_class.setdefault(d.region_class, []).append(d)
    class_mean: dict[str, float] = {}
    for rclass, ds in by_class.items():
        tot_len = sum(d.length for d in ds)
        class_mean[rclass] = sum(d.mean_depth * d.length for d in ds) / tot_len

    tsv = out / f"{args.sample_name}_telomere_region_depth.tsv"
    with tsv.open("w") as fh:
        fh.write("contig\tstart\tstop\tclass\tlength\tmean_depth\tn_reads\n")
        for d in depths:
            fh.write(f"{d.contig}\t{d.start}\t{d.stop}\t{d.region_class}\t{d.length}\t"
                     f"{d.mean_depth:.4f}\t{d.n_reads}\n")

    stats = out / f"{args.sample_name}_telomere_region_depth_stats.txt"
    lines = [f"Sample: {args.sample_name}",
             f"BAM: {args.bam}",
             f"Regions: {len(depths)}", ""]
    lines.append("Mean aligned depth by region class:")
    for rclass, mv in sorted(class_mean.items()):
        lines.append(f"  {rclass:16s} {mv:8.3f}x")
    if control_depth is not None:
        lines.append(f"  {'control':16s} {control_depth:8.3f}x   (euchromatic baseline)")
    lines.append("")

    tel_class = None
    for cand in ("telomere_core", "tel_core", "core"):
        if cand in class_mean:
            tel_class = cand
            break
    if tel_class is None and class_mean:
        tel_class = min(class_mean, key=lambda k: class_mean[k])

    verdict: list[str] = []
    if control_depth and tel_class:
        recovery = class_mean[tel_class] / control_depth
        lines.append(f"Alignment recovery at '{tel_class}': "
                     f"{class_mean[tel_class]:.3f}x / {control_depth:.3f}x = {recovery:.4f}")
        lines.append("")
        if recovery < args.min_recovery:
            verdict.append(
                f"判定:不可取(NOT viable)。纯端粒区('{tel_class}')比对深度仅为常染色质的 "
                f"{recovery*100:.2f}%,大量 (TTAGGG)n 读段未比对/多比对而丢失。")
            verdict.append(
                "→ 区域比对深度法低估端粒,且残留读段是带唯一锚点的偏倚子集,"
                "不能用于端粒长度(ATL)估计。")
            verdict.append(
                "→ 端粒长度请用参考无关方法:read_density_telomere.py(全读段 k-mer 密度)"
                "或 telogator2(逐臂);本脚本仅作比对偏差诊断。")
        else:
            verdict.append(f"判定:比对恢复率 {recovery:.2f} 达标,区域深度尚可参考。")
    else:
        verdict.append("判定:需 --control 提供常染色质基线方能给出恢复率与不可取结论。")

    lines.extend(["", "结论 / Verdict:"] + verdict)
    stats.write_text("\n".join(lines) + "\n")
    _pb = str(out / f"{args.sample_name}_telomere_region_depth_diagnostic")
    for _p in plot_depth_diagnostic(depths, class_mean, control_depth, tel_class,
                                    args.sample_name, _pb, args.formats):
        log.info(f"Wrote: {_p}")


    for ln in lines:
        log.info(ln)
    log.info(f"Wrote: {tsv}")
    log.info(f"Wrote: {stats}")


if __name__ == "__main__":
    main()
