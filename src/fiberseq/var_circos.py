#!/usr/bin/env python3
"""Genome-wide variant + signal circos plot from Fiber-seq WGS products.

Rings (outer -> inner):
  1. cytoband-style chromosome ideograms (from fai lengths)
  2. SNV/indel density (deepvariant VCF, 1 Mb bins)
  3. SV arcs (pbsv VCF: DEL/INS/DUP/INV/BND), colored by type
  4. FIRE peak density (1 Mb bins)
  5. m6A density (ft m6a bed.gz, 1 Mb bins; sampled for speed)
  6. mCpG density (ft cpg bed.gz, 1 Mb bins; sampled)

Self-contained matplotlib circos (no circos/perl dependency). VCFs are read
with gzip (plain text parse of bcftools-compatible VCF).

Usage:
  python var_circos.py --fai genome.fa.fai --vcf s.vcf.gz --sv-vcf s.sv.vcf.gz \
      --fire-peaks s.fire_peaks.bed --m6a s.m6a.bed.gz --cpg s.cpg.bed.gz \
      -o circos_out/ --sample DMSO
"""
from __future__ import annotations

import argparse
import gzip
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("VarCircos", log_file=str(log_file))
    else:
        lg = setup_logger("VarCircos")
    for h in lg.handlers:
        if type(h) is logging.StreamHandler:
            h.setStream(sys.stdout)
    return lg


logger = _get_logger()


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_fai(path: str) -> List[tuple]:
    """Primary chromosomes in fai order (skip unplaced/alt < 5 Mb)."""
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 2:
                continue
            ln = int(f[1])
            if ln < 5_000_000:
                continue
            out.append((f[0], ln))
    return out


def chrom_bins(chroms, bin_size: int):
    """Assign each chromosome a contiguous bin range on a linear axis."""
    layout, total = {}, 0
    for c, ln in chroms:
        nb = int(np.ceil(ln / bin_size)) + 1  # +1 spacer bin
        layout[c] = (total, nb, ln)
        total += nb
    return layout, total


def binned_density(path: str, layout: Dict, bin_size: int, weight_col=None,
                   chrom_col=0, pos_mode="start", sample_rate: int = 1) -> np.ndarray:
    """Density of BED/VCF-like records in 1-D bins across all chromosomes."""
    arr = np.zeros(sum(nb for _, nb, _ in layout.values()), dtype=np.float64)
    n_seen = 0
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            chrom = f[chrom_col]
            if chrom not in layout:
                continue
            n_seen += 1
            if sample_rate > 1 and (n_seen % sample_rate):
                continue
            try:
                pos = int(f[1]) if pos_mode == "start" else int(f[1])
            except (ValueError, IndexError):
                continue
            base, nb, ln = layout[chrom]
            b = min(int(pos // bin_size), nb - 2)
            w = 1.0
            if weight_col is not None and len(f) > weight_col:
                try:
                    w = float(f[weight_col])
                except ValueError:
                    w = 1.0
            arr[base + b] += w
    return arr


def vcf_density(path: str, layout: Dict, bin_size: int) -> np.ndarray:
    return binned_density(path, layout, bin_size, chrom_col=0, pos_mode="start")


def load_svs(path: str) -> List[tuple]:
    """SV records (chrom, pos, chrom2, pos2, svtype) from a pbsv VCF."""
    svs = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 8:
                continue
            info = f[7]
            svtype = "BND"
            for kv in info.split(";"):
                if kv.startswith("SVTYPE="):
                    svtype = kv.split("=", 1)[1]
                    break
            chrom, pos = f[0], int(f[1])
            chrom2, pos2 = chrom, pos
            if svtype == "BND" and "]" in f[4] or "[" in f[4]:
                alt = f[4]
                for tok in alt.replace("]", "[").split("["):
                    if ":" in tok and any(ch.isalpha() for ch in tok):
                        try:
                            c2, p2 = tok.rsplit(":", 1)
                            chrom2, pos2 = c2, int(p2)
                            break
                        except ValueError:
                            pass
            if svtype in ("DEL", "INS", "DUP", "INV"):
                for kv in info.split(";"):
                    if kv.startswith("END="):
                        try:
                            chrom2, pos2 = chrom, int(kv.split("=", 1)[1])
                        except ValueError:
                            pass
                        break
            svs.append((chrom, pos, chrom2, pos2, svtype))
    return svs


def plot_circos(chroms, layout, total_bins, densities, svs, sample: str,
                outdir: Path, formats: Sequence[str], bin_size: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Wedge, Arc

    fig = plt.figure(figsize=(11, 11))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_direction(-1)  # type: ignore[attr-defined]
    ax.set_theta_offset(np.pi / 2)  # type: ignore[attr-defined]
    ax.set_axis_off()

    T = 2 * np.pi
    gap = 0.015 * T
    usable = T - gap * len(chroms)
    n_real = sum(nb - 1 for _, nb, _ in layout.values())
    bp_bin = usable / n_real

    def bin_to_theta(bin_i: int) -> float:
        return bin_i * bp_bin

    def chrom_theta(chrom: str, pos: int):
        base, nb, ln = layout[chrom]
        b = base + min(int(pos // bin_size), nb - 2)
        return bin_to_theta(b)

    # ideogram
    r_ideo = (0.92, 0.98)
    for i, (c, ln) in enumerate(chroms):
        base, nb, _ = layout[c]
        t0 = bin_to_theta(base) + gap / 2
        t1 = bin_to_theta(base + nb - 1) - gap / 2
        ax.add_patch(Wedge((0, 0), r_ideo[1], np.degrees(t0), np.degrees(t1),
                           width=r_ideo[1] - r_ideo[0], facecolor="#888888",
                           edgecolor="white", lw=0.5, alpha=0.9))
        mid = (t0 + t1) / 2
        ax.text(mid, 1.02, c.replace("CM", "").replace("NC_", "").split(".")[0],
                ha="center", va="center", fontsize=7)

    # density rings
    ring_specs = [(0.86, densities.get("snv"), "SNV/indel density", "Blues"),
                  (0.76, densities.get("fire"), "FIRE peak density", "Oranges"),
                  (0.66, densities.get("m6a"), "m6A density", "Greens"),
                  (0.56, densities.get("mcpg"), "mCpG density", "Purples")]
    for r_out, arr, label, cmap_name in ring_specs:
        if arr is None:
            continue
        arr = arr / max(arr.max(), 1e-9)
        cmap = matplotlib.colormaps[cmap_name]
        for i, (c, ln) in enumerate(chroms):
            base, nb, _ = layout[c]
            for b in range(nb - 1):
                v = arr[base + b]
                if v <= 0:
                    continue
                t0, t1 = bin_to_theta(base + b), bin_to_theta(base + b + 1)
                ax.add_patch(Wedge((0, 0), r_out, np.degrees(t0), np.degrees(t1),
                                   width=0.09, facecolor=cmap(0.25 + 0.75 * v),
                                   edgecolor="none"))
        ax.text(0, r_out, label, ha="center", va="center", fontsize=7)

    # SV arcs
    sv_colors = {"DEL": "#1f77b4", "INS": "#d62728", "DUP": "#2ca02c",
                 "INV": "#9467bd", "BND": "#ff7f0e"}
    r_sv = 0.50
    seen_types = set()
    for chrom, pos, chrom2, pos2, svtype in svs:
        if chrom not in layout or chrom2 not in layout:
            continue
        t1c, t2c = chrom_theta(chrom, pos), chrom_theta(chrom2, pos2)
        if abs(t1c - t2c) < 0.01:
            continue
        col = sv_colors.get(svtype, "gray")
        # route the arc through the shorter side
        diff = (t2c - t1c) % T
        if diff > np.pi:
            t1c, t2c = t2c, t1c
            diff = T - diff
        mid = t1c + diff / 2
        ax.plot([t1c, mid, t2c], [r_ideo[0] - 0.01, r_sv, r_ideo[0] - 0.01],
                color=col, lw=0.6, alpha=0.5, solid_capstyle="round")
        seen_types.add(svtype)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=sv_colors[t], lw=2, label=t)
               for t in sorted(seen_types)]
    if handles:
        ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(1.05, 0.05),
                  fontsize=8, title="SV")

    ax.set_title(f"{sample}: variants, SVs and Fiber-seq signals", fontsize=12, pad=20)
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(outdir / f"{sample}.circos.{fmt}", dpi=300)
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="var_circos.py",
        description="Self-contained circos of SNV/SV/FIRE/methylation densities.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fai", required=True, help="genome .fai")
    p.add_argument("--vcf", default=None, help="small-variant VCF(.gz)")
    p.add_argument("--sv-vcf", default=None, help="SV VCF(.gz) for arcs")
    p.add_argument("--fire-peaks", default=None, help="ft fire_peaks.bed")
    p.add_argument("--m6a", default=None, help="ft m6a.bed.gz")
    p.add_argument("--cpg", default=None, help="ft cpg.bed.gz")
    p.add_argument("--sample", required=True, help="sample label")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("--bin-size", type=int, default=1_000_000)
    p.add_argument("--meth-sample-rate", type=int, default=20,
                   help="keep 1/N methylation rows for density (speed)")
    p.add_argument("-f", "--formats", nargs="+", default=["png"],
                   choices=["png", "pdf", "svg", "eps"])
    p.add_argument("-l", "--log", type=Path, default=None)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    chroms = load_fai(args.fai)
    logger.info(f"chromosomes: {len(chroms)}")
    layout, total = chrom_bins(chroms, args.bin_size)

    densities: Dict[str, np.ndarray] = {}
    if args.vcf:
        densities["snv"] = vcf_density(args.vcf, layout, args.bin_size)
        logger.info(f"snv bins with signal: {(densities['snv'] > 0).sum()}")
    if args.fire_peaks:
        densities["fire"] = binned_density(args.fire_peaks, layout, args.bin_size)
    if args.m6a:
        densities["m6a"] = binned_density(args.m6a, layout, args.bin_size,
                                          sample_rate=args.meth_sample_rate)
    if args.cpg:
        densities["mcpg"] = binned_density(args.cpg, layout, args.bin_size,
                                           sample_rate=args.meth_sample_rate)

    svs = load_svs(args.sv_vcf) if args.sv_vcf else []
    logger.info(f"SVs: {len(svs)}")

    plot_circos(chroms, layout, total, densities, svs, args.sample,
                args.outdir, args.formats, args.bin_size)
    write_density_table(densities, layout, chroms, args.bin_size, args.outdir)
    logger.info(f"done -> {args.outdir}")
    return 0


def write_density_table(densities, layout, chroms, bin_size, outdir: Path) -> None:
    """Dump the binned densities behind the rings for reuse."""
    if not densities:
        return
    with open(outdir / "circos_bin_density.tsv", "w") as fh:
        fh.write("chrom\tbin_start\t" + "\t".join(densities.keys()) + "\n")
        for c, ln in chroms:
            base, nb, _ = layout[c]
            for b in range(nb - 1):
                vals = [f"{d[base + b]:.6f}" for d in densities.values()]
                fh.write(f"{c}\t{b * bin_size}\t" + "\t".join(vals) + "\n")


if __name__ == "__main__":
    sys.exit(main())
