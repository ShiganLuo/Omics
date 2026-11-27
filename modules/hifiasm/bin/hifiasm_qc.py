#!/usr/bin/env python3
"""Assembly QC for a hifiasm contig FASTA: statistics + telomere-end capture.

One streaming pass computes per-contig composition (length / GC / N) and scans
the first/last window of every contig for terminal telomeric repeats
(CCCTAA/TTAGGG). Writes per-contig TSVs, summaries, and three figures:
assembly overview (length distribution / Nx / GC / rank-size), telomere QC
(status counts / length distribution / landscape / completeness), and a
telomere capture map over the top contigs.

Usage:
    python hifiasm_qc.py --fasta asm.bp.p_ctg.fa --sample_name s1 \
        --output_dir out/s1/assembly/qc --expected-chrom-ends 42
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import numpy as np

try:
    from .LogUtil import setup_logger
    from .plot_util import normalize_formats, save_fig
except ImportError:
    from LogUtil import setup_logger
    from plot_util import normalize_formats, save_fig

import matplotlib.pyplot as plt

logger = setup_logger("hifiasm_qc", stream=sys.stdout)

TEL_MOTIF_FWD = "CCCTAA"
TEL_MOTIF_REV = "TTAGGG"
DEFAULT_SCAN_LENGTH = 100_000
DEFAULT_WINDOW_SIZE = 200
DEFAULT_STEP_SIZE = 50
DEFAULT_MIN_TEL_FRAC = 0.5
DEFAULT_MIN_TEL_LEN = 1_000
DEFAULT_END_MARGIN = 1_000
DEFAULT_TOP_CONTIGS = 30
MAJOR_CONTIG_MIN = 1_000_000


@dataclass
class ContigRecord:
    """Per-contig composition plus both terminal windows."""
    name: str
    length: int
    gc_bases: int
    n_bases: int
    head: str
    tail: str


@dataclass
class EndCall:
    """Telomere call for one contig end."""
    tel_length: int
    tel_start: int
    tel_end: int
    end_gap: int
    motif: str
    supporting_windows: int
    total_windows: int
    end_gap_margin: int = DEFAULT_END_MARGIN
    min_tel_len: int = DEFAULT_MIN_TEL_LEN

    @property
    def captured(self) -> bool:
        """True when a terminal telomere passes length and proximity cutoffs."""
        return (self.tel_length >= self.min_tel_len
                and self.end_gap <= self.end_gap_margin)


@dataclass
class ContigQC:
    """QC record for one contig."""
    record: ContigRecord
    end5: EndCall
    end3: EndCall

    @property
    def status(self) -> str:
        """Terminal capture class: both / 5prime_only / 3prime_only / none."""
        has5 = self.end5.captured
        has3 = self.end3.captured
        if has5 and has3:
            return "both"
        if has5:
            return "5prime_only"
        if has3:
            return "3prime_only"
        return "none"


def iter_contig_records(fasta_path: str,
                        scan_length: int) -> Iterator[ContigRecord]:
    """Stream a FASTA, yielding per-contig composition plus end windows.

    Sequences are never retained whole: only running character counts and the
    first/last ``scan_length`` bases of the current contig are kept.
    """
    name: Optional[str] = None
    length = 0
    gc = 0
    ns = 0
    head_parts: List[str] = []
    head_len = 0
    tail = ""
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield ContigRecord(name, length, gc, ns,
                                       "".join(head_parts), tail)
                name = line[1:].split()[0]
                length, gc, ns, head_parts, head_len, tail = 0, 0, 0, [], 0, ""
            else:
                s = line.strip().upper()
                length += len(s)
                gc += s.count("G") + s.count("C")
                ns += s.count("N")
                if head_len < scan_length:
                    take = s[:scan_length - head_len]
                    head_parts.append(take)
                    head_len += len(take)
                tail = (tail + s)[-scan_length:]
    if name is not None:
        yield ContigRecord(name, length, gc, ns, "".join(head_parts), tail)


def find_terminal_telomere(seq: str, from_start: bool, window_size: int,
                           step_size: int, min_tel_frac: float,
                           min_tel_len: int,
                           end_margin: int) -> EndCall:
    """Call the terminal telomeric region on one contig end window.

    Motif hits are clustered with gap tolerance; among clusters within
    ``end_gap`` <= end_margin of the physical end the largest is chosen
    (true terminal telomere). If no cluster is end-proximate the largest
    cluster is still reported (likely internal telomeric sequence) and the
    capture check fails via the recorded end_gap.

    Positions are relative to the window; end_gap is the distance from the
    cluster to the physical contig end.
    """
    seg_len = len(seq)
    hits = [(m.start(), m.end())
            for m in re.finditer(r"(CCCTAA|TTAGGG)", seq)]
    if not hits:
        return EndCall(0, 0, 0, seg_len, "", 0, 0,
                       end_gap_margin=end_margin, min_tel_len=min_tel_len)

    clusters: List[List[Tuple[int, int]]] = [[hits[0]]]
    for start, end in hits[1:]:
        if start - clusters[-1][-1][1] < 30:
            clusters[-1].append((start, end))
        else:
            clusters.append([(start, end)])

    def cluster_gap(cluster: List[Tuple[int, int]]) -> int:
        if from_start:
            return cluster[0][0]
        return seg_len - cluster[-1][1]

    def cluster_span(cluster: List[Tuple[int, int]]) -> int:
        return cluster[-1][1] - cluster[0][0]

    proximate = [c for c in clusters if cluster_gap(c) <= end_margin]
    best = max(proximate or clusters, key=cluster_span)
    cluster_start = best[0][0]
    cluster_end = best[-1][1]
    region = seq[cluster_start:cluster_end]
    motif = (TEL_MOTIF_FWD if region.count(TEL_MOTIF_FWD) >= region.count(TEL_MOTIF_REV)
             else TEL_MOTIF_REV)

    supporting = 0
    total_windows = 0
    for start in range(max(cluster_start, 0), cluster_end - window_size + 1,
                       step_size):
        window = seq[start:start + window_size]
        tel_frac = (window.count(TEL_MOTIF_FWD)
                    + window.count(TEL_MOTIF_REV)) * 6 / window_size
        total_windows += 1
        if tel_frac >= min_tel_frac:
            supporting += 1

    return EndCall(
        tel_length=cluster_span(best),
        tel_start=cluster_start,
        tel_end=cluster_end,
        end_gap=cluster_gap(best),
        motif=motif,
        supporting_windows=supporting,
        total_windows=total_windows,
        end_gap_margin=end_margin,
        min_tel_len=min_tel_len,
    )


def scan_contig(record: ContigRecord, window_size: int, step_size: int,
                min_tel_frac: float, min_tel_len: int,
                end_margin: int) -> ContigQC:
    """Run the terminal telomere call on both ends of one contig."""
    end5 = find_terminal_telomere(record.head, True, window_size, step_size,
                                  min_tel_frac, min_tel_len, end_margin)
    end3 = find_terminal_telomere(record.tail, False, window_size, step_size,
                                  min_tel_frac, min_tel_len, end_margin)
    return ContigQC(record, end5, end3)


def compute_assembly_stats(records: List[ContigRecord]) -> dict:
    """Compute assembly-level metrics from per-contig lengths."""
    lengths = np.array([r.length for r in records], dtype=np.int64)
    total = int(lengths.sum())
    if total == 0:
        raise ValueError("assembly has no sequence; cannot compute statistics")
    order = np.sort(lengths)[::-1]
    cum = np.cumsum(order)

    def nx(frac: float) -> Tuple[int, int]:
        idx = min(int(np.searchsorted(cum, total * frac)), len(order) - 1)
        return int(order[idx]), idx + 1

    n50, l50 = nx(0.5)
    n90, l90 = nx(0.9)
    gc_total = sum(r.gc_bases for r in records)
    n_total = sum(r.n_bases for r in records)
    non_n = total - n_total
    stats = {
        "n_contigs": len(records),
        "total_span": total,
        "longest": int(order[0]),
        "mean": float(order.mean()),
        "median": float(np.median(order)),
        "n50": n50,
        "n90": n90,
        "l50": l50,
        "l90": l90,
        "aun": float((order.astype(np.float64) ** 2).sum() / total),
        "gc_frac": (gc_total / non_n) if non_n else 0.0,
        "n_frac": n_total / total,
    }
    for threshold in (1_000, 10_000, 100_000, 1_000_000):
        mask = lengths >= threshold
        stats[f"n_contigs_ge_{threshold}"] = int(mask.sum())
        stats[f"span_ge_{threshold}"] = int(lengths[mask].sum())
    return stats


def write_contig_tsv(records: List[ContigRecord], sample_name: str,
                     output_dir: Path) -> Path:
    """Write per-contig composition table sorted by length descending."""
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / f"{sample_name}_assembly_stats.tsv"
    with open(tsv_path, "w") as f:
        f.write("contig\tlength\tgc_bases\tgc_frac\tn_bases\tn_frac\n")
        for r in sorted(records, key=lambda x: x.length, reverse=True):
            denom = r.length - r.n_bases
            gc_frac = f"{r.gc_bases / denom:.4f}" if denom else "NA"
            f.write(f"{r.name}\t{r.length}\t{r.gc_bases}\t{gc_frac}\t"
                    f"{r.n_bases}\t{r.n_bases / r.length:.4f}\n")
    logger.info(f"Wrote: {tsv_path}")
    return tsv_path


def write_assembly_summary(stats: dict, sample_name: str,
                           output_dir: Path) -> Path:
    """Write the human-readable assembly statistics summary."""
    txt_path = output_dir / f"{sample_name}_assembly_stats.txt"
    with open(txt_path, "w") as f:
        f.write(f"Sample: {sample_name}\n")
        f.write("Method: streaming contig composition scan (hifiasm_qc)\n\n")
        f.write("=== Assembly Statistics ===\n")
        f.write(f"  Contigs:            {int(stats['n_contigs']):,}\n")
        f.write(f"  Total span:         {int(stats['total_span']):,} bp\n")
        f.write(f"  Longest contig:     {int(stats['longest']):,} bp\n")
        f.write(f"  Mean contig:        {stats['mean']:,.0f} bp\n")
        f.write(f"  Median contig:      {stats['median']:,.0f} bp\n")
        f.write(f"  N50:                {int(stats['n50']):,} bp (L50 = {int(stats['l50'])})\n")
        f.write(f"  N90:                {int(stats['n90']):,} bp (L90 = {int(stats['l90'])})\n")
        f.write(f"  auN:                {stats['aun']:,.0f} bp\n")
        f.write(f"  GC (non-N bases):   {stats['gc_frac'] * 100:.2f} %\n")
        f.write(f"  N bases:            {stats['n_frac'] * 100:.3f} %\n\n")
        f.write("=== Contigs by size threshold ===\n")
        for threshold in (1_000, 10_000, 100_000, 1_000_000):
            f.write(f"  >= {threshold:>10,} bp:  "
                    f"{int(stats[f'n_contigs_ge_{threshold}']):>6,} contigs,  "
                    f"{int(stats[f'span_ge_{threshold}']):>15,} bp\n")
    logger.info(f"Wrote: {txt_path}")
    return txt_path


def write_telomere_tsv(qc_rows: List[ContigQC], sample_name: str,
                       output_dir: Path) -> Path:
    """Write per-contig terminal telomere calls sorted by length descending."""
    tsv_path = output_dir / f"{sample_name}_telomere_ends.tsv"
    with open(tsv_path, "w") as f:
        f.write("contig\tcontig_length\t"
                "tel_5p_length\ttel_5p_end_gap\ttel_5p_motif\ttel_5p_supporting_windows\t"
                "tel_3p_length\ttel_3p_end_gap\ttel_3p_motif\ttel_3p_supporting_windows\t"
                "status\n")
        for row in sorted(qc_rows, key=lambda x: x.record.length, reverse=True):
            f.write(f"{row.record.name}\t{row.record.length}\t"
                    f"{row.end5.tel_length}\t{row.end5.end_gap}\t{row.end5.motif}\t"
                    f"{row.end5.supporting_windows}\t"
                    f"{row.end3.tel_length}\t{row.end3.end_gap}\t{row.end3.motif}\t"
                    f"{row.end3.supporting_windows}\t{row.status}\n")
    logger.info(f"Wrote: {tsv_path}")
    return tsv_path


def write_telomere_summary(qc_rows: List[ContigQC], sample_name: str,
                           scan_length: int, min_tel_len: int,
                           end_margin: int, expected_chrom_ends: Optional[int],
                           output_dir: Path) -> Path:
    """Write the telomere-region assembly quality summary."""
    txt_path = output_dir / f"{sample_name}_telomere_qc.txt"
    major = [r for r in qc_rows if r.record.length >= MAJOR_CONTIG_MIN] or qc_rows
    counts = {k: 0 for k in ("both", "5prime_only", "3prime_only", "none")}
    for row in major:
        counts[row.status] += 1
    captured_ends = [(row, row.end5) for row in qc_rows if row.end5.captured] + \
                    [(row, row.end3) for row in qc_rows if row.end3.captured]
    lengths = np.array([e.tel_length for _, e in captured_ends], dtype=np.int64) \
        if captured_ends else np.array([], dtype=np.int64)

    with open(txt_path, "w") as f:
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Method: terminal telomeric repeat capture scan "
                f"(scan_length={scan_length}bp, min_tel_len={min_tel_len}bp, "
                f"end_margin={end_margin}bp)\n\n")
        f.write("=== Contig terminal status "
                f"(contigs >= {MAJOR_CONTIG_MIN:,} bp, n={len(major)}) ===\n")
        f.write(f"  Both ends telomeric (putative complete chromosome): {counts['both']}\n")
        f.write(f"  5' end only:                                        {counts['5prime_only']}\n")
        f.write(f"  3' end only:                                        {counts['3prime_only']}\n")
        f.write(f"  No terminal telomere:                               {counts['none']}\n\n")
        f.write("=== Captured telomere ends (all contigs) ===\n")
        f.write(f"  5' ends captured:  {sum(1 for row in qc_rows if row.end5.captured)}\n")
        f.write(f"  3' ends captured:  {sum(1 for row in qc_rows if row.end3.captured)}\n")
        f.write(f"  Total captured:    {len(captured_ends)} / {2 * len(qc_rows)} contig ends\n")
        if expected_chrom_ends:
            pct = len(captured_ends) / expected_chrom_ends * 100.0
            f.write(f"  Expected chromosome ends: {expected_chrom_ends}\n")
            f.write(f"  Telomere completeness:    {pct:.1f} % "
                    f"({len(captured_ends)}/{expected_chrom_ends})\n")
        if len(lengths):
            f.write("\n=== Captured telomere length distribution ===\n")
            f.write(f"  Count:   {len(lengths)}\n")
            f.write(f"  Mean:    {np.mean(lengths):.0f} bp\n")
            f.write(f"  Median:  {np.median(lengths):.0f} bp\n")
            f.write(f"  Min:     {np.min(lengths):.0f} bp\n")
            f.write(f"  Max:     {np.max(lengths):.0f} bp\n")
        else:
            f.write("\nNo terminal telomeric repeats captured. Possible causes:\n")
            f.write("  - assembly does not span chromosome ends\n")
            f.write("  - telomeric arrays truncated during assembly\n")
            f.write("  - contigs are internal fragments\n")
    logger.info(f"Wrote: {txt_path}")
    return txt_path


def plot_assembly_overview(records: List[ContigRecord], stats: dict,
                           sample_name: str, out_base: str,
                           formats: Optional[List[str]] = None) -> List[str]:
    """4-panel assembly overview figure (multi-format).

    Panels: (1) contig length distribution; (2) Nx curve; (3) GC vs length;
    (4) rank-size cumulative span curve.
    """
    lengths = np.sort(np.array([r.length for r in records], dtype=np.int64))[::-1]
    total = int(lengths.sum())

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0][0]
    ax.hist(np.log10(lengths), bins=40, color="#3C6E9F", edgecolor="none")
    ax.axvline(np.log10(stats["n50"]), color="#4C9F70", linestyle="--",
               label=f"N50 = {int(stats['n50']):,} bp")
    ax.axvline(np.log10(stats["median"]), color="#C0504D", linestyle="--",
               label=f"median = {int(stats['median']):,} bp")
    ax.set_xlabel("log10 contig length (bp)")
    ax.set_ylabel("Contigs")
    ax.set_title("Contig length distribution")
    ax.legend(fontsize=8)

    ax = axes[0][1]
    cum_frac = np.cumsum(lengths) / total * 100.0
    ax.plot(cum_frac, lengths, color="#3C6E9F", linewidth=1.5)
    for frac, value, tag in ((50, stats["n50"], "N50"), (90, stats["n90"], "N90")):
        ax.plot([frac], [value], "o", color="#4C9F70")
        ax.annotate(f"{tag} = {int(value):,} bp", (frac, value),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("Cumulative assembly span (%)")
    ax.set_ylabel("Nx (bp)")
    ax.set_title("Nx curve")

    ax = axes[1][0]
    plot_records = records
    if len(plot_records) > 50_000:
        plot_records = plot_records[::len(plot_records) // 50_000]
    gc_vals = []
    len_vals = []
    for r in plot_records:
        denom = r.length - r.n_bases
        if denom > 0:
            gc_vals.append(r.gc_bases / denom * 100.0)
            len_vals.append(r.length)
    ax.scatter(np.log10(len_vals), gc_vals, s=6, alpha=0.4,
               color="#3C6E9F", edgecolors="none")
    ax.set_xlabel("log10 contig length (bp)")
    ax.set_ylabel("GC (%)")
    ax.set_title("GC content vs contig length")
    ax.set_ylim(0, 100)

    ax = axes[1][1]
    ranks = np.arange(1, len(lengths) + 1)
    ax.plot(ranks, np.cumsum(lengths) / total * 100.0, color="#3C6E9F",
            linewidth=1.5)
    ax.set_xscale("log")
    ax.set_xlabel("Contig rank (by length)")
    ax.set_ylabel("Cumulative span (%)")
    ax.set_title("Rank-size curve (assembly continuity)")

    fig.suptitle(f"Assembly quality overview — {sample_name}", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    written = save_fig(fig, out_base, formats)
    plt.close(fig)
    return written


def plot_telomere_qc(qc_rows: List[ContigQC], sample_name: str,
                     min_tel_len: int, end_margin: int,
                     expected_chrom_ends: Optional[int], out_base: str,
                     formats: Optional[List[str]] = None) -> List[str]:
    """4-panel telomere-region assembly quality figure (multi-format).

    Panels: (1) contig terminal status counts; (2) captured telomere length
    distribution; (3) per-contig telomere landscape; (4) telomere completeness.
    """
    major = [r for r in qc_rows if r.record.length >= MAJOR_CONTIG_MIN] or qc_rows
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0][0]
    cats = ["both", "5prime_only", "3prime_only", "none"]
    labels = ["both ends", "5' only", "3' only", "no telomere"]
    colors = ["#4C9F70", "#3C6E9F", "#C0504D", "#999999"]
    counts = [sum(1 for r in major if r.status == c) for c in cats]
    ax.bar(labels, counts, color=colors, edgecolor="none")
    for i, count in enumerate(counts):
        ax.text(i, count, str(count), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Contigs")
    ax.set_title(f"Contig terminal status (>= {MAJOR_CONTIG_MIN // 1_000_000} Mb, n={len(major)})")

    ax = axes[0][1]
    cap5 = [r.end5.tel_length for r in qc_rows if r.end5.captured]
    cap3 = [r.end3.tel_length for r in qc_rows if r.end3.captured]
    if cap5 or cap3:
        ax.hist([cap5, cap3], bins=30, stacked=True,
                color=["#3C6E9F", "#C0504D"], edgecolor="none",
                label=["5' ends", "3' ends"])
        all_len = np.array(cap5 + cap3)
        ax.axvline(np.median(all_len), color="#4C9F70", linestyle="--",
                   label=f"median = {int(np.median(all_len)):,} bp")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no terminal telomeres captured",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_xlabel("Captured telomere length (bp)")
    ax.set_ylabel("Contig ends")
    ax.set_title(f"Telomere length at contig ends (>= {min_tel_len:,} bp)")

    ax = axes[1][0]
    landscape = [r for r in qc_rows if r.record.length >= 100_000][:200]
    if landscape:
        ranks = range(1, len(landscape) + 1)
        xs5 = [i for i, r in zip(ranks, landscape) if r.end5.captured]
        ys5 = [r.end5.tel_length for r in landscape if r.end5.captured]
        xs3 = [i for i, r in zip(ranks, landscape) if r.end3.captured]
        ys3 = [r.end3.tel_length for r in landscape if r.end3.captured]
        if xs5:
            ax.scatter(xs5, ys5, s=14, c="#3C6E9F", marker="o",
                       label="5' end", alpha=0.8)
        if xs3:
            ax.scatter(xs3, ys3, s=14, c="#C0504D", marker="^",
                       label="3' end", alpha=0.8)
        ax.set_xlabel("Contig rank by length (>= 100 kb)")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no contigs >= 100 kb",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_ylabel("Telomere length (bp)")
    ax.set_title("Telomere landscape across contigs")

    ax = axes[1][1]
    captured = sum(1 for r in qc_rows if r.end5.captured) + \
               sum(1 for r in qc_rows if r.end3.captured)
    reference = expected_chrom_ends if expected_chrom_ends else 2 * len(qc_rows)
    missing = max(reference - captured, 0)
    ax.barh(["telomere capture"], [captured], color="#4C9F70",
            edgecolor="none", label=f"captured ({captured})")
    ax.barh(["telomere capture"], [missing], left=[captured], color="#E0E0E0",
            edgecolor="none",
            label=f"{'missing' if expected_chrom_ends else 'uncaptured'} ({missing})")
    ax.set_xlim(0, reference)
    ax.set_xlabel("Chromosome ends")
    pct = captured / reference * 100.0 if reference else 0.0
    ax.set_title(f"Telomere completeness: {pct:.1f}% "
                 f"(end margin {end_margin:,} bp)")
    ax.legend(fontsize=8, loc="lower right")

    fig.suptitle(f"Telomere-region assembly quality — {sample_name}", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    written = save_fig(fig, out_base, formats)
    plt.close(fig)
    return written


def plot_telomere_map(qc_rows: List[ContigQC], sample_name: str,
                      top_contigs: int, out_base: str,
                      formats: Optional[List[str]] = None) -> List[str]:
    """Telomere capture map: top contigs as bars with colored end caps.

    Green caps mark captured terminal telomeres (with length labels); red
    caps mark contig ends without a captured telomere.
    """
    top = sorted(qc_rows, key=lambda x: x.record.length, reverse=True)[:top_contigs]
    if not top:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "no contigs to draw", ha="center", va="center")
        ax.set_axis_off()
        written = save_fig(fig, out_base, formats)
        plt.close(fig)
        return written

    max_len = top[0].record.length
    cap_frac = 0.015
    fig, ax = plt.subplots(figsize=(13, max(3.5, 0.42 * len(top) + 1.5)))
    pad = max_len * 0.10

    for i, row in enumerate(top):
        y = len(top) - i - 1
        ax.barh(y, row.record.length, height=0.55, color="#D6E0F0",
                edgecolor="none")
        for x_pos, call, align in ((0, row.end5, "right"),
                                   (row.record.length, row.end3, "left")):
            cap_w = max_len * cap_frac
            offset = 0 if x_pos == 0 else x_pos - cap_w
            if call.captured:
                ax.barh(y, cap_w, left=offset, height=0.7, color="#4C9F70",
                        edgecolor="none")
                label = f"{call.tel_length / 1000:.1f}kb"
                label_color = "#2E7D32"
            else:
                ax.barh(y, cap_w, left=offset, height=0.7, color="#C0504D",
                        edgecolor="none")
                label = "none"
                label_color = "#B71C1C"
            label_x = -pad * 0.08 if align == "right" else row.record.length + pad * 0.08
            ax.text(label_x, y, label, ha=align, va="center", fontsize=7,
                    color=label_color)

    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([f"{r.record.name} ({r.record.length / 1e6:.1f} Mb)"
                        for r in reversed(top)], fontsize=8)
    ax.set_xlim(-pad, max_len + pad)
    ax.set_xlabel("Contig length (bp)")
    ax.set_title(f"Telomere capture map — top {len(top)} contigs "
                 f"({sample_name})")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#4C9F70", label="captured telomere"),
                       Patch(color="#C0504D", label="no telomere"),
                       Patch(color="#D6E0F0", label="contig body")],
              fontsize=8, loc="lower right")
    fig.tight_layout()
    written = save_fig(fig, out_base, formats)
    plt.close(fig)
    return written


def main() -> None:
    """Parse CLI arguments and run scan / statistics / plots in order."""
    global logger
    parser = argparse.ArgumentParser(
        description="Assembly QC: statistics + terminal telomere capture "
                    "from a hifiasm contig FASTA")
    parser.add_argument("-f", "--fasta", required=True,
                        help="Input assembly FASTA (e.g. asm.bp.p_ctg.fa)")
    parser.add_argument("-s", "--sample_name", required=True,
                        help="Sample name for output files")
    parser.add_argument("-o", "--output_dir", default="results",
                        help="Output directory")
    parser.add_argument("-L", "--scan_length", type=int,
                        default=DEFAULT_SCAN_LENGTH,
                        help=f"bp scanned from each contig end "
                             f"(default: {DEFAULT_SCAN_LENGTH})")
    parser.add_argument("-w", "--window_size", type=int,
                        default=DEFAULT_WINDOW_SIZE,
                        help=f"sliding window size (default: {DEFAULT_WINDOW_SIZE})")
    parser.add_argument("-e", "--step_size", type=int,
                        default=DEFAULT_STEP_SIZE,
                        help=f"sliding window step (default: {DEFAULT_STEP_SIZE})")
    parser.add_argument("-t", "--min_tel_frac", type=float,
                        default=DEFAULT_MIN_TEL_FRAC,
                        help=f"min telomeric fraction per window "
                             f"(default: {DEFAULT_MIN_TEL_FRAC})")
    parser.add_argument("-m", "--min_tel_len", type=int,
                        default=DEFAULT_MIN_TEL_LEN,
                        help=f"min captured telomere length in bp "
                             f"(default: {DEFAULT_MIN_TEL_LEN})")
    parser.add_argument("-g", "--end_margin", type=int,
                        default=DEFAULT_END_MARGIN,
                        help=f"max gap from contig end for a terminal call "
                             f"(default: {DEFAULT_END_MARGIN})")
    parser.add_argument("--expected-chrom-ends", type=int, default=0,
                        help="expected chromosome ends (e.g. 42 for a "
                             "haploid 21-chromosome assembly); 0 disables "
                             "the completeness comparison")
    parser.add_argument("--top-contigs", type=int, default=DEFAULT_TOP_CONTIGS,
                        help=f"contigs drawn in the capture map "
                             f"(default: {DEFAULT_TOP_CONTIGS})")
    parser.add_argument("-l", "--log", default="",
                        help="Optional log file path (default: log to stdout only)")
    parser.add_argument("-F", "--format", action="append", dest="formats",
                        metavar="FMT",
                        help="Output figure format(s), repeatable (default: png)")
    args = parser.parse_args()

    if args.log:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        logger = setup_logger("hifiasm_qc", log_file=args.log, stream=sys.stdout)

    formats = normalize_formats(args.formats)
    output_dir = Path(args.output_dir)

    records: List[ContigRecord] = []
    qc_rows: List[ContigQC] = []
    for record in iter_contig_records(args.fasta, args.scan_length):
        records.append(record)
        qc_rows.append(scan_contig(record, args.window_size, args.step_size,
                                   args.min_tel_frac, args.min_tel_len,
                                   args.end_margin))
    if not records:
        raise ValueError(f"no contigs parsed from {args.fasta}; check the FASTA")
    logger.info(f"Scanned {len(records)} contigs "
                f"(scan_length={args.scan_length}bp, "
                f"min_tel_len={args.min_tel_len}bp, "
                f"end_margin={args.end_margin}bp)")

    stats = compute_assembly_stats(records)
    write_contig_tsv(records, args.sample_name, output_dir)
    write_assembly_summary(stats, args.sample_name, output_dir)
    write_telomere_tsv(qc_rows, args.sample_name, output_dir)
    write_telomere_summary(qc_rows, args.sample_name, args.scan_length,
                           args.min_tel_len, args.end_margin,
                           args.expected_chrom_ends or None, output_dir)

    for path in plot_assembly_overview(
            records, stats, args.sample_name,
            str(output_dir / f"{args.sample_name}_assembly_overview"), formats):
        logger.info(f"Wrote: {path}")
    for path in plot_telomere_qc(
            qc_rows, args.sample_name, args.min_tel_len, args.end_margin,
            args.expected_chrom_ends or None,
            str(output_dir / f"{args.sample_name}_telomere_qc"), formats):
        logger.info(f"Wrote: {path}")
    for path in plot_telomere_map(
            qc_rows, args.sample_name, args.top_contigs,
            str(output_dir / f"{args.sample_name}_telomere_map"), formats):
        logger.info(f"Wrote: {path}")


if __name__ == "__main__":
    main()
