#!/usr/bin/env python3
"""Telomere length + repeat-composition analysis for Fiber-seq samples.

Per sample, reads from a fire.sorted.bam are scanned (no sequence needed:
uses the m6a/msp BED12 fiber spans for read names plus BAM SEQ):

  - terminal telomeric repeat length (5'/3') via the sliding-window
    (TTAGGG|CCCTAA) detector (same rules as src/telomere/analyze_telomere.py)
  - repeat-composition variants: counts of canonical TTAGGG, variant motifs
    (TCAGGG, TTGGGG, TGTGGG, ...), and total telomeric bp per read -> per
    sample composition table (answers "telomere composition changes")
  - internal telomeric sequence (ITS) length
  - per-sample summary (median/mean telomere length, fraction of tel+ reads)
    + two-sample comparison table + length-distribution violin plot

Usage:
  python telomere_profile.py -s DMSO=DMSO.fire.sorted.bam \
      -s PlaB=PlaB.fire.sorted.bam -o tel_out/ [--region chr1:0-500000]
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pysam

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402
from src.common.util.FiberlibUtil import write_tsv  # noqa: E402

CANONICAL = "TTAGGG"
VARIANTS = ["TCAGGG", "TTGGGG", "TGTGGG", "TTAGGC", "TTAGAG", "CTAGGG",
            "TTGGGC", "TTCGGG", "TTAAGG", "TTACGG"]
WINDOW = 100
STEP = 20
MIN_TEL_CONTENT = 0.6
MOTIF_RE = re.compile("|".join([CANONICAL] + VARIANTS), re.IGNORECASE)


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("TelomereProfile", log_file=str(log_file))
    else:
        lg = setup_logger("TelomereProfile")
    for h in lg.handlers:
        if type(h) is logging.StreamHandler:
            h.setStream(sys.stdout)
    return lg


logger = _get_logger()


def tel_windows(seq: str) -> List[Tuple[int, int]]:
    """Terminal-only sliding windows with >= MIN_TEL_CONTENT telomeric 6-mers.

    Only the first/last --edge windows of the read are scanned (terminal
    telomeres); a whole-read regex prefilter rejects tel-negative reads
    cheaply.
    """
    hits: List[Tuple[int, int]] = []
    edge = min(len(seq), 30000)  # telomeres never exceed this in mouse
    tel_motifs = {CANONICAL, *VARIANTS}
    # cheap prefilter: any 3 canonical/variant motifs anywhere near the ends
    head, tail = seq[:edge].upper(), seq[-edge:].upper()
    if sum(head.count(m) for m in tel_motifs) + sum(tail.count(m) for m in tel_motifs) < 3:
        return hits
    for region in (head, tail):
        for start in range(0, max(0, len(region) - WINDOW + 1), STEP):
            win = region[start:start + WINDOW]
            n_motif = sum(1 for i in range(0, len(win) - 5, 3)
                          if win[i:i + 6] in tel_motifs)
            if n_motif * 6 >= MIN_TEL_CONTENT * WINDOW:
                hits.append((start, start + WINDOW))
    return hits


def merge_windows(hits: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged = []
    for s, e in hits:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def composition(seq: str) -> Dict[str, int]:
    """Count canonical + variant motifs in a telomeric stretch."""
    counts: Dict[str, int] = defaultdict(int)
    for m in MOTIF_RE.finditer(seq.upper()):
        counts[m.group().upper()] += 1
    return counts


def analyze_read(seq: str, min_len: int) -> Optional[dict]:
    """Terminal telomere lengths + composition for one read."""
    if len(seq) < min_len:
        return None
    hits = merge_windows(tel_windows(seq))
    if not hits:
        return None
    rlen = len(seq)
    p5 = p3 = internal = 0
    comp5: Dict[str, int] = defaultdict(int)
    comp3: Dict[str, int] = defaultdict(int)
    for s, e in hits:
        stretch = seq[s:e]
        if s <= rlen * 0.05:
            p5 = max(p5, e)
            for k, v in composition(stretch).items():
                comp5[k] += v
        elif e >= rlen * 0.95:
            p3 = max(p3, e - s)
            for k, v in composition(stretch).items():
                comp3[k] += v
        else:
            internal += e - s
    if p5 + p3 == 0:
        return None
    return {"len5": p5, "len3": p3, "internal": internal,
            "comp5": dict(comp5), "comp3": dict(comp3)}


def analyze_bam(bam_path: str, label: str, min_len: int, max_reads: int,
                region: Optional[str]) -> List[dict]:
    """Terminal-telomere stats for all (or region-restricted) reads."""
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        fetch = bam.fetch(region=region) if region else bam.fetch(until_eof=True)
        for n, read in enumerate(fetch):
            if n >= max_reads:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            seq = read.query_sequence
            if seq is None:
                continue
            res = analyze_read(seq, min_len)
            if res is None:
                continue
            rows.append({"sample": label, "read": read.query_name,
                         "read_len": len(seq), "chrom": read.reference_name,
                         **{k: v for k, v in res.items() if not k.startswith("comp")},
                         "comp5": res["comp5"], "comp3": res["comp3"]})
    logger.info(f"{label}: {len(rows)} telomere-positive reads")
    return rows


def summarize(rows: List[dict]) -> List[dict]:
    """Per-sample length summary + motif composition fractions."""
    per = defaultdict(list)
    for r in rows:
        per[r["sample"]].append(r)
    out = []
    for sample, rs in sorted(per.items()):
        t5 = np.array([r["len5"] for r in rs])
        t3 = np.array([r["len3"] for r in rs])
        tall = np.concatenate([t5, t3])
        comp = defaultdict(int)
        for r in rs:
            for k, v in r["comp5"].items():
                comp[k] += v
            for k, v in r["comp3"].items():
                comp[k] += v
        tot = sum(comp.values()) or 1
        base = {"sample": sample, "n_reads": len(rs),
                "tel5_median": float(np.median(t5)) if len(t5) else 0.0,
                "tel3_median": float(np.median(t3)) if len(t3) else 0.0,
                "tel_all_median": float(np.median(tall)) if len(tall) else 0.0,
                "tel_all_mean": float(np.mean(tall)) if len(tall) else 0.0,
                "tel_all_q10": float(np.percentile(tall, 10)) if len(tall) else 0.0,
                "tel_all_q90": float(np.percentile(tall, 90)) if len(tall) else 0.0,
                "internal_median": float(np.median([r["internal"] for r in rs]))}
        for motif in [CANONICAL] + VARIANTS:
            base[f"motif_{motif}_frac"] = round(comp.get(motif, 0) / tot, 5)
        out.append(base)
    return out


def plot_lengths(rows: List[dict], outdir: Path, formats: Sequence[str],
                 cap: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    per = defaultdict(list)
    for r in rows:
        per[r["sample"]].extend([r["len5"], r["len3"]])
    samples = sorted(per)
    fig, ax = plt.subplots(figsize=(1.8 * max(1, len(samples)) + 2, 4))
    data = [np.clip(np.array(per[s]), 0, cap) for s in samples]
    parts = ax.violinplot(data, showmedians=True)
    ax.set_xticks(range(1, len(samples) + 1))
    ax.set_xticklabels(samples)
    ax.set_ylabel(f"terminal telomere length (bp, capped {cap})")
    for i, d in enumerate(data):
        ax.text(i + 1, np.median(d), f"{np.median(d):.0f}", ha="center",
                va="bottom", fontsize=8)
    ax.set_title("Telomere length by sample")
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(outdir / f"telomere_length_violin.{fmt}", dpi=300)
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="telomere_profile.py",
        description="Telomere length + repeat composition from fire BAMs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-s", "--sample", action="append", required=True,
                   metavar="LABEL=BAM", help="sample label + fire.sorted.bam")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("--min-read-len", type=int, default=5000)
    p.add_argument("--max-reads", type=int, default=2_000_000)
    p.add_argument("--region", default=None,
                   help="restrict analysis (e.g. chr1:0-1000000); default whole BAM")
    p.add_argument("--cap", type=int, default=30000, help="violin cap (bp)")
    p.add_argument("--no-plot", action="store_true")
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

    rows: List[dict] = []
    for spec in args.sample:
        label, _, bam = spec.partition("=")
        rows.extend(analyze_bam(bam, label, args.min_read_len, args.max_reads,
                                args.region))

    per_read_header = ["sample", "read", "chrom", "read_len", "len5", "len3", "internal"]
    write_tsv(str(args.outdir / "telomere_per_read.tsv"), per_read_header,
              [[{k: r.get(k, "") for k in per_read_header}[c] for c in per_read_header]
               for r in rows])

    summary = summarize(rows)
    sheader = sorted({k for r in summary for k in r})
    write_tsv(str(args.outdir / "telomere_summary.tsv"), sheader,
              [[r.get(c, "") for c in sheader] for r in summary])

    # composition per sample side table
    comp_rows = []
    for r in summary:
        for motif in [CANONICAL] + VARIANTS:
            comp_rows.append((r["sample"], motif, r.get(f"motif_{motif}_frac", 0)))
    write_tsv(str(args.outdir / "telomere_composition.tsv"),
              ["sample", "motif", "fraction"], comp_rows)

    if not args.no_plot and rows:
        plot_lengths(rows, args.outdir, args.formats, args.cap)
    logger.info(f"done -> {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
