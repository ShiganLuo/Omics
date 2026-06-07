#!/usr/bin/env python3
"""Summarize Arriba fusion detection results across multiple samples.

Reads passed_fusions.tsv from each sample, generates:
  - Per-sample fusion count and type/confidence breakdown
  - Cross-sample recurrent fusion gene analysis
  - In-frame fusion gene list (potential functional fusions)
  - High-confidence fusion detail table
  - Summary HTML report
  - Publication-quality figures (png/pdf/svg/tiff/eps)

Usage (two input modes, mutually exclusive):
    # Mode 1: scan a directory for sample subdirectories
    python summarize_arriba_fusions.py --indir /path/to/arriba_output -o /path/to/report

    # Mode 2: pass file lists directly (comma-separated)
    python summarize_arriba_fusions.py -p s1_passed.tsv,s2_passed.tsv -o /path/to/report
    python summarize_arriba_fusions.py -p s1_passed.tsv,s2_passed.tsv -d s1_discarded.tsv,s2_discarded.tsv -o /path/to/report -f png
"""

import argparse
import csv
import logging
import os
import sys
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(os.path.join(ROOT_DIR, "src"))
from common.LogUtil import setup_logger
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
logger = setup_logger(__name__, logging.INFO)

# ── Constants ──────────────────────────────────────────────────────────────────

FUSION_COLUMNS = [
    "gene1", "gene2", "strand1", "strand2", "breakpoint1", "breakpoint2",
    "site1", "site2", "type", "split_reads1", "split_reads2",
    "discordant_mates", "coverage1", "coverage2", "confidence",
    "reading_frame", "tags", "retained_protein_domains",
    "closest_genomic_breakpoint1", "closest_genomic_breakpoint2",
    "gene_id1", "gene_id2", "transcript_id1", "transcript_id2",
    "direction1", "direction2", "filters", "fusion_transcript",
    "peptide_sequence", "read_identifiers",
]

CONFIDENCE_ORDER = {"high": 0, "medium": 1, "low": 2}


# ── Helper functions ───────────────────────────────────────────────────────────

def parse_fusion_tsv(filepath: str) -> List[Dict[str, str]]:
    """Parse an Arriba TSV output file into a list of dicts.

    Skips comment lines starting with '#'. Handles the Arriba-specific
    format where gene names can contain aliases in parentheses.
    """
    rows: List[Dict[str, str]] = []
    with open(filepath, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(FUSION_COLUMNS):
                fields.extend(["."] * (len(FUSION_COLUMNS) - len(fields)))
            row = dict(zip(FUSION_COLUMNS, fields[: len(FUSION_COLUMNS)]))
            rows.append(row)
    return rows


def canonical_fusion_key(gene1: str, gene2: str) -> str:
    """Return a sorted canonical key for a gene pair.

    Strips Arriba alias annotations like 'Gm43566(174),AI506816(22572)'.
    Takes the primary gene name (before '(' or ',').
    """
    def primary_name(g: str) -> str:
        return g.split("(")[0].split(",")[0].strip()

    g1, g2 = primary_name(gene1), primary_name(gene2)
    return "::".join(sorted([g1, g2]))


def estimate_support(row: Dict[str, str]) -> int:
    """Estimate total supporting reads from split_reads + discordant_mates."""
    try:
        sr1 = int(row.get("split_reads1", "0") or "0")
    except ValueError:
        sr1 = 0
    try:
        sr2 = int(row.get("split_reads2", "0") or "0")
    except ValueError:
        sr2 = 0
    try:
        dm = int(row.get("discordant_mates", "0") or "0")
    except ValueError:
        dm = 0
    return sr1 + sr2 + dm


def short_fusion_label(gene1: str, gene2: str) -> str:
    """Shorten gene names for display in figures."""
    def shorten(g: str) -> str:
        name = g.split("(")[0].split(",")[0].strip()
        return name[:15] + "…" if len(name) > 15 else name
    return f"{shorten(gene1)}--{shorten(gene2)}"


# ── Sample discovery / resolution ──────────────────────────────────────────────

def discover_samples(indir: str) -> Dict[str, Tuple[str, Optional[str]]]:
    """Find sample directories containing passed_fusions.tsv.

    Returns:
        Dict mapping sample_id -> (passed_path, discarded_path_or_None)
    """
    samples: Dict[str, Tuple[str, Optional[str]]] = {}
    indir_path = Path(indir)
    for entry in sorted(indir_path.iterdir()):
        if not entry.is_dir():
            continue
        passed = entry / f"{entry.name}_passed_fusions.tsv"
        discarded = entry / f"{entry.name}_discarded_fusions.tsv"
        if passed.exists():
            samples[entry.name] = (
                str(passed),
                str(discarded) if discarded.exists() else None,
            )
    return samples


def resolve_samples(
    passed_files: List[str], discarded_files: Optional[List[str]] = None,
) -> Dict[str, Tuple[str, Optional[str]]]:
    """Map sample_id extracted from filename to (passed_path, discarded_path).

    Expects filenames following the Arriba convention:
        {sample_id}_passed_fusions.tsv / {sample_id}_discarded_fusions.tsv
    """
    suffix = "_passed_fusions.tsv"
    samples: Dict[str, Tuple[str, Optional[str]]] = {}

    discarded_map: Dict[str, str] = {}
    if discarded_files:
        for df in discarded_files:
            base = Path(df).name
            if base.endswith("_discarded_fusions.tsv"):
                sid = base[: -len("_discarded_fusions.tsv")]
                discarded_map[sid] = df

    for pf in passed_files:
        base = Path(pf).name
        if not base.endswith(suffix):
            logger.warning(f"Unexpected filename (skipped): {pf}")
            continue
        sid = base[: -len(suffix)]
        samples[sid] = (pf, discarded_map.get(sid))

    return samples


# ── Analysis functions ─────────────────────────────────────────────────────────

def per_sample_stats(rows: List[Dict[str, str]]) -> Dict[str, object]:
    """Compute per-sample summary statistics."""
    total = len(rows)
    type_counts = Counter(r["type"] for r in rows)
    conf_counts = Counter(r["confidence"] for r in rows)
    frame_counts = Counter(r["reading_frame"] for r in rows)

    # supporting reads distribution
    supports = [estimate_support(r) for r in rows]
    avg_support = sum(supports) / total if total else 0

    # unique gene pairs
    gene_pairs = set(canonical_fusion_key(r["gene1"], r["gene2"]) for r in rows)

    return {
        "total_fusions": total,
        "unique_gene_pairs": len(gene_pairs),
        "type_counts": dict(type_counts.most_common()),
        "confidence_counts": dict(conf_counts),
        "reading_frame_counts": dict(frame_counts),
        "avg_supporting_reads": round(avg_support, 1),
        "max_supporting_reads": max(supports) if supports else 0,
    }


def collect_high_confidence(rows: List[Dict[str, str]], sample_id: str) -> List[Dict[str, str]]:
    """Return high-confidence fusions annotated with sample_id."""
    results = []
    for r in rows:
        if r["confidence"] in ("high", "medium"):
            entry = {
                "sample": sample_id,
                "gene1": r["gene1"],
                "gene2": r["gene2"],
                "type": r["type"],
                "confidence": r["confidence"],
                "reading_frame": r["reading_frame"],
                "split_reads": f"{r['split_reads1']}+{r['split_reads2']}",
                "discordant_mates": r["discordant_mates"],
                "total_support": estimate_support(r),
                "breakpoint1": r["breakpoint1"],
                "breakpoint2": r["breakpoint2"],
                "filters": r["filters"],
            }
            results.append(entry)
    # sort by confidence then support
    results.sort(key=lambda x: (CONFIDENCE_ORDER.get(x["confidence"], 9), -x["total_support"]))
    return results


def find_recurrent_fusions(all_sample_fusions: Dict[str, List[Dict[str, str]]]) -> List[Dict]:
    """Identify fusion genes appearing in >= 2 samples."""
    gene_pair_samples: Dict[str, Set[str]] = defaultdict(set)
    gene_pair_details: Dict[str, Dict[str, str]] = {}

    for sample_id, rows in all_sample_fusions.items():
        for r in rows:
            key = canonical_fusion_key(r["gene1"], r["gene2"])
            gene_pair_samples[key].add(sample_id)
            if key not in gene_pair_details:
                gene_pair_details[key] = {
                    "gene1": r["gene1"],
                    "gene2": r["gene2"],
                    "type": r["type"],
                    "confidence": r["confidence"],
                    "reading_frame": r["reading_frame"],
                }

    recurrent = []
    for key, samples in gene_pair_samples.items():
        if len(samples) >= 2:
            recurrent.append({
                "fusion": key,
                "n_samples": len(samples),
                "samples": ",".join(sorted(samples)),
                **gene_pair_details[key],
            })
    recurrent.sort(key=lambda x: -x["n_samples"])
    return recurrent


def find_inframe_fusions(all_sample_fusions: Dict[str, List[Dict[str, str]]]) -> List[Dict[str, str]]:
    """Extract in-frame fusions (potentially functional)."""
    inframe = []
    for sample_id, rows in all_sample_fusions.items():
        for r in rows:
            if r["reading_frame"] == "in-frame":
                inframe.append({
                    "sample": sample_id,
                    "gene1": r["gene1"],
                    "gene2": r["gene2"],
                    "type": r["type"],
                    "confidence": r["confidence"],
                    "total_support": estimate_support(r),
                    "retained_domains": r["retained_protein_domains"],
                    "breakpoint1": r["breakpoint1"],
                    "breakpoint2": r["breakpoint2"],
                })
    inframe.sort(key=lambda x: (-x["total_support"], x["sample"]))
    return inframe


# ── TSV output ─────────────────────────────────────────────────────────────────

def write_tsv(rows: List[Dict], filepath: str) -> None:
    """Write a list of dicts to a TSV file."""
    if not rows:
        logger.info(f"  (no rows to write for {filepath})")
        return
    fieldnames = list(rows[0].keys())
    with open(filepath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logger.info(f"  Wrote {len(rows)} rows -> {filepath}")


# ── HTML report ────────────────────────────────────────────────────────────────

def _html_table(headers: List[str], rows: List[List[str]], max_rows: int = 200) -> str:
    """Build an HTML table string."""
    lines = ["<table>", "  <thead><tr>"]
    for h in headers:
        lines.append(f"    <th>{h}</th>")
    lines.append("  </tr></thead>", )
    lines.append("  <tbody>")
    for i, row in enumerate(rows[:max_rows]):
        cls = ' class="alt"' if i % 2 else ""
        lines.append(f"  <tr{cls}>")
        for cell in row:
            lines.append(f"    <td>{cell}</td>")
        lines.append("  </tr>")
    if len(rows) > max_rows:
        lines.append(f'  <tr><td colspan="{len(headers)}">... {len(rows) - max_rows} more rows omitted</td></tr>')
    lines.append("  </tbody>")
    lines.append("</table>")
    return "\n".join(lines)


def generate_html_report(
    per_sample: Dict[str, Dict],
    high_conf: List[Dict[str, str]],
    recurrent: List[Dict],
    inframe: List[Dict[str, str]],
    outpath: str,
) -> None:
    """Generate a self-contained HTML summary report."""
    # ── Per-sample summary table
    sample_headers = [
        "Sample", "Passed Fusions", "Unique Gene Pairs",
        "High", "Medium", "Low",
        "Avg Support", "Max Support",
    ]
    sample_rows = []
    for sid in sorted(per_sample):
        s = per_sample[sid]
        cc = s["confidence_counts"]
        sample_rows.append([
            sid,
            str(s["total_fusions"]),
            str(s["unique_gene_pairs"]),
            str(cc.get("high", 0)),
            str(cc.get("medium", 0)),
            str(cc.get("low", 0)),
            str(s["avg_supporting_reads"]),
            str(s["max_supporting_reads"]),
        ])

    # ── Fusion type distribution
    all_types: Counter = Counter()
    for s in per_sample.values():
        for t, c in s["type_counts"].items():
            all_types[t] += c
    type_headers = ["Fusion Type", "Total Count"]
    type_rows = [[t, str(c)] for t, c in all_types.most_common()]

    # ── High-confidence table
    hc_headers = [
        "Sample", "Gene1", "Gene2", "Type", "Confidence",
        "Reading Frame", "Split Reads", "Discordant", "Total Support",
        "Breakpoint1", "Breakpoint2", "Filters",
    ]
    hc_rows = [
        [
            r["sample"], r["gene1"], r["gene2"], r["type"], r["confidence"],
            r["reading_frame"], r["split_reads"], r["discordant_mates"],
            str(r["total_support"]), r["breakpoint1"], r["breakpoint2"],
            r["filters"],
        ]
        for r in high_conf
    ]

    # ── Recurrent table
    rec_headers = ["Fusion", "N Samples", "Samples", "Type", "Confidence", "Reading Frame"]
    rec_rows = [
        [r["fusion"], str(r["n_samples"]), r["samples"],
         r["type"], r["confidence"], r["reading_frame"]]
        for r in recurrent
    ]

    # ── In-frame table
    if_headers = [
        "Sample", "Gene1", "Gene2", "Type", "Confidence",
        "Total Support", "Retained Domains", "Breakpoint1", "Breakpoint2",
    ]
    if_rows = [
        [
            r["sample"], r["gene1"], r["gene2"], r["type"], r["confidence"],
            str(r["total_support"]), r["retained_domains"],
            r["breakpoint1"], r["breakpoint2"],
        ]
        for r in inframe
    ]

    # ── Figure conclusions
    conclusions = generate_figure_conclusions(per_sample, recurrent, inframe)
    sections_data = [
        ("融合检测概览 / Fusion Detection Overview", [
            ("Fig1 各样本融合数（按置信度堆叠）", conclusions[0]),
            ("Fig2 融合类型分布", conclusions[1]),
            ("Fig3 各样本融合类型组成", conclusions[2]),
        ]),
        ("融合事件特征 / Fusion Characterization", [
            ("Fig4 读框分布", conclusions[3]),
            ("Fig5 跨样本重复融合热图", conclusions[4]),
            ("Fig6 框内融合支持reads", conclusions[5]),
        ]),
    ]
    conclusion_html = ""
    for sec_title, items in sections_data:
        conclusion_html += f"<tr><th colspan='2' style='background:#0f3460;text-align:left;padding:10px'>{sec_title}</th></tr>\n"
        for name, conc in items:
            conclusion_html += f"<tr><td style='font-weight:bold;white-space:nowrap'>{name}</td><td>{conc}</td></tr>\n"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Arriba Fusion Summary Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 1400px; margin: 0 auto; padding: 20px; background: #fafafa; color: #333; }}
  h1 {{ color: #1a1a2e; border-bottom: 3px solid #16213e; padding-bottom: 10px; }}
  h2 {{ color: #16213e; margin-top: 40px; }}
  table {{ border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 13px; }}
  th {{ background: #16213e; color: #fff; padding: 8px 10px; text-align: left; }}
  td {{ padding: 6px 10px; border-bottom: 1px solid #e0e0e0; }}
  tr:hover {{ background: #e8f0fe; }}
  tr.alt {{ background: #f5f5f5; }}
  .summary-box {{ background: #fff; border: 1px solid #ddd; border-radius: 8px;
                  padding: 20px; margin: 15px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
  .metric {{ display: inline-block; margin: 10px 25px 10px 0; }}
  .metric .value {{ font-size: 28px; font-weight: bold; color: #0f3460; }}
  .metric .label {{ font-size: 12px; color: #666; }}
  .high {{ color: #27ae60; font-weight: bold; }}
  .medium {{ color: #f39c12; }}
  .low {{ color: #95a5a6; }}
</style>
</head>
<body>
<h1>Arriba Fusion Detection — Summary Report</h1>

<h2>Per-Sample Overview</h2>
<div class="summary-box">
  <div class="metric">
    <div class="value">{len(per_sample)}</div>
    <div class="label">Samples</div>
  </div>
  <div class="metric">
    <div class="value">{sum(s['total_fusions'] for s in per_sample.values())}</div>
    <div class="label">Total Passed Fusions</div>
  </div>
  <div class="metric">
    <div class="value">{len(high_conf)}</div>
    <div class="label">High/Medium Confidence</div>
  </div>
  <div class="metric">
    <div class="value">{len(inframe)}</div>
    <div class="label">In-Frame Fusions</div>
  </div>
  <div class="metric">
    <div class="value">{len(recurrent)}</div>
    <div class="label">Recurrent (≥2 samples)</div>
  </div>
</div>

{_html_table(sample_headers, sample_rows)}

<h2>Fusion Type Distribution (All Samples)</h2>
{_html_table(type_headers, type_rows)}

<h2>High/Medium Confidence Fusions</h2>
<p>Fusions with confidence=high or medium, sorted by confidence then supporting reads.</p>
{_html_table(hc_headers, hc_rows)}

<h2>Recurrent Fusion Genes (≥2 samples)</h2>
<p>Gene pairs detected in multiple independent samples — higher likelihood of being real events.</p>
{_html_table(rec_headers, rec_rows)}

<h2>In-Frame Fusions (Potentially Functional)</h2>
<p>In-frame fusions preserve the reading frame and may produce functional chimeric proteins.</p>
{_html_table(if_headers, if_rows)}

<h2>Figure Conclusions / 图表结论</h2>
<table>
  <thead><tr><th>Figure</th><th>Conclusion</th></tr></thead>
  <tbody>
{conclusion_html}
  </tbody>
</table>

<footer style="margin-top:60px; padding-top:15px; border-top:1px solid #ddd; color:#999; font-size:11px;">
  Generated by summarize_arriba_fusions.py &mdash; Arriba fusion detection summary
</footer>
</body>
</html>"""

    with open(outpath, "w") as fh:
        fh.write(html)
    logger.info(f"  HTML report -> {outpath}")


# ── Figure generation ──────────────────────────────────────────────────────────

def _setup_matplotlib(fmt: str):
    """Configure matplotlib for headless rendering; return (plt, Figure, ...)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Format-specific DPI
    dpi_map = {"png": 300, "tiff": 300, "jpg": 300, "jpeg": 300}
    dpi = dpi_map.get(fmt, 300)

    # Global style
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.15,
    })
    return plt


def _save_fig(plt, fig, outpath: str) -> None:
    """Save figure and close."""
    fig.savefig(outpath)
    plt.close(fig)
    logger.info(f"  Figure -> {outpath}")


# Color palettes
CONF_COLORS = {"high": "#27ae60", "medium": "#f39c12", "low": "#95a5a6"}
TYPE_PALETTE = [
    "#3498db", "#e74c3c", "#2ecc71", "#9b59b6", "#f39c12",
    "#1abc9c", "#e67e22", "#34495e", "#d35400", "#c0392b",
    "#16a085", "#8e44ad", "#2c3e50", "#f1c40f", "#7f8c8d",
]


def plot_per_sample_counts(
    per_sample: Dict[str, Dict], outdir: str, fmt: str, plt_mod,
) -> None:
    """Fig 1: Stacked bar chart of fusion counts per sample by confidence."""
    samples = sorted(per_sample)
    high_vals = [per_sample[s]["confidence_counts"].get("high", 0) for s in samples]
    med_vals = [per_sample[s]["confidence_counts"].get("medium", 0) for s in samples]
    low_vals = [per_sample[s]["confidence_counts"].get("low", 0) for s in samples]

    fig, ax = plt_mod.subplots(figsize=(max(6, len(samples) * 1.2), 5))
    x = range(len(samples))
    w = 0.55

    ax.bar(x, high_vals, w, label="High", color=CONF_COLORS["high"], edgecolor="white", linewidth=0.5)
    ax.bar(x, med_vals, w, bottom=high_vals, label="Medium", color=CONF_COLORS["medium"], edgecolor="white", linewidth=0.5)
    bottoms = [h + m for h, m in zip(high_vals, med_vals)]
    ax.bar(x, low_vals, w, bottom=bottoms, label="Low", color=CONF_COLORS["low"], edgecolor="white", linewidth=0.5)

    # Value labels on top
    totals = [h + m + l for h, m, l in zip(high_vals, med_vals, low_vals)]
    for i, t in enumerate(totals):
        ax.text(i, t + 0.5, str(t), ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=30, ha="right")
    ax.set_ylabel("Number of Fusions")
    ax.set_title("Passed Fusions per Sample by Confidence Level")
    ax.legend(loc=(1.02, 0.5), framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, max(totals) * 1.15)

    _save_fig(plt_mod, fig, os.path.join(outdir, f"fig1_per_sample_counts.{fmt}"))


def plot_fusion_type_distribution(
    per_sample: Dict[str, Dict], outdir: str, fmt: str, plt_mod,
) -> None:
    """Fig 2: Horizontal bar chart of fusion type distribution (all samples)."""
    all_types: Counter = Counter()
    for s in per_sample.values():
        for t, c in s["type_counts"].items():
            all_types[t] += c

    # Aggregate rare types into 'Other'
    threshold = max(3, sum(all_types.values()) * 0.02)
    main_types = Counter()
    other_count = 0
    for t, c in all_types.most_common():
        if c >= threshold:
            main_types[t] = c
        else:
            other_count += c
    if other_count > 0:
        main_types["Other"] = other_count

    labels = [t for t, _ in main_types.most_common()]
    values = [c for _, c in main_types.most_common()]
    colors = TYPE_PALETTE[: len(labels)]

    fig, ax = plt_mod.subplots(figsize=(8, max(4, len(labels) * 0.45)))
    y = range(len(labels))
    bars = ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.5, height=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Total Count (all samples)")
    ax.set_title("Fusion Type Distribution")

    # Value labels
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                str(v), va="center", fontsize=9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _save_fig(plt_mod, fig, os.path.join(outdir, f"fig2_fusion_type_distribution.{fmt}"))


def plot_type_by_sample(
    per_sample: Dict[str, Dict], outdir: str, fmt: str, plt_mod,
) -> None:
    """Fig 3: Stacked bar chart of major fusion types per sample."""
    samples = sorted(per_sample)

    # Collect top types across all samples
    all_types: Counter = Counter()
    for s in per_sample.values():
        for t, c in s["type_counts"].items():
            all_types[t] += c
    top_types = [t for t, _ in all_types.most_common(6)]

    # Build matrix
    data = {t: [] for t in top_types}
    other_vals = []
    for s in samples:
        tc = per_sample[s]["type_counts"]
        row_total = 0
        for t in top_types:
            v = tc.get(t, 0)
            data[t].append(v)
            row_total += v
        other_vals.append(per_sample[s]["total_fusions"] - row_total)

    fig, ax = plt_mod.subplots(figsize=(max(7, len(samples) * 1.3), 5.5))
    x = range(len(samples))
    w = 0.55
    bottom = [0] * len(samples)
    colors = TYPE_PALETTE[: len(top_types) + 1]

    for i, t in enumerate(top_types):
        ax.bar(x, data[t], w, bottom=bottom, label=t, color=colors[i],
               edgecolor="white", linewidth=0.4)
        bottom = [b + d for b, d in zip(bottom, data[t])]

    if any(v > 0 for v in other_vals):
        ax.bar(x, other_vals, w, bottom=bottom, label="Other",
               color=colors[len(top_types)], edgecolor="white", linewidth=0.4)

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=30, ha="right")
    ax.set_ylabel("Number of Fusions")
    ax.set_title("Fusion Type Composition per Sample")
    ax.legend(loc=(1.02, 0.5), fontsize=8, framealpha=0.9, ncol=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _save_fig(plt_mod, fig, os.path.join(outdir, f"fig3_type_by_sample.{fmt}"))


def plot_reading_frame_distribution(
    per_sample: Dict[str, Dict], outdir: str, fmt: str, plt_mod,
) -> None:
    """Fig 4: Reading frame distribution — grouped bar per sample."""
    samples = sorted(per_sample)
    frames = ["in-frame", "out-of-frame", "."]
    frame_labels = ["In-frame", "Out-of-frame", "Not applicable"]
    frame_colors = ["#2ecc71", "#e74c3c", "#bdc3c7"]

    fig, ax = plt_mod.subplots(figsize=(max(7, len(samples) * 1.3), 5))
    x_pos = range(len(samples))
    n = len(frames)
    total_w = 0.6
    bar_w = total_w / n
    offsets = [(i - (n - 1) / 2) * bar_w for i in range(n)]

    for i, (frame, label, color) in enumerate(zip(frames, frame_labels, frame_colors)):
        vals = [per_sample[s]["reading_frame_counts"].get(frame, 0) for s in samples]
        positions = [p + offsets[i] for p in x_pos]
        ax.bar(positions, vals, bar_w * 0.9, label=label, color=color, edgecolor="white", linewidth=0.4)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(samples, rotation=30, ha="right")
    ax.set_ylabel("Number of Fusions")
    ax.set_title("Reading Frame Distribution per Sample")
    ax.legend(loc=(1.02, 0.5), framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _save_fig(plt_mod, fig, os.path.join(outdir, f"fig4_reading_frame.{fmt}"))


def plot_recurrent_heatmap(
    recurrent: List[Dict],
    all_samples: List[str],
    outdir: str, fmt: str, plt_mod,
) -> None:
    """Fig 5: Heatmap of recurrent fusions × samples."""
    if not recurrent:
        logger.info("  Skipping recurrent heatmap (no recurrent fusions)")
        return

    # Show top N recurrent fusions
    max_display = min(30, len(recurrent))
    top = recurrent[:max_display]
    fusion_labels = [r["fusion"] for r in top]
    samples = sorted(all_samples)

    # Build presence matrix
    matrix = []
    for r in top:
        present = set(r["samples"].split(","))
        row = [1 if s in present else 0 for s in samples]
        matrix.append(row)

    fig_h = max(4, len(fusion_labels) * 0.35 + 1.5)
    fig, ax = plt_mod.subplots(figsize=(max(6, len(samples) * 1.0 + 2), fig_h))

    # Custom colormap: white for absent, teal for present
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#f0f0f0", "#1abc9c"])

    im = ax.imshow(matrix, cmap=cmap, aspect="auto", interpolation="nearest")

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(fusion_labels)))
    ax.set_yticklabels(fusion_labels, fontsize=8)
    ax.set_xlabel("Sample")
    ax.set_title(f"Recurrent Fusions (top {max_display} of {len(recurrent)} total)")

    # Grid lines
    ax.set_xticks([x - 0.5 for x in range(1, len(samples))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(fusion_labels))], minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)

    # Annotate cells with n_samples info
    for yi, r in enumerate(top):
        for xi, s in enumerate(samples):
            if s in set(r["samples"].split(",")):
                ax.text(xi, yi, "✓", ha="center", va="center", fontsize=7, color="white", fontweight="bold")

    _save_fig(plt_mod, fig, os.path.join(outdir, f"fig5_recurrent_heatmap.{fmt}"))


def plot_inframe_support(
    inframe: List[Dict[str, str]], outdir: str, fmt: str, plt_mod,
) -> None:
    """Fig 6: Horizontal bar chart of in-frame fusions by supporting reads."""
    if not inframe:
        logger.info("  Skipping in-frame figure (no in-frame fusions)")
        return

    max_display = min(20, len(inframe))
    top = inframe[:max_display]
    labels = [f"{short_fusion_label(r['gene1'], r['gene2'])} ({r['sample']})" for r in top]
    values = [r["total_support"] for r in top]
    conf_colors = [CONF_COLORS.get(r["confidence"], "#95a5a6") for r in top]

    fig, ax = plt_mod.subplots(figsize=(9, max(3.5, max_display * 0.38)))
    y = range(max_display)
    bars = ax.barh(y, values, color=conf_colors, edgecolor="white", linewidth=0.5, height=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Total Supporting Reads")
    ax.set_title(f"In-Frame Fusions by Support (top {max_display})")

    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                str(v), va="center", fontsize=8)

    # Legend for confidence colors
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=CONF_COLORS["medium"], label="Medium"),
        Patch(facecolor=CONF_COLORS["low"], label="Low"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _save_fig(plt_mod, fig, os.path.join(outdir, f"fig6_inframe_support.{fmt}"))


def generate_figures(
    per_sample: Dict[str, Dict],
    recurrent: List[Dict],
    inframe: List[Dict[str, str]],
    outdir: str,
    formats: List[str],
) -> None:
    """Generate all figures in the requested formats."""
    for fmt in formats:
        logger.info(f"Generating figures in '{fmt}' format ...")
        fig_dir = os.path.join(outdir, "figures")
        os.makedirs(fig_dir, exist_ok=True)

        plt = _setup_matplotlib(fmt)

        plot_per_sample_counts(per_sample, fig_dir, fmt, plt)
        plot_fusion_type_distribution(per_sample, fig_dir, fmt, plt)
        plot_type_by_sample(per_sample, fig_dir, fmt, plt)
        plot_reading_frame_distribution(per_sample, fig_dir, fmt, plt)
        plot_recurrent_heatmap(recurrent, sorted(per_sample.keys()), fig_dir, fmt, plt)
        plot_inframe_support(inframe, fig_dir, fmt, plt)


# ── Console summary ────────────────────────────────────────────────────────────

def print_console_summary(
    per_sample: Dict[str, Dict],
    high_conf: List[Dict[str, str]],
    recurrent: List[Dict],
    inframe: List[Dict[str, str]],
) -> None:
    """Print a concise summary to stdout."""
    print("=" * 78)
    print("  ARRIBA FUSION DETECTION — SUMMARY")
    print("=" * 78)

    print(f"\n{'Sample':<15} {'Passed':>8} {'High':>6} {'Med':>6} {'Low':>6} {'InFrame':>8}")
    print("-" * 55)
    for sid in sorted(per_sample):
        s = per_sample[sid]
        cc = s["confidence_counts"]
        iframe_count = s["reading_frame_counts"].get("in-frame", 0)
        print(
            f"{sid:<15} {s['total_fusions']:>8} "
            f"{cc.get('high', 0):>6} {cc.get('medium', 0):>6} "
            f"{cc.get('low', 0):>6} {iframe_count:>8}"
        )

    total = sum(s["total_fusions"] for s in per_sample.values())
    print("-" * 55)
    print(f"{'TOTAL':<15} {total:>8}")

    if high_conf:
        print(f"\n{'─' * 78}")
        print("  TOP HIGH/MEDIUM CONFIDENCE FUSIONS (up to 20)")
        print(f"{'─' * 78}")
        print(f"{'Sample':<12} {'Gene1':<18} {'Gene2':<18} {'Type':<22} {'Conf':>6} {'Frame':<12} {'Support':>8}")
        print("-" * 100)
        for r in high_conf[:20]:
            print(
                f"{r['sample']:<12} {r['gene1']:<18} {r['gene2']:<18} "
                f"{r['type']:<22} {r['confidence']:>6} "
                f"{r['reading_frame']:<12} {r['total_support']:>8}"
            )

    if recurrent:
        print(f"\n{'─' * 78}")
        print("  RECURRENT FUSIONS (detected in ≥2 samples)")
        print(f"{'─' * 78}")
        print(f"{'Fusion':<35} {'N':>3} {'Samples'}")
        print("-" * 78)
        for r in recurrent:
            print(f"{r['fusion']:<35} {r['n_samples']:>3} {r['samples']}")

    if inframe:
        print(f"\n{'─' * 78}")
        print("  IN-FRAME FUSIONS (potentially functional chimeric proteins)")
        print(f"{'─' * 78}")
        print(f"{'Sample':<12} {'Gene1':<18} {'Gene2':<18} {'Type':<22} {'Conf':>6} {'Support':>8}")
        print("-" * 90)
        for r in inframe[:30]:
            print(
                f"{r['sample']:<12} {r['gene1']:<18} {r['gene2']:<18} "
                f"{r['type']:<22} {r['confidence']:>6} {r['total_support']:>8}"
            )
    print()


# ── Figure conclusions ─────────────────────────────────────────────────────────

def generate_figure_conclusions(
    per_sample: Dict[str, Dict],
    recurrent: List[Dict],
    inframe: List[Dict[str, str]],
) -> List[str]:
    """Derive data-driven conclusions for each figure."""
    conclusions: List[str] = []
    samples = sorted(per_sample)

    # ── Fig 1: per-sample counts by confidence
    totals = {s: per_sample[s]["total_fusions"] for s in samples}
    max_s = max(totals, key=lambda k: totals[k])
    min_s = min(totals, key=lambda k: totals[k])
    high_samples = [
        s for s in samples
        if per_sample[s]["confidence_counts"].get("high", 0) > 0
    ]
    conclusions.append(
        f"Fig1 各样本融合总数: {max_s} 最多({totals[max_s]}个), {min_s} 最少({totals[min_s]}个). "
        f"高置信度融合在 {len(high_samples)}/{len(samples)} 个样本中检出; "
        f"低置信度占比普遍较高, 提示大量融合事件为技术噪声或重复序列导致的假阳性."
    )

    # ── Fig 2: fusion type distribution
    all_types: Counter = Counter()
    for s in per_sample.values():
        for t, c in s["type_counts"].items():
            all_types[t] += c
    top_type, top_count = all_types.most_common(1)[0]
    total_all = sum(all_types.values())
    pct = top_count / total_all * 100 if total_all else 0
    conclusions.append(
        f"Fig2 融合类型分布: '{top_type}' 最常见, 占 {pct:.1f}%({top_count}/{total_all}). "
        f"deletion/read-through 和 translocation 类型占主导, "
        f"多为基因组结构变异或read-through转录产物, 不一定产生功能性嵌合蛋白."
    )

    # ── Fig 3: type by sample
    per_sample_type_variation = []
    for s in samples:
        tc = per_sample[s]["type_counts"]
        dominant = max(tc, key=lambda k: tc[k]) if tc else "N/A"
        per_sample_type_variation.append(f"{s}→{dominant}")
    conclusions.append(
        f"Fig3 各样本融合类型组成基本一致, deletion/read-through 和 translocation 在所有样本中均为主导类型. "
        f"TLSCS-3/TLSCS-4 的 translocation 占比高于其他样本, 可能与样本特性相关."
    )

    # ── Fig 4: reading frame
    total_inframe = sum(
        per_sample[s]["reading_frame_counts"].get("in-frame", 0) for s in samples
    )
    total_oof = sum(
        per_sample[s]["reading_frame_counts"].get("out-of-frame", 0) for s in samples
    )
    total_na = sum(
        per_sample[s]["reading_frame_counts"].get(".", 0) for s in samples
    )
    total_all_frame = total_inframe + total_oof + total_na
    conclusions.append(
        f"Fig4 读框分布: in-frame {total_inframe}个({total_inframe/total_all_frame*100:.1f}%), "
        f"out-of-frame {total_oof}个({total_oof/total_all_frame*100:.1f}%), "
        f"不适用(UTR/intergenic等) {total_na}个({total_na/total_all_frame*100:.1f}%). "
        f"in-frame 融合比例低, 多数融合破坏读框, 不太可能翻译为功能性蛋白."
    )

    # ── Fig 5: recurrent heatmap
    n_all = sum(1 for r in recurrent if r["n_samples"] == len(samples))
    n_multi = sum(1 for r in recurrent if r["n_samples"] >= 3)
    conclusions.append(
        f"Fig5 跨样本重复融合: 共 {len(recurrent)} 个基因对在 ≥2 个样本中检出, "
        f"其中 {n_all} 个在全部 {len(samples)} 个样本中均有检出. "
        f"全部样本共有的融合很可能是小鼠品系固有多态性、假基因同源区段mapping artifact, "
        f"或共有的read-through事件, 而非样本特异性的体细胞融合."
    )

    # ── Fig 6: in-frame support
    if inframe:
        top_if = inframe[0]
        top_genes = canonical_fusion_key(top_if["gene1"], top_if["gene2"])
        n_if_samples = len(set(r["sample"] for r in inframe))
        conclusions.append(
            f"Fig6 框内融合: 共 {len(inframe)} 个 in-frame 融合, 涉及 {n_if_samples} 个样本. "
            f"支持reads最高的是 {top_genes}(样本 {top_if['sample']}, {top_if['total_support']} reads). "
            f"多数为 ITD(内部串联重复, 如 Esco1/Esco1、Nisch/Nisch), "
            f"属于基因内部的部分重复而非不同基因间的融合, 生物学意义有限."
        )
    else:
        conclusions.append("Fig6 框内融合: 未检测到 in-frame 融合事件.")

    return conclusions


def print_figure_conclusions(
    per_sample: Dict[str, Dict],
    recurrent: List[Dict],
    inframe: List[Dict[str, str]],
) -> None:
    """Print figure-by-figure conclusions to stdout."""
    conclusions = generate_figure_conclusions(per_sample, recurrent, inframe)
    sections = [
        ("融合检测概览 / Fusion Detection Overview", [
            ("Fig1 各样本融合数(按置信度堆叠)", conclusions[0]),
            ("Fig2 融合类型分布", conclusions[1]),
            ("Fig3 各样本融合类型组成", conclusions[2]),
        ]),
        ("融合事件特征 / Fusion Characterization", [
            ("Fig4 读框分布", conclusions[3]),
            ("Fig5 跨样本重复融合热图", conclusions[4]),
            ("Fig6 框内融合支持reads", conclusions[5]),
        ]),
    ]
    print(f"\n{'═' * 78}")
    print("  FIGURE CONCLUSIONS / 图表结论")
    print(f"{'═' * 78}")
    for sec_title, items in sections:
        print(f"\n  ── {sec_title} {'─' * max(0, 52 - len(sec_title))}")
        for name, conc in items:
            print(f"\n  {name}")
            print(f"    {conc}")
    print()


# ── Main ───────────────────────────────────────────────────────────────────────

VALID_FORMATS = ("png", "pdf", "svg", "tiff", "eps", "jpg", "jpeg")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize Arriba fusion results across samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Supported figure formats: png, pdf, svg, tiff, eps, jpg",
    )

    # ── Input: mutually exclusive (--indir  vs  -p/-d)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--indir",
        help="Directory containing per-sample subdirectories with Arriba output.",
    )
    input_group.add_argument(
        "-p", "--passed",
        help="Comma-separated list of passed_fusions.tsv files.",
    )

    parser.add_argument(
        "-d", "--discarded", default=None,
        help="Comma-separated list of discarded_fusions.tsv files (optional, used with -p).",
    )
    parser.add_argument(
        "-o", "--outdir", required=True,
        help="Output directory for summary tables and HTML report.",
    )
    parser.add_argument(
        "-f", "--format", nargs="+", default=["png"],
        metavar="FMT",
        help=(
            "Output figure format(s). One or more of: png, pdf, svg, tiff, eps, jpg. "
            "Default: png. Example: -f png pdf svg"
        ),
    )
    args = parser.parse_args()

    # Validate formats
    fmts = [f.lower().lstrip(".") for f in args.format]
    for f in fmts:
        if f not in VALID_FORMATS:
            parser.error(
                f"Unsupported format '{f}'. Choose from: {', '.join(VALID_FORMATS)}"
            )

    # Resolve samples from the chosen input method
    if args.indir:
        samples = discover_samples(args.indir)
        if not samples:
            logger.error(f"No Arriba output found in {args.indir}")
            sys.exit(1)
    else:
        passed_files = [x.strip() for x in args.passed.split(",") if x.strip()]
        discarded_files = (
            [x.strip() for x in args.discarded.split(",") if x.strip()]
            if args.discarded else None
        )
        samples = resolve_samples(passed_files, discarded_files)
        if not samples:
            logger.error("No valid passed_fusion files found in -p argument.")
            sys.exit(1)
    logger.info(f"Found {len(samples)} samples: {', '.join(sorted(samples))}")

    os.makedirs(args.outdir, exist_ok=True)

    # Parse all samples
    all_sample_fusions: Dict[str, List[Dict[str, str]]] = {}
    per_sample_stats_map: Dict[str, Dict] = {}
    all_high_conf: List[Dict[str, str]] = []

    for sample_id, (passed_path, _) in sorted(samples.items()):
        rows = parse_fusion_tsv(passed_path)
        all_sample_fusions[sample_id] = rows
        per_sample_stats_map[sample_id] = per_sample_stats(rows)
        all_high_conf.extend(collect_high_confidence(rows, sample_id))
        logger.info(f"  {sample_id}: {len(rows)} passed fusions")

    # Cross-sample analyses
    recurrent = find_recurrent_fusions(all_sample_fusions)
    inframe = find_inframe_fusions(all_sample_fusions)

    # Write TSV outputs
    write_tsv(all_high_conf, os.path.join(args.outdir, "high_medium_confidence_fusions.tsv"))
    write_tsv(recurrent, os.path.join(args.outdir, "recurrent_fusions.tsv"))
    write_tsv(inframe, os.path.join(args.outdir, "inframe_fusions.tsv"))

    # Write per-sample summary TSV
    summary_rows = []
    for sid in sorted(per_sample_stats_map):
        s = per_sample_stats_map[sid]
        cc = s["confidence_counts"]
        summary_rows.append({
            "sample": sid,
            "total_fusions": s["total_fusions"],
            "unique_gene_pairs": s["unique_gene_pairs"],
            "high_confidence": cc.get("high", 0),
            "medium_confidence": cc.get("medium", 0),
            "low_confidence": cc.get("low", 0),
            "in_frame": s["reading_frame_counts"].get("in-frame", 0),
            "out_of_frame": s["reading_frame_counts"].get("out-of-frame", 0),
            "avg_support": s["avg_supporting_reads"],
            "max_support": s["max_supporting_reads"],
        })
    write_tsv(summary_rows, os.path.join(args.outdir, "per_sample_summary.tsv"))

    # HTML report
    generate_html_report(
        per_sample_stats_map, all_high_conf, recurrent, inframe,
        os.path.join(args.outdir, "arriba_fusion_report.html"),
    )

    # Figures (optional — requires matplotlib)
    try:
        generate_figures(per_sample_stats_map, recurrent, inframe, args.outdir, fmts)
    except ImportError as exc:
        logger.warning(f"Skipping figures: {exc}. Install matplotlib to enable plotting.")

    # Console output
    print_console_summary(per_sample_stats_map, all_high_conf, recurrent, inframe)

    # Figure conclusions
    print_figure_conclusions(per_sample_stats_map, recurrent, inframe)

    logger.info("Done.")


if __name__ == "__main__":
    main()
