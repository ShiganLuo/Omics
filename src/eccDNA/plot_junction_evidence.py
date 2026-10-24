#!/usr/bin/env python3
"""Plot method-2 junction evidence across the four ecDNA samples.

Figure: two panels.
  A. Stacked bars of nonsplice / backsplice junction reads per ecDNA circle,
     grouped by sample (both replicates side by side), log-scaled internal-read
     background as a dot row (transcription-permissive coverage control).
  B. Ring-seam soft-clip verified read counts per sample (the direct
     ring-closure evidence), annotated with the 3-read threshold.

The data come from the SE subset realignment (see align_subset_se.sh) and the
independent soft-clip seam scan; counts are hard-coded here from
junction_results/*.se.tsv + the seam scan output recorded in RESULTS.md.
"""
from __future__ import annotations

import argparse
import logging
from typing import List, Literal, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("plot_junction_evidence")

PlotFormat = Literal["png", "pdf", "svg"]

# ---------------------------------------------------------------------------
# Data (from junction_results/*.se.tsv and the soft-clip seam scan)
# ---------------------------------------------------------------------------
# (sample, rep): {circ_id: (nonsplice, backsplice, internal)}
JUNCTIONS = {
    ("DMSO_P06", 1): {"33": (2, 0, 11665), "105": (0, 0, 22723), "104": (0, 0, 28),
                       "24": (0, 0, 14), "47": (0, 0, 33012), "66": (0, 0, 6), "9": (0, 0, 2091)},
    ("DMSO_P06", 2): {"33": (8, 0, 11319), "105": (0, 0, 23055), "104": (0, 0, 51),
                       "24": (0, 0, 8), "47": (0, 0, 33604), "66": (0, 0, 4), "9": (0, 0, 2068)},
    ("DMSO_P20", 1): {"16": (0, 0, 88), "22": (0, 0, 21132), "32": (0, 0, 12003),
                       "55": (0, 0, 15), "81": (0, 0, 39), "82": (0, 0, 21565), "83": (0, 0, 390)},
    ("DMSO_P20", 2): {"16": (0, 0, 92), "22": (0, 0, 20962), "32": (0, 0, 11914),
                       "55": (0, 0, 0), "81": (0, 0, 27), "82": (0, 0, 22801), "83": (0, 0, 406)},
    ("PlaB_P06", 1): {"35": (6, 0, 18897), "100": (0, 0, 16524), "21": (0, 0, 92),
                       "69": (0, 0, 15), "99": (0, 0, 363)},
    ("PlaB_P06", 2): {"35": (5, 1, 19816), "100": (0, 0, 15204), "21": (0, 0, 98),
                       "69": (0, 0, 8), "99": (0, 0, 373)},
    ("PlaB_P20", 1): {"31": (8, 1, 24454), "18": (0, 0, 114), "43": (0, 0, 31485),
                       "58": (0, 0, 14), "81": (0, 0, 663), "82": (0, 0, 15755), "9": (0, 0, 2564)},
    ("PlaB_P20", 2): {"31": (5, 0, 23382), "18": (0, 0, 144), "43": (0, 0, 30614),
                       "58": (0, 0, 28), "81": (0, 0, 715), "82": (0, 0, 15037), "9": (0, 0, 2419)},
}
# Ring-seam (linearization-point) soft-clip verified reads, per (sample, rep):
SEAM = {
    ("DMSO_P06", 1): {"105": 1},
    ("DMSO_P06", 2): {"33": 4, "105": 1},
    ("DMSO_P20", 1): {},
    ("DMSO_P20", 2): {},
    ("PlaB_P06", 1): {"35": 1},
    ("PlaB_P06", 2): {},
    ("PlaB_P20", 1): {"31": 1},
    ("PlaB_P20", 2): {"31": 1},
}
# Junction opportunities per ring: single-fragment = 1 (seam), 2-frag = 2 (1 fj + seam),
# 3-frag = 3 (2 fj + seam). Fragment-junction threshold must be judged per opportunity.
N_JUNCTIONS = {
    "DMSO_P06": {"33": 2, "105": 1, "104": 1, "24": 1, "47": 1, "66": 1, "9": 1},
    "DMSO_P20": {"82": 1, "22": 2, "32": 3, "83": 3, "16": 1, "81": 1, "55": 1},
    "PlaB_P06": {"35": 3, "100": 1, "21": 1, "69": 1, "99": 1},
    "PlaB_P20": {"31": 3, "43": 1, "18": 1, "58": 1, "81": 1, "82": 1, "9": 1},
}
SAMPLES = ["DMSO_P06", "DMSO_P20", "PlaB_P06", "PlaB_P20"]
# The chr4 homologous ring in each sample (circ id differs per decoil run)
HOMOLOG = {"DMSO_P06": "33", "DMSO_P20": "22", "PlaB_P06": "35", "PlaB_P20": "31"}
THRESHOLD = 3


def _aggregate() -> dict:
    """Per (sample, rep): seam / fj nonsplice / backsplice totals for the homologous ring.

    Seam counts come from the soft-clip verified scan (SEAM), homologous ring
    only; fragment-junction and backsplice counts from align_ecdna_junction.py
    (JUNCTIONS). Non-homolog seam hits (e.g. circ105) are not aggregated here.
    """
    out = {}
    for key, data in JUNCTIONS.items():
        sample, rep = key
        hom = HOMOLOG[sample]
        ns, bs, _ = data[hom]
        out[key] = {"seam": SEAM[key].get(hom, 0),
                    "fj": ns,
                    "bs": bs}
    return out


def plot_panel_a(ax: plt.Axes, agg: dict) -> None:
    """Grouped bars per sample/rep: seam vs fragment-junction vs backsplice.

    Geometry: each sample group is centered exactly at x = gi. Within a group,
    rep1 occupies the left half, rep2 the right half; each half holds 3 bars.
    """
    cats = [("seam", "#b2182b", "junction_spanning_nonsplice\n(ring seam, decisive)"),
            ("fj", "#35618f", "junction_spanning_nonsplice\n(fragment junction)"),
            ("bs", "#c9a227", "junction_spanning_backsplice\n(circRNA-suspect control)")]
    group_half = 0.38          # half-width of one sample's group
    bar_w = group_half / 3     # 3 bars per rep-half
    for gi, sample in enumerate(SAMPLES):
        for rep in (1, 2):
            sign = -1 if rep == 1 else 1
            d = agg[(sample, rep)]
            for ci, (cat, color, _) in enumerate(cats):
                x = gi + sign * (group_half - (ci + 0.5) * bar_w)
                v = d[cat]
                ax.bar(x, max(v, 0.04), width=bar_w * 0.88, color=color,
                       edgecolor="white", lw=0.4, zorder=2)
                if v:
                    ax.text(x, v + 0.12, str(v), ha="center", fontsize=6.5, zorder=4)
    # threshold line + label at far left, clear of bars (max bar there is small)
    ax.axhline(THRESHOLD, color="#555555", lw=0.9, ls="--", alpha=0.8, zorder=1)
    ax.annotate(f"expression threshold = {THRESHOLD} seam reads",
                xy=(len(SAMPLES) - 0.45, THRESHOLD), xytext=(6, 0),
                textcoords="offset points", fontsize=6.5, color="#555555",
                ha="left", va="center", annotation_clip=False)
    # aligned multi-row tick labels: all centered exactly at gi
    for gi, sample in enumerate(SAMPLES):
        ax.text(gi - group_half / 2, -0.35, "rep1", ha="center", va="top",
                fontsize=6, color="#999999")
        ax.text(gi + group_half / 2, -0.35, "rep2", ha="center", va="top",
                fontsize=6, color="#999999")
        ax.text(gi, -0.95, f"circ{HOMOLOG[sample]}", ha="center", va="top",
                fontsize=7.5, color="#35618f")
        ax.text(gi, -1.75, sample, ha="center", va="top", fontsize=8.5)
    ax.set_ylabel("Junction-spanning reads (homologous ring)")
    ax.set_ylim(0, 10)
    ax.set_xlim(-0.55, len(SAMPLES) - 0.45)
    ax.set_xticks([])
    ax.spines["bottom"].set_position(("data", 0))
    ax.tick_params(axis="x", length=0)
    ax.set_title("A  Method-2 junction evidence per ecDNA circle", fontsize=10,
                 loc="left", pad=6)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in cats]
    labels = [lab.replace("\n", " ") for _, _, lab in cats]
    ax.legend(handles, labels, fontsize=6.8, frameon=False,
              loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.14))


def plot_panel_b(ax: plt.Axes, agg: dict) -> None:
    """Seam nonsplice (the decisive count) as the headline panel."""
    bar_w = 0.30
    gap = 0.06                 # explicit gap between rep1 and rep2 bars
    for gi, sample in enumerate(SAMPLES):
        for rep in (1, 2):
            n = agg[(sample, rep)]["seam"]
            sign = -1 if rep == 1 else 1
            x = gi + sign * (bar_w / 2 + gap / 2)
            hit = n >= THRESHOLD
            ax.bar(x, max(n, 0.05), width=bar_w,
                   color="#b2182b" if hit else "#92a5bd",
                   edgecolor="white", lw=0.5, zorder=2)
            if n:
                ax.text(x, n + 0.1, str(n), ha="center", fontsize=8,
                        fontweight="bold" if hit else "normal", zorder=3)
            else:
                ax.text(x, 0.12, "0", ha="center", fontsize=6.5, color="#aaaaaa")
    ax.axhline(THRESHOLD, color="#b2182b", lw=0.9, ls="--", alpha=0.8)
    ax.text(len(SAMPLES) - 0.5, THRESHOLD + 0.12, f"threshold = {THRESHOLD}",
            fontsize=7, color="#b2182b", ha="right", va="bottom")
    ax.set_xticks(range(len(SAMPLES)))
    ax.set_xticklabels(SAMPLES, fontsize=8.5)
    ax.set_ylabel("Ring-seam nonsplice reads\n(decisive evidence)")
    ax.set_ylim(0, 5.2)
    ax.set_xlim(-0.55, len(SAMPLES) - 0.45)
    ax.set_title("B  Ring-seam spanning reads (decisive)", fontsize=10,
                 loc="left", pad=6)


def run(out_prefix: str, formats: Optional[List[PlotFormat]] = None) -> None:
    """Render the two-panel junction-evidence figure.

    Parameters
    ----------
    out_prefix : str
        Output path prefix; each format appends its extension.
    formats : list of str, optional
        Image formats. Default ``["png"]``.
    """
    if formats is None:
        formats = ["png"]
    agg = _aggregate()
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [2.1, 1]},
        constrained_layout=True)
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.88))
    plot_panel_a(ax_a, agg)
    plot_panel_b(ax_b, agg)
    for fmt in formats:
        path = f"{out_prefix}.{fmt}"
        fig.savefig(path, dpi=300)
        logger.info("saved %s", path)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("-o", "--out-prefix", default="junction_evidence",
                        help="Output path prefix without extension.")
    parser.add_argument("-f", "--format", action="append", dest="formats",
                        metavar="FMT",
                        help="Image format, repeatable. Default: png.")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    run(args.out_prefix, args.formats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
