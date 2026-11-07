#!/usr/bin/env python
"""Visualization suite for mei_analysis.py results.

Reads the mei/ outputs (mei_events.tsv etc.) and renders the publication-
style figure set:

  mei_family_by_group.png    family x group stacked bar + group totals
  mei_insert_length_by_family.png      insertion length distribution per family
  mei_genomic_context.png            genomic context (exon/intron/intergenic)
                                  per family x group (100% stacked)
  mei_te_region_ratio.png          TE-into-TE vs unique-site insertion ratio
  mei_full_length_L1.png     full-length L1 candidates: length hist +
                                  per-group counts
  mei_polya_tsd_fingerprint.png          MEI fingerprint rates (polyA, TSD) per family
  mei_report.html                 self-contained summary report embedding
                                  all figures and key tables

Usage:
  python mei_plot.py --mei-dir <mei_outdir> [--title "youth vs aged"]

Layout: every figure is composed first, then de-collided by the layout
engine at the bottom of this file (_stack_labels + _fix_overlaps): all
text is measured with the real renderer (pixel-accurate boxes) and
nudged, stacked or hidden until no two boxes intersect. Figures therefore
render collision-free for ANY dataset, label length, tick count or group
names — no per-figure hand tuning.
"""

import argparse
import base64
import json
import logging
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.text
import numpy as np
import pandas as pd

logger = logging.getLogger("mei_plot")

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# default group vocabulary; overridden in main() from --ctrl-name/--exp-name
# so custom designs (serum_vs_KSR) get readable labels instead of ctrl/exp.
GROUP_ORDER = ["ctrl_specific", "exp_specific", "shared"]
GROUP_BASE = None  # (ctrl, exp) set at runtime
GROUP_LABEL = {"ctrl_specific": "ctrl-specific",
               "exp_specific": "exp-specific",
               "shared": "shared"}
GROUP_COLOR = {"ctrl_specific": "#4A72B7", "exp_specific": "#D73027",
               "shared": "#7F7F7F"}

# Fallback bucket vocabulary for mei dirs produced before te_buckets was
# recorded in mei_summary.json (plot-only; analysis always writes its own).
FAM_ORDER = ["Alu", "L1", "SVA", "LTR", "DNA", "microsatellite", "other",
             "unclassified", "short_unclassified"]
FAM_COLOR = {"Alu": "#E08214", "L1": "#4A72B7", "SVA": "#8E44AD",
             "LTR": "#27AE60", "DNA": "#16A085", "microsatellite": "#95A5A6",
             "B1": "#F0A04B", "B2/ID": "#E4B363", "IAP/ETn": "#C0392B",
             "other": "#BDC3C7", "unclassified": "#7F8C8D",
             "short_unclassified": "#D5D8DC"}
FAM_DISPLAY = {"microsatellite": "microsatellite"}

# Real group names from the workflow config (Params.sv.comparisons.<key>),
# injected via --ctrl-name/--exp-name. Falls back to generic labels.
GROUP_NAMES = {"ctrl": "ctrl", "exp": "exp"}


def group_label(key: str) -> str:
    """Display label for a group key: '<name>-specific' or 'shared'.

    Name-agnostic: works both on the default ctrl_/exp_ vocabulary and on
    relabeled keys produced by --ctrl-name/--exp-name (main() renames the
    group column before any figure is drawn).
    """
    if key == "shared":
        return "shared"
    if key.endswith("_specific"):
        return key[:-len("_specific")].replace("_", " ") + "-specific"
    return key


def _clean_str(s):
    """NaN/None/'nan' -> '' so empty-string semantics survive a CSV round-trip."""
    out = pd.Series(s.astype("object").where(s.notna(), ""))
    return out.apply(lambda v: "" if v is None or str(v) == "nan" else str(v))


def b64(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def bucket_vocab(summary, df):
    """(slot_buckets, fingerprint_buckets) from the recorded te_buckets
    profile; falls back to FAM_ORDER for pre-profile mei dirs.

    slot_buckets: ordered figure slots — every profile bucket with
    zero_slot != false (zero-count buckets KEEP their slot, labelled 0),
    plus short_unclassified when present in the data.
    fingerprint_buckets: buckets flagged retro_fingerprint (fig6).
    """
    prof = (summary or {}).get("te_buckets")
    if isinstance(prof, list) and prof:
        names = [b.get("name") for b in prof
                 if isinstance(b, dict) and b.get("name")]
        slots = [n for b, n in zip(prof, names) if b.get("zero_slot", True)]
        fp = [n for b, n in zip(prof, names) if b.get("retro_fingerprint")]
    else:
        present = set(df["family"].dropna())
        slots = [f for f in FAM_ORDER if f in present]
        fp = [f for f in slots if f not in
              ("microsatellite", "other", "unclassified",
               "short_unclassified")]
    if "short_unclassified" in set(df["family"].dropna()) \
            and "short_unclassified" not in slots:
        slots.append("short_unclassified")
    return slots, fp


def fam_color(f):
    return FAM_COLOR.get(f, "#BDC3C7")


def _empty_state(fig, out, msg="no MEI events"):
    """Degenerate-data figure: explicit message, declared output exists."""
    for ax in fig.axes:
        ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=11,
                color="#7F7C8D", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    _save(fig, out)


def fig1_family_by_group(df, out, fams=None):
    fams = fams or FAM_ORDER
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4),
                             gridspec_kw={"width_ratios": [2, 1]})
    groups = [g for g in GROUP_ORDER if (df["group"] == g).any()]
    if not groups or df.empty:
        _empty_state(fig, out)
        return

    # left: grouped bars — every profile bucket keeps its slot; EVERY bar
    # is labelled with its count (zero-count groups get a gray "0" at the
    # axis floor, so zero is explicit, not ambiguous)
    ax = axes[0]
    x = np.arange(len(fams))
    w = 0.8 / max(len(groups), 1)
    vtexts = []
    for gi, g in enumerate(groups):
        vals = [(df[(df.family == f) & (df.group == g)]).shape[0]
                for f in fams]
        pos = x + gi * w - 0.4 + w / 2
        ax.bar(pos, vals, width=w * 0.92,
               color=GROUP_COLOR[g], label=group_label(g))
        for xi, v in zip(pos, vals):
            if v:
                vtexts.append(ax.text(xi, v, f"{v:,}", ha="center",
                                      va="bottom", fontsize=7,
                                      color="#333333"))
            else:
                vtexts.append(ax.text(xi, 0.56, "0", ha="center",
                                      va="center", fontsize=7,
                                      color="#95A5A6"))
    ax.set_xticks(x)
    labels = [FAM_DISPLAY.get(f, f) for f in fams]
    rot, ha = fit_ticks(ax, labels)
    ax.set_xticklabels(labels, rotation=rot, ha=ha,
                       rotation_mode="anchor", fontsize=8.5)
    ax.set_ylabel("insertions")
    ax.set_yscale("log")
    ax.set_ylim(0.5, None)

    # right: totals per group (same family subset as panel A)
    sub = df[df.family.isin(fams)]
    tot = [int((sub.group == g).sum()) for g in groups]
    axes[1].bar(range(len(groups)), tot,
                color=[GROUP_COLOR[g] for g in groups])
    for i, v in enumerate(tot):
        axes[1].text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=9)
    axes[1].set_xticks(range(len(groups)))
    glabels = [group_label(g) for g in groups]
    grot, gha = fit_ticks(axes[1], glabels)
    axes[1].set_xticklabels(glabels, rotation=grot, ha=gha,
                            rotation_mode="anchor", fontsize=8.5)
    axes[1].set_ylabel("total insertions")
    axes[1].set_title("B  group totals")
    if tot:
        axes[1].set_ylim(0, max(tot) * 1.18)

    # geometry is final NOW: fit the legend band against the real axes
    # width, then the title pad-measured above it; a second tight_layout
    # reserves room for the padded title so it cannot clip at the figure
    # edge; value labels are de-collided only after that
    fig.tight_layout()
    leg, leg_h = legend_above(ax, ncol=max(len(groups), 1))
    ax.set_title("A  TE family x group (log scale)",
                 pad=title_pad(leg_h))
    fig.tight_layout()
    stack_value_labels(fig, vtexts)
    raise_ylim_for_labels(fig, ax, vtexts)
    _save(fig, out)


def fig2_length_family(df, out, fams=None):
    fams = fams or ["Alu", "L1", "SVA", "LTR"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.logspace(1.7, 4.3, 45)
    for f in fams:
        vals = df.loc[df.family == f, "ins_len"].clip(lower=50)
        if len(vals) == 0:
            continue
        ax.hist(vals, bins=bins, alpha=0.55, label=f"{f} (n={len(vals):,})",
                color=fam_color(f))
    ax.set_xscale("log")
    ax.set_xlabel("inserted sequence length (bp, log)")
    ax.set_ylabel("count")
    ax.set_title("insertion length distribution per TE family")
    ax.legend(frameon=False, fontsize=9)
    ax.axvline(6000, color="k", ls="--", lw=0.8)
    ax.text(6300, ax.get_ylim()[1] * 0.85, "full-length L1\n(>=6kb)",
            fontsize=8)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig3_context(df, out, fams=None):
    """Insertion count by genomic context, one panel per group.

    All three contexts are drawn on ONE log-count axis — the same unit and
    scale for exon/intron/intergenic — so the rare exon category (9-12
    events vs tens of thousands) is a real, visible bar instead of an
    invisible sliver on a percentage axis (which silently read as 'no exon
    insertions'). Zero-count buckets keep their slot (y row) with no bar.
    """
    fams = fams or FAM_ORDER
    groups = [g for g in GROUP_ORDER if (df["group"] == g).any()]
    if not groups:
        fig, _ax = plt.subplots(figsize=(6, 4))
        _empty_state(fig, out)
        return
    ctxs = ["exon", "intron", "intergenic"]
    ctx_colors = {"exon": "#C0392B", "intron": "#E67E22",
                  "intergenic": "#BDC3C7"}

    fig, axes = plt.subplots(1, len(groups),
                             figsize=(4.2 * len(groups), 4.5), sharey=True)
    if len(groups) == 1:
        axes = [axes]
    y = np.arange(len(fams))
    h = 0.8 / len(ctxs)
    for ax, g in zip(axes, groups):
        sub = df[(df.group == g) & (df.family.isin(fams))]
        tab = (sub.groupby(["family", "context"]).size()
               .unstack(fill_value=0).reindex(fams, fill_value=0))
        for ci, c in enumerate(ctxs):
            vals = tab[c].values if c in tab.columns else np.zeros(len(fams))
            pos = y + 0.4 - ci * h
            ax.barh(pos, vals, height=h * 0.92, color=ctx_colors[c])
            for fi, v in enumerate(vals):
                if v > 0:  # label every bar with its count
                    ax.text(v * 1.18, pos[fi], f"{int(v):,}", va="center",
                            ha="left", fontsize=7.5, color="#333333")
        ax.set_xscale("log")
        ax.set_xlim(1, 2e5)
        ax.set_yticks(y)
        ax.set_yticklabels(
            [f"{FAM_DISPLAY.get(f, f)}\n(n={int(tab.sum(axis=1)[f]):,})"
             for f in fams], fontsize=8)
        ax.set_xlabel("insertion count (log)")
        ax.set_title(group_label(g), fontsize=10)
    axes[0].set_ylabel("TE family")
    handles = [plt.Rectangle((0, 0), 1, 1, color=ctx_colors[c]) for c in ctxs]
    fig.legend(handles, ctxs, loc="upper center", ncol=3, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 1.00))
    fig.suptitle("insertions by genomic context", fontweight="bold", y=1.08)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig4_te_region(df, out, fams=None):
    fams = fams or FAM_ORDER
    groups = [g for g in GROUP_ORDER if (df["group"] == g).any()]
    te_flag = _clean_str(df["in_te_region"]) if "in_te_region" in df else \
        pd.Series([""] * len(df))
    df = df.assign(_in_te=te_flag != "")
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    if not groups or df.empty:
        _empty_state(fig, out)
        return
    x = np.arange(len(fams))
    w = 0.8 / len(groups)
    vtexts = []
    for gi, g in enumerate(groups):
        vals = []
        for f in fams:
            sub = df[(df.family == f) & (df.group == g)]
            vals.append(sub["_in_te"].mean() * 100 if len(sub) else 0.0)
        pos = x + gi * w - 0.4 + w / 2
        ax.bar(pos, vals, width=w * 0.92,
               color=GROUP_COLOR[g], label=group_label(g))
        for xi, v in zip(pos, vals):
            # "0" not "0.0" for zero bars: three adjacent decimal zeros
            # read as one fused string at this font size
            lab = "0" if v == 0 else f"{v:.1f}"
            vtexts.append(ax.text(xi, v, lab, ha="center",
                                  va="bottom", fontsize=7,
                                  color="#333333"))
    ax.set_xticks(x)
    labels = [FAM_DISPLAY.get(f, f) for f in fams]
    rot, ha = fit_ticks(ax, labels)
    ax.set_xticklabels(labels, rotation=rot, ha=ha,
                       rotation_mode="anchor", fontsize=8.5)
    ax.set_ylabel("% insertions inside reference TE region")
    # axis grows with the data instead of a fixed 0-105 that wastes 90%
    # of the panel when ratios are single-digit
    ymax = max([v for v in (float(t.get_text() or 0) for t in vtexts)] or [0.0])
    ax.set_ylim(0, max(ymax, 5.0) * 1.25)
    # geometry final -> legend band fit against real width -> title above
    # (second tight_layout reserves the padded title's room; no clipping)
    fig.tight_layout()
    leg, leg_h = legend_above(ax, ncol=max(len(groups), 1))
    ax.set_title("TE-into-TE vs unique-site insertions",
                 pad=title_pad(leg_h))
    fig.tight_layout()
    stack_value_labels(fig, vtexts)
    raise_ylim_for_labels(fig, ax, vtexts)
    _save(fig, out)


def fig5_full_length_l1(df, out, fl_file):
    fl = df[df["full_length_L1"] == True] if "full_length_L1" in df else \
        pd.read_csv(fl_file, sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    groups = [g for g in GROUP_ORDER if (fl.group == g).any()]
    if not groups:
        # zero candidates: emit the figure in an explicit empty state so the
        # stage never dies on degenerate data (skill: every declared output)
        for ax in axes:
            ax.text(0.5, 0.5, "no full-length L1 candidates",
                    ha="center", va="center", fontsize=11, color="#7F8C8D")
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle("full-length L1 candidates", fontweight="bold")
        fig.tight_layout()
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        return
    axes[0].hist([fl.loc[fl.group == g, "svlen"] for g in groups],
                 bins=30, stacked=True,
                 color=[GROUP_COLOR[g] for g in groups],
                 label=[group_label(g) for g in groups])
    axes[0].set_xlabel("SV length (bp)")
    axes[0].set_ylabel("full-length L1 candidates")
    axes[0].set_title("A  length distribution")
    axes[0].legend(frameon=False, fontsize=9)

    counts = [(fl.group == g).sum() for g in groups]
    axes[1].bar(range(len(groups)), counts,
                color=[GROUP_COLOR[g] for g in groups])
    for i, v in enumerate(counts):
        axes[1].text(i, v * 1.02, f"{v:,}", ha="center", fontsize=9)
    axes[1].set_xticks(range(len(groups)))
    axes[1].set_xticklabels([group_label(g) for g in groups], fontsize=9)
    axes[1].set_ylim(0, max(counts) * 1.18 if counts else 1)
    axes[1].set_title("B  per-group counts")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig6_polya_tsd(df, out, fams=None):
    # Fingerprints (polyA/TSD) are hallmarks of L1-mediated retrotransposition
    # — only retrotransposon MEI families carry the mechanism (profile flag
    # retro_fingerprint). LTR/DNA transpose by other mechanisms.
    fams = fams or [f for f in ["Alu", "L1", "SVA"]
                    if (df["family"] == f).any()]
    sub = df[df.family.isin(fams)].copy()
    sub["polyA_b"] = sub["polyA"].astype(str) == "True"
    sub["tsd_b"] = _clean_str(sub["tsd"]) != ""
    fig, axes = plt.subplots(1, 2, figsize=(max(10, 2.2 * len(fams) + 5),
                                            4.8), sharey=True)
    x = np.arange(len(fams))
    for ax, col, ttl in [(axes[0], "polyA_b", "poly(A) tail"),
                         (axes[1], "tsd_b", "target site duplication")]:
        vals = [sub.loc[sub.family == f, col].mean() * 100
                for f in fams]
        ax.bar(x, vals, color=[fam_color(f) for f in fams], width=0.6)
        for i, v in enumerate(vals):
            ax.text(i, min(v + 2.5, 96), f"{v:.1f}%", ha="center", fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(fams)
        ax.set_title(ttl)
        ax.set_ylim(0, 108)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
    axes[0].set_ylabel("% of insertions")
    fig.suptitle("MEI fingerprints (L1-mediated insertion hallmarks)",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out)
    plt.close(fig)


def _repel_labels(fig, labels, step: float = 6.0,
                  max_iter: int = 60) -> None:
    """Nudge overlapping point labels apart (vertical repulsion).

    Labels are ``annotate`` texts positioned in offset points; repeated
    pairwise passes push overlapping pairs apart until stable.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for _ in range(max_iter):
        boxes = [t.get_window_extent(renderer) for t in labels]
        moved = False
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if boxes[i].overlaps(boxes[j]):
                    xi, yi = labels[i].get_position()
                    labels[i].set_position((xi, yi + step / 2))
                    xj, yj = labels[j].get_position()
                    labels[j].set_position((xj, yj - step / 2))
                    moved = True
        if not moved:
            break
        fig.canvas.draw()


def fig7_diff_volcano(diff_tsv: str, out: str,
                      ctrl: str = "ctrl", exp: str = "exp") -> None:
    """Volcano of TE-repeat-level differences (mei_diff_by_repeat.tsv).

    x: Haldane-corrected log2 fold change of the exp vs ctrl call share;
    y: -log10 Fisher exact p. Significant points (p < 0.05) are black-
    edged and labelled names carry a " *" suffix for them. Labels: the
    top 5 per side (ctrl-/exp-enriched) by the joint score
    |log2FC| * -log10(p) — the "double-high" rule skips points extreme
    in only one dimension. Labels are auto-repelled. Exploratory: p
    values reflect call-count imbalance (no biological replicates).
    """
    tsv = diff_tsv.replace("mei_diff.tsv", "mei_diff_by_repeat.tsv")
    dd = pd.read_csv(tsv, sep="\t") if os.path.exists(tsv) else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    if dd.empty or "log2fc_exp_vs_ctrl" not in dd or "te_repeat" not in dd:
        _empty_state(fig, out)
        return
    x = dd["log2fc_exp_vs_ctrl"].astype(float)
    p = dd["p_value"].astype(float).clip(lower=1e-300)
    y = -np.log10(p)
    sig = dd["p_value"].astype(float) < 0.05
    up = x > 0
    ax.scatter(x[~up & ~sig], y[~up & ~sig], s=18, color="#4A72B7",
               alpha=0.8, label="ctrl-enriched", zorder=3)
    ax.scatter(x[up & ~sig], y[up & ~sig], s=18, color="#D73027",
               alpha=0.8, label="exp-enriched", zorder=3)
    ax.scatter(x[sig], y[sig], s=26,
               c=np.where(up[sig], "#D73027", "#4A72B7"),
               edgecolors="black", linewidths=0.6, zorder=4,
               label="p < 0.05")
    # labels: top 5 per side (ctrl-/exp-enriched) by the joint score
    # |log2FC| * -log10(p) — a "double-high" rule that skips points
    # extreme in only one dimension (huge n with tiny effect, or huge
    # effect backed by ~no evidence)
    score = x.abs() * y
    labels = []
    for sel in (x < 0, x >= 0):
        for i in score[sel].sort_values(ascending=False).head(5).index:
            name = str(dd.at[i, "te_repeat"])
            if float(dd.at[i, "p_value"]) < 0.05:
                name += " *"
            labels.append(ax.annotate(name, (float(x[i]), float(y[i])),
                                      fontsize=7, xytext=(5, 5),
                                      textcoords="offset points",
                                      color="#333333", zorder=5))
    _repel_labels(fig, labels)
    ax.set_xlabel(f"log2({exp}-specific / {ctrl}-specific call share)")
    ax.set_ylabel("-log10(Fisher exact p)")
    leg, leg_h = legend_above(ax, ncol=3)
    ax.set_title("TE repeat differential ranking (exploratory)",
                 pad=title_pad(leg_h))
    fig.tight_layout()
    _save(fig, out)


def html_report(meidir: str, title: str, figs, df, summary, out_html):
    fam_group = df.groupby(["family", "group"]).size().unstack(fill_value=0)
    prof = (summary or {}).get("te_buckets") or []
    te_names = [b.get("name") for b in prof
                if isinstance(b, dict) and b.get("name")
                and not b.get("catch_all") and b.get("name") != "microsatellite"]
    if not te_names:  # pre-profile mei dir
        te_names = [f for f in ["Alu", "L1", "SVA", "LTR", "DNA", "B1",
                                "B2/ID", "IAP/ETn"]
                    if f in set(df["family"].dropna())]
    ctx = (df[df.family.isin(te_names)]
           .groupby(["family", "context"]).size().unstack(fill_value=0))
    n_exon = int((df["context"] == "exon").sum())
    n_fl = int((df["full_length_L1"] == True).sum()) \
        if "full_length_L1" in df else 0
    rows_fam = "".join(
        f"<tr><td>{f}</td>" + "".join(f"<td>{int(fam_group.loc[f, g]):,}</td>"
                                      for g in fam_group.columns if g in fam_group.columns)
        + "</tr>" for f in fam_group.index)
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>MEI report — {title}</title>
<style>
 body {{ font-family: -apple-system, 'Noto Sans', sans-serif; margin: 2em; color:#222; }}
 h1 {{ font-size: 1.5em; border-bottom: 2px solid #4A72B7; padding-bottom:.3em; }}
 h2 {{ font-size: 1.15em; margin-top: 1.6em; color:#4A72B7; }}
 table {{ border-collapse: collapse; margin: .8em 0; font-size: .9em; }}
 th, td {{ border: 1px solid #ccc; padding: 4px 10px; text-align: right; }}
 th {{ background: #f0f3f8; }}
 td:first-child, th:first-child {{ text-align: left; }}
 img {{ max-width: 100%; border: 1px solid #eee; margin: .5em 0; }}
 .note {{ background:#fff8e6; border-left:4px solid #E08214; padding:.6em 1em; font-size:.9em; }}
 .kpis {{ display:flex; gap:1.2em; flex-wrap:wrap; }}
 .kpi {{ background:#f0f3f8; border-radius:8px; padding:.8em 1.2em; min-width:9em; }}
 .kpi b {{ display:block; font-size:1.6em; color:#4A72B7; }}
</style></head><body>
<h1>Mobile Element Insertion (MEI) analysis — {title}</h1>
<div class="kpis">
 <div class="kpi"><b>{summary['n_ins_total']:,}</b>insertions analyzed</div>
 <div class="kpi"><b>{int((df.family.isin(['Alu','L1','SVA','LTR','DNA'])).sum()):,}</b>TE-classified</div>
 <div class="kpi"><b>{n_fl:,}</b>full-length L1 candidates</div>
 <div class="kpi"><b>{n_exon}</b>exon insertions</div>
 <div class="kpi"><b>{summary.get('n_with_polyA',0):,}</b>with poly(A) tail</div>
</div>
<div class="note"><b>Caveat (1-vs-1 design):</b> each group has a single sample.
"Group-specific" insertions are dominated by individual polymorphism, NOT
group/age effect. Quantitative group comparisons require &ge;3 samples per
group. The near-symmetric family counts (see Fig 1) are consistent with
polymorphism, not directional activation.</div>
<h2>Family &times; group</h2>
<table><tr><th>family</th>{''.join(f'<th>{g}</th>' for g in fam_group.columns)}</tr>
{rows_fam}</table>
<h2>Figures</h2>
""" + "".join(f'<h2>{os.path.basename(f)}</h2><img src="data:image/png;base64,{b64(f)}">' for f in figs) + \
f"""
<h2>Method (reproducibility)</h2>
<ol style="font-size:.92em">
<li>Input: SURVIVOR-merged SV VCF (pbsv calls; INS ALT carries inserted sequence).</li>
<li>TE classification: RepeatMasker (-s -no_is -lib rhesus_te_lib.fa) on inserted
sequences &ge;100 bp; best hit by covered span; family bucket Alu/L1/SVA/LTR/DNA.</li>
<li>MEI fingerprints: poly(A) tail = &ge;12 A-run or &ge;90% A in last 60 bp;
TSD = identical 7–20 bp k-mer immediately 5' and 3' of the breakpoint (ref FASTA via pysam).</li>
<li>Context: gene/exon from GTF; pre-existing TE regions from rmsk lift-over BED.</li>
<li>Full-length L1 candidate: SVLEN &ge;6000 and L1 family.</li>
<li>Script: <code>modules/sv/bin/utils/mei_analysis.py</code> +
<code>modules/sv/bin/utils/mei_plot.py</code>; env: sv.sif container.</li>
</ol>
</body></html>"""
    with open(out_html, "w") as fh:
        fh.write(html)


def _relabeled(group, ctrl, exp):
    """Map internal ctrl_/exp_ prefixes to caller-provided branch names."""
    if group is None:
        return group
    if group.startswith("ctrl_"):
        return group.replace("ctrl_", ctrl + "_", 1)
    if group.startswith("exp_"):
        return group.replace("exp_", exp + "_", 1)
    return group


# ---------------------------------------------------------------------------
# Layout engine: compose first, then de-collide with measured geometry.
#
# Problem being solved (2026-11): one-off per-figure placements (legend
# bbox_to_anchor, fixed y headroom) broke whenever the data changed — long
# branch names pushed the fig1 legend into its title, unrotated fig4 tick
# labels fused into "DNAMicrosatellite…", zero markers crowded the axis
# floor while real bars went unlabeled. The fix is generic and measured,
# not per-dataset:
#
#   fit_ticks             pick tick rotation from rendered label widths
#   legend_above          park a legend in its own band above the axes,
#                         then place the title pad-measured above THAT
#   stack_value_labels    column-detect overlapping value labels and
#                         stack them with pixel-measured separation
#   raise_ylim_for_labels grow the y headroom until stacked labels fit
#   text_collisions       report any remaining text-text intersections
#
# All measurements use the real Agg renderer (pixel truth, no estimates),
# so the figures adapt to any branch names, tick counts or value ranges.
# ---------------------------------------------------------------------------

def _renderer(fig):
    fig.canvas.draw()
    return fig.canvas.get_renderer()


def _bbox(artist, rend):
    return artist.get_window_extent(renderer=rend)


def _inter(a, b, pad=1.0):
    """True if pixel bboxes a and b intersect (grown by pad px)."""
    return not (a.x1 + pad < b.x0 or b.x1 + pad < a.x0 or
                a.y1 + pad < b.y0 or b.y1 + pad < a.y0)


def _text_pt(fig, s, fontsize):
    """(width, height) of string s in points, measured by the renderer."""
    t = fig.text(0, 0, str(s), fontsize=fontsize)
    try:
        bb = t.get_window_extent(renderer=_renderer(fig))
        return bb.width * 72.0 / fig.dpi, bb.height * 72.0 / fig.dpi
    finally:
        t.remove()


def fit_ticks(ax, labels, fontsize=8.5):
    """(rotation, ha) so horizontal tick labels fit side by side.

    Measures the rendered width of every label and picks the first
    rotation in (0, 35, 50, 70, 90)deg whose projected per-tick advance
    fits the axis slot; 90deg terminates the ladder with zero horizontal
    advance, so a fit always exists for any label set.
    """
    if len(ax.get_xticks()) < 2:
        return 0, "center"
    fig = ax.figure
    trans = ax.transData.transform
    slot = abs(trans((1, 0))[0] - trans((0, 0))[0]) * 72.0 / fig.dpi
    for rot in (0, 35, 50, 70, 90):
        ha = "center" if rot == 0 else "right"
        widest = max((_text_pt(fig, s, fontsize)[0] for s in labels),
                     default=0.0)
        if widest * np.cos(np.radians(rot)) + 10 <= slot * 0.95:
            return rot, ha
    return 90, "right"


def _raise_px(t, dy_px):
    """Move Text t up by dy_px, preserving its coordinate system."""
    pos = t.get_position()
    if len(pos) == 3:  # offset-text 3-tuple: move the plain part only
        pos = pos[:2]
    tr = t.get_transform()
    x, y = tr.transform(pos)
    nx, ny = tr.inverted().transform((x, y + dy_px))
    t.set_position((float(nx), float(ny)))


def stack_value_labels(fig, texts, pad_px=2.5, rounds=4):
    """De-collide value labels by raising the upper one of each colliding
    pair (measured in pixels, cascading over rounds until stable).

    Labels only ever move up and only when their rendered boxes actually
    intersect, so non-overlapping neighbours are untouched and every
    label stays visually attached to its own bar. Works for any number
    of groups, label widths or value ranges. Returns moves made.
    """
    moved = 0
    for _ in range(rounds):
        rend = _renderer(fig)
        items = [(t, _bbox(t, rend)) for t in texts
                 if t.get_visible() and t.get_text()]
        items.sort(key=lambda p: p[1].y0)  # bottom-up
        moved = 0
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                ti, bi = items[i]
                tj, bj = items[j]
                if _inter(bi, bj, pad_px):
                    _raise_px(tj, bi.y1 + pad_px - bj.y0)
                    items[j] = (tj, _bbox(tj, rend))
                    moved += 1
        if not moved:
            break
    return moved


def raise_ylim_for_labels(fig, ax, texts, pad_px=3.0, max_iter=4):
    """Grow y headroom (log- or linear-safely) until all labels fit.

    Stacking can push labels above the axes box; each iteration expands
    the upper y-limit and re-stacks, so the two converge: more headroom
    lowers the bars, fewer labels need raising.
    """
    for _ in range(max_iter):
        rend = _renderer(fig)
        top_px = _bbox(ax, rend).y1
        worst = max((_bbox(t, rend).y1 for t in texts
                     if t.get_visible() and t.get_text()), default=top_px)
        if worst <= top_px - pad_px:
            return
        lo, hi = ax.get_ylim()
        if ax.get_yscale() == "log":
            ax.set_ylim(lo, hi * 1.4)
        else:
            ax.set_ylim(lo, lo + (hi - lo) * 1.4)
        stack_value_labels(fig, texts)


def legend_above(ax, ncol=None, fontsize=9, **kw):
    """Legend in its own band above the axes. Returns (legend, height_pt).

    The band [axes-top, axes-top + legend height] is dedicated to the
    legend; nothing else is placed there (the title is set pad-measured
    above it via title_pad), so legend/title collisions are impossible
    by construction, whatever the branch names or entry count. If the
    requested column count makes the legend wider than the axes, it is
    re-laid-out with fewer columns (growing upward instead of sideways)
    until it fits — measured, not guessed.
    """
    fig = ax.figure
    ax_w_pt = ax.get_window_extent().width * 72.0 / fig.dpi \
        if ax.get_window_extent().width else 0
    n = max(ncol or 1, 1)
    while True:
        leg = ax.legend(frameon=False, fontsize=fontsize, loc="lower left",
                        bbox_to_anchor=(0.0, 1.0), ncol=n,
                        handlelength=1.2, columnspacing=1.2, **kw)
        rend = _renderer(fig)
        w_pt = _bbox(leg, rend).width * 72.0 / fig.dpi
        if n == 1 or not ax_w_pt or w_pt <= ax_w_pt * 0.98:
            break
        n -= 1
        leg.remove()  # detach the too-wide legend from the artist tree
    rend = _renderer(fig)
    h_pt = _bbox(leg, rend).height * 72.0 / fig.dpi
    return leg, h_pt


def title_pad(leg_h_pt, gap_pt=5.0):
    """set_title pad that clears the legend band: legend height + gap."""
    return leg_h_pt + gap_pt


def text_collisions(fig, pad=0.5):
    """Remaining text-text intersections as (text_a, text_b) pairs.

    Renderer-measured audit of every visible Text in the figure,
    including tick labels, legend entries and annotation strings — the
    acceptance check for the layout engine (logged at debug level, and
    used by the regression harness to prove figures are collision-free).
    """
    rend = _renderer(fig)
    texts = [a for a in fig.findobj(matplotlib.text.Text)
             if a.get_visible() and a.get_text()]
    hits = []
    for i, a in enumerate(texts):
        for b in texts[i + 1:]:
            if _inter(_bbox(a, rend), _bbox(b, rend), pad):
                hits.append((a.get_text(), b.get_text()))
    return hits


def _save(fig, out):
    """Write the figure after a final collision audit (debug-logged)."""
    hits = text_collisions(fig)
    if hits:
        logger.debug("layout audit %s: %d text collision(s) remain: %s",
                     os.path.basename(out), len(hits),
                     "; ".join(f"[{a!r} vs {b!r}]" for a, b in hits[:6]))
    fig.savefig(out)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mei-dir", required=True,
                    help="output dir of mei_analysis.py")
    ap.add_argument("--title", default="youth_vs_aged")
    ap.add_argument("--ctrl-name", default=None,
                    help="control group display name (from "
                         "Params.sv.comparisons.<key>.control_group_name)")
    ap.add_argument("--exp-name", default=None,
                    help="experimental group display name (from "
                         "Params.sv.comparisons.<key>.experimental_group_name)")
    args = ap.parse_args()

    d = args.mei_dir
    df = pd.read_csv(os.path.join(d, "mei_events.tsv"), sep="\t",
                     dtype={"tsd": str, "in_te_region": str})
    # empty-string columns survive a TSV round-trip as NaN -> 'nan' strings;
    # normalize before any != '' test (fig4/fig6 read these columns)
    for col in ("tsd", "in_te_region", "context"):
        if col in df:
            df[col] = _clean_str(df[col])
    with open(os.path.join(d, "mei_summary.json")) as fh:
        summary = json.load(fh)
    # branch names: CLI wins, else fall back to what mei_analysis recorded
    branches = summary.get("branches") or {}
    ctrl = args.ctrl_name or branches.get("control") or "ctrl"
    exp = args.exp_name or branches.get("experiment") or "exp"
    if args.ctrl_name:
        GROUP_NAMES["ctrl"] = args.ctrl_name
    if args.exp_name:
        GROUP_NAMES["exp"] = args.exp_name
    # relabel the group column (ctrl_specific -> <ctrl>_specific) and keep the
    # default GROUP_ORDER/GROUP_LABEL machinery working on the new prefixes
    global GROUP_ORDER
    df["group"] = df["group"].apply(lambda g: _relabeled(g, ctrl, exp))
    GROUP_ORDER = [f"{ctrl}_specific", f"{exp}_specific", "shared"]
    GROUP_LABEL[f"{ctrl}_specific"] = f"{ctrl}-specific"
    GROUP_LABEL[f"{exp}_specific"] = f"{exp}-specific"
    GROUP_COLOR[f"{ctrl}_specific"] = GROUP_COLOR.pop("ctrl_specific", "#4A72B7")
    GROUP_COLOR[f"{exp}_specific"] = GROUP_COLOR.pop("exp_specific", "#D73027")

    figs = []
    slots, fp_fams = bucket_vocab(summary, df)
    jobs = [
        ("mei_family_by_group.png",
         lambda f: fig1_family_by_group(df, f, fams=slots)),
        ("mei_insert_length_by_family.png",
         lambda f: fig2_length_family(df, f, fams=slots)),
        ("mei_genomic_context.png",
         lambda f: fig3_context(df, f, fams=slots)),
        ("mei_te_region_ratio.png",
         lambda f: fig4_te_region(df, f, fams=slots)),
        ("mei_full_length_L1.png",
         lambda f: fig5_full_length_l1(df, f,
                                       os.path.join(d, "mei_full_length_L1.tsv"))),
        ("mei_polya_tsd_fingerprint.png",
         lambda f: fig6_polya_tsd(df, f, fams=fp_fams)),
        ("mei_diff_volcano.png",
         lambda f: fig7_diff_volcano(os.path.join(d, "mei_diff.tsv"), f,
                                     ctrl, exp)),
    ]
    for name, fn in jobs:
        out = os.path.join(d, name)
        fn(out)
        figs.append(out)
        print(f"[MEI-plot] wrote {out}")

    html_report(d, args.title, figs, df, summary,
                os.path.join(d, "mei_report.html"))
    print(f"[MEI-plot] wrote {os.path.join(d, 'mei_report.html')}")


if __name__ == "__main__":
    main()
