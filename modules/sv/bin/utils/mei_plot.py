#!/usr/bin/env python
"""Visualization suite for mei_analysis.py results.

Reads the mei/ outputs (mei_events.tsv etc.) and renders the publication-
style figure set:

  mei_fig1_family_by_group.png    family x group stacked bar + group totals
  mei_fig2_length_family.png      insertion length distribution per family
  mei_fig3_context.png            genomic context (exon/intron/intergenic)
                                  per family x group (100% stacked)
  mei_fig4_te_region.png          TE-into-TE vs unique-site insertion ratio
  mei_fig5_full_length_L1.png     full-length L1 candidates: length hist +
                                  per-group counts
  mei_fig6_polya_tsd.png          MEI fingerprint rates (polyA, TSD) per family
  mei_report.html                 self-contained summary report embedding
                                  all figures and key tables

Usage:
  python mei_plot.py --mei-dir <mei_outdir> [--title "youth vs aged"]
"""

import argparse
import base64
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

GROUP_ORDER = ["ctrl_specific", "exp_specific", "shared"]
GROUP_LABEL = {"ctrl_specific": "group1-specific\n(ctrl)",
               "exp_specific": "group2-specific\n(exp)",
               "shared": "shared"}
GROUP_COLOR = {"ctrl_specific": "#4A72B7", "exp_specific": "#D73027",
               "shared": "#7F7F7F"}
FAM_ORDER = ["Alu", "L1", "SVA", "LTR", "DNA", "simple_repeat", "other",
             "unclassified", "short_unclassified"]
FAM_COLOR = {"Alu": "#E08214", "L1": "#4A72B7", "SVA": "#8E44AD",
             "LTR": "#27AE60", "DNA": "#16A085", "simple_repeat": "#95A5A6",
             "other": "#BDC3C7", "unclassified": "#7F8C8D",
             "short_unclassified": "#D5D8DC"}


def b64(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def fig1_family_by_group(df, out):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5),
                             gridspec_kw={"width_ratios": [2, 1]})
    fams = [f for f in FAM_ORDER if (df["family"] == f).any()]
    groups = [g for g in GROUP_ORDER if (df["group"] == g).any()]

    # left: grouped bars
    x = np.arange(len(fams))
    w = 0.8 / len(groups)
    for gi, g in enumerate(groups):
        vals = [(df[(df.family == f) & (df.group == g)]).shape[0] for f in fams]
        axes[0].bar(x + gi * w - 0.4 + w / 2, vals, width=w * 0.92,
                    color=GROUP_COLOR[g], label=GROUP_LABEL[g].replace("\n", " "))
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(fams, rotation=30, ha="right")
    axes[0].set_ylabel("insertions")
    axes[0].set_yscale("log")
    axes[0].legend(frameon=False)
    axes[0].set_title("A  TE family x group (log scale)")

    # right: totals per group
    tot = [ (df.group == g).sum() for g in groups]
    axes[1].bar(range(len(groups)), tot,
                color=[GROUP_COLOR[g] for g in groups])
    for i, v in enumerate(tot):
        axes[1].text(i, v * 1.02, f"{v:,}", ha="center", fontsize=9)
    axes[1].set_xticks(range(len(groups)))
    axes[1].set_xticklabels([GROUP_LABEL[g] for g in groups], fontsize=8)
    axes[1].set_ylabel("total insertions")
    axes[1].set_title("B  group totals")
    ymax = max(tot)
    axes[1].set_ylim(0, ymax * 1.15)

    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig2_length_family(df, out):
    fams = ["Alu", "L1", "SVA", "LTR"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.logspace(1.7, 4.3, 45)
    for f in fams:
        vals = df.loc[df.family == f, "ins_len"].clip(lower=50)
        if len(vals) == 0:
            continue
        ax.hist(vals, bins=bins, alpha=0.55, label=f"{f} (n={len(vals):,})",
                color=FAM_COLOR[f])
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


def fig3_context(df, out):
    fams = [f for f in ["Alu", "L1", "SVA", "LTR", "DNA"]
            if (df["family"] == f).any()]
    groups = [g for g in GROUP_ORDER if (df["group"] == g).any()]
    ctxs = ["exon", "intron", "intergenic"]
    ctx_colors = {"exon": "#C0392B", "intron": "#E67E22",
                  "intergenic": "#BDC3C7"}

    fig, axes = plt.subplots(1, len(groups), figsize=(4.2 * len(groups), 4.5),
                             sharey=True)
    if len(groups) == 1:
        axes = [axes]
    for ax, g in zip(axes, groups):
        sub = df[(df.group == g) & (df.family.isin(fams))]
        tab = (sub.groupby(["family", "context"]).size()
               .unstack(fill_value=0).reindex(fams, fill_value=0))
        pct = tab.div(tab.sum(axis=1), axis=0) * 100
        bottom = np.zeros(len(fams))
        y = np.arange(len(fams))
        for c in ctxs:
            if c in pct.columns:
                vals = pct[c].values if c in pct.columns else np.zeros(len(fams))
            else:
                vals = np.zeros(len(fams))
            ax.barh(y, vals, left=bottom, color=ctx_colors[c], label=c,
                    height=0.7)
            bottom += vals
        ax.set_yticks(y)
        ax.set_yticklabels([f"{f}\n(n={int(tab.sum(axis=1)[f]):,})" for f in fams],
                           fontsize=8)
        ax.set_xlabel("% of insertions")
        ax.set_title(GROUP_LABEL[g], fontsize=10)
        ax.set_xlim(0, 100)
    axes[0].set_ylabel("TE family")
    handles = [plt.Rectangle((0, 0), 1, 1, color=ctx_colors[c]) for c in ctxs]
    axes[-1].legend(handles, ctxs, frameon=False, fontsize=8,
                    loc="lower right")
    fig.suptitle("genomic context of TE insertions", fontweight="bold")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig4_te_region(df, out):
    fams = [f for f in ["Alu", "L1", "SVA", "LTR", "DNA"]
            if (df["family"] == f).any()]
    groups = [g for g in GROUP_ORDER if (df["group"] == g).any()]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    x = np.arange(len(fams))
    w = 0.8 / len(groups)
    for gi, g in enumerate(groups):
        vals = []
        for f in fams:
            sub = df[(df.family == f) & (df.group == g)]
            if len(sub) == 0 or "in_te_region" not in sub:
                vals.append(0.0)
                continue
            in_te = (sub["in_te_region"].astype(str) != "").mean() * 100
            vals.append(in_te)
        ax.bar(x + gi * w - 0.4 + w / 2, vals, width=w * 0.92,
               color=GROUP_COLOR[g], label=GROUP_LABEL[g].replace("\n", " "))
    ax.set_xticks(x)
    ax.set_xticklabels(fams)
    ax.set_ylabel("% insertions inside reference TE region")
    ax.set_title("TE-into-TE vs unique-site insertions")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig5_full_length_l1(df, out, fl_file):
    fl = df[df["full_length_L1"] == True] if "full_length_L1" in df else \
        pd.read_csv(fl_file, sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    groups = [g for g in GROUP_ORDER if (fl.group == g).any()]
    axes[0].hist([fl.loc[fl.group == g, "svlen"] for g in groups],
                 bins=30, stacked=True,
                 color=[GROUP_COLOR[g] for g in groups],
                 label=[GROUP_LABEL[g].replace("\n", " ") for g in groups])
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
    axes[1].set_xticklabels([GROUP_LABEL[g] for g in groups], fontsize=8)
    axes[1].set_ylim(0, max(counts) * 1.18 if counts else 1)
    axes[1].set_title("B  per-group counts")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig6_polya_tsd(df, out):
    fams = [f for f in ["Alu", "L1", "SVA", "LTR", "DNA"]
            if (df["family"] == f).any()]
    sub = df[df.family.isin(fams)].copy()
    sub["polyA_b"] = sub["polyA"].astype(str) == "True"
    sub["tsd_b"] = sub["tsd"].astype(str) != ""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    x = np.arange(len(fams))
    for ax, col, ttl in [(axes[0], "polyA_b", "poly(A) tail"),
                         (axes[1], "tsd_b", "target site duplication")]:
        vals = [sub.loc[sub.family == f, col].mean() * 100 for f in fams]
        ax.bar(x, vals, color=[FAM_COLOR[f] for f in fams], width=0.65)
        for i, v in enumerate(vals):
            ax.text(i, v + 1.5, f"{v:.1f}%", ha="center", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(fams, rotation=20, ha="right")
        ax.set_title(ttl)
        ax.set_ylim(0, 100)
    axes[0].set_ylabel("% of insertions")
    fig.suptitle("MEI fingerprints (L1-mediated insertion hallmarks)",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def html_report(meidir, title, figs, df, summary, out_html):
    fam_group = df.groupby(["family", "group"]).size().unstack(fill_value=0)
    ctx = (df[df.family.isin(["Alu", "L1", "SVA", "LTR", "DNA"])]
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mei-dir", required=True,
                    help="output dir of mei_analysis.py")
    ap.add_argument("--title", default="youth_vs_aged")
    args = ap.parse_args()

    d = args.mei_dir
    df = pd.read_csv(os.path.join(d, "mei_events.tsv"), sep="\t",
                     dtype={"tsd": str, "in_te_region": str})
    with open(os.path.join(d, "mei_summary.json")) as fh:
        summary = json.load(fh)

    figs = []
    jobs = [
        ("mei_fig1_family_by_group.png", lambda f: fig1_family_by_group(df, f)),
        ("mei_fig2_length_family.png", lambda f: fig2_length_family(df, f)),
        ("mei_fig3_context.png", lambda f: fig3_context(df, f)),
        ("mei_fig4_te_region.png", lambda f: fig4_te_region(df, f)),
        ("mei_fig5_full_length_L1.png",
         lambda f: fig5_full_length_l1(df, f,
                                       os.path.join(d, "mei_full_length_L1.tsv"))),
        ("mei_fig6_polya_tsd.png", lambda f: fig6_polya_tsd(df, f)),
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
