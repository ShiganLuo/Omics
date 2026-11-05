#!/usr/bin/env python3
"""TE insertion dating by divergence-from-consensus (molecular clock).

Promoted from a one-off script that lived inside the mei/te_age output dir into
the shared SV utils (2026-11). Now takes CLI args instead of hard-coded
relative paths, but the method is unchanged, so re-running on an existing
mei/ output reproduces the previous te_age_* files.

Method: per insertion, take RepeatMasker hits matching the best-hit repeat
name, compute span-weighted mean divergence (% substitutions vs consensus),
apply Jukes-Cantor correction, convert to age via T = d / (2*mu).
mu (per-lineage neutral substitution rate) is an ASSUMPTION: default
1.0e-9 /site/year, sensitivity 0.8 / 1.5 e-9. Group comparisons
(youth-specific vs aged-specific vs shared) are rate-independent.

Usage:
  python te_age_analysis.py --mei-dir <mei_outdir> [-o <te_age_outdir>]
                            [--families Alu,L1,LTR,SVA,DNA]

Input (all under --mei-dir):
  mei_events.tsv          per-event table from mei_analysis.py
  repeatmasker/ins.fa.out RepeatMasker raw output (divergence column)
  repeatmasker/ins.fa     inserted-sequence FASTA (id order = events >=min_len)

Output (in --outdir, default <mei-dir>/te_age):
  te_divergence_per_ins.tsv  per-insertion divergence + ages at 3 mu values
  te_age_by_group.tsv        family x group stats + MWU/KS p-values
  te_age_subfamily.tsv       per-subfamily divergence table
  te_full_length_L1_age.tsv  full-length L1 (SVLEN>=6000) with ages
  te_age figures (divergence_hist/landscape/divergence_box/full_length_L1)         figs (hist/landscape/box/full-length-L1)
  te_age_report.md           digest report
"""
import argparse
import math
import os
import sys
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MU_LIST = [("mu0.8e-9", 0.8e-9), ("mu1.0e-9", 1.0e-9), ("mu1.5e-9", 1.5e-9)]
GROUP_LABEL = {"ctrl_specific": "youth-specific", "exp_specific": "aged-specific",
               "shared": "shared"}


# ---------------- parse RepeatMasker .out (same right-anchored walk as
# mei_analysis.py, + div) ----------------
def parse_out(path):
    hits = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split()
            if len(f) < 15 or not f[0].isdigit():
                continue
            try:
                sw = int(f[0]); div = float(f[1])
                sid = f[4]; begin, end = int(f[5]), int(f[6])
            except (ValueError, IndexError):
                continue

            def _num(t):
                return t.isdigit()

            def _pnum(t):
                body = t.strip("()")
                return t.startswith("(") and body.lstrip("-").isdigit()

            k = len(f) - 1
            if f[k] == "*":
                k -= 1
            if k < 8 or not f[k].rstrip("*").isdigit():
                continue
            j = k - 1
            while j >= 8 and (_num(f[j]) or _pnum(f[j])):
                j -= 1
            if j < 8:
                continue
            cls = f[j]; rep = f[j - 1].rstrip("?")
            hits[sid].append({"begin": begin, "end": end, "repeat": rep,
                              "class": cls, "div": div, "sw": sw})
    return hits


def jc(p):
    """Jukes-Cantor corrected substitutions/site from observed p (0-1)."""
    if p is None or p >= 0.7 or p < 0:
        return None
    return -0.75 * math.log(1.0 - 4.0 * p / 3.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mei-dir", required=True,
                    help="mei_analysis.py output dir (contains mei_events.tsv "
                         "and repeatmasker/)")
    ap.add_argument("--min-len", type=int, default=100,
                    help="same insertion-length cutoff mei_analysis used "
                         "(events >= this fill repeatmasker/ins.fa)")
    ap.add_argument("-o", "--outdir", default=None,
                    help="output dir (default: <mei-dir>/te_age)")
    ap.add_argument("--families", default="Alu,L1,LTR,SVA,DNA",
                    help="comma-separated families to analyze, in priority "
                         "order; fig1/fig2 panels take the first 3 with data")
    args = ap.parse_args()

    mei = os.path.abspath(args.mei_dir)
    outdir = os.path.abspath(args.outdir or os.path.join(mei, "te_age"))
    os.makedirs(outdir, exist_ok=True)
    events = os.path.join(mei, "mei_events.tsv")
    rm_out = os.path.join(mei, "repeatmasker", "ins.fa.out")
    rm_fa = os.path.join(mei, "repeatmasker", "ins.fa")
    for p in (events, rm_out, rm_fa):
        if not os.path.exists(p):
            sys.exit(f"FATAL: missing input {p}")
    fams = [f.strip() for f in args.families.split(",") if f.strip()]
    top3 = [f for f in fams if f != "SVA"][:3]  # SVA: usually too few events
    if len(top3) < 3:  # fall back to raw order if families were removed
        top3 = fams[:3]

    ev = pd.read_csv(events, sep="\t", keep_default_na=False)
    hits = parse_out(rm_out)

    # ---- map insN ids: i-th event row with ins_len>=min_len (VCF record order)
    # = ins{i} ----
    n_fa = 0
    with open(rm_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                n_fa += 1
    sub = ev[ev["ins_len"] >= args.min_len].copy().reset_index(drop=True)
    if len(sub) != n_fa:
        sys.exit(f"FATAL: events>=100bp ({len(sub)}) != ins.fa seqs ({n_fa})")
    sub["ins_id"] = [f"ins{i}" for i in range(len(sub))]
    print(f"[ok] joined {len(sub)} classified events to ins.fa ids")

    # ---- per-insertion divergence ----
    rows, mism = [], 0
    for _, r in sub.iterrows():
        hs = hits.get(r["ins_id"], [])
        if not hs:
            rows.append({"ins_id": r["ins_id"], "n_hits": 0,
                         "best_repeat": "", "w_div": None})
            if r["te_repeat"] != "none":
                mism += 1
            continue
        hs_sorted = sorted(hs, key=lambda h: h["end"] - h["begin"], reverse=True)
        best = hs_sorted[0]
        if best["repeat"] != r["te_repeat"]:
            mism += 1
        same = [h for h in hs if h["repeat"] == best["repeat"]]
        wsum = sum(h["end"] - h["begin"] + 1 for h in same)
        wdiv = sum(h["div"] * (h["end"] - h["begin"] + 1) for h in same) / wsum
        rows.append({"ins_id": r["ins_id"], "n_hits": len(hs),
                     "best_repeat": best["repeat"], "w_div": round(wdiv, 3),
                     "best_span_div": best["div"], "n_frag": len(same)})
    dd = pd.DataFrame(rows)
    print(f"[ok] best-repeat mismatch vs mei_events te_repeat: {mism}/{len(sub)} "
          f"({100*mism/len(sub):.2f}%)")

    sub = sub.merge(dd, on="ins_id", how="left")
    sub["p_div"] = sub["w_div"] / 100.0
    sub["jc_div"] = sub["p_div"].apply(jc)
    for col, mu in MU_LIST:
        sub[col + "_age_myr"] = sub["jc_div"].apply(
            lambda d: round(d / (2 * mu) * 1e-6, 3) if d is not None and pd.notna(d) else None)

    keep = ["chrom", "pos", "id", "svlen", "ins_len", "supp_vec", "group",
            "family", "te_repeat", "ins_id", "n_hits", "n_frag", "w_div",
            "best_span_div", "jc_div"] + [c + "_age_myr" for c, _ in MU_LIST]
    if "polyA" in sub.columns:
        keep += ["polyA", "tsd"]
    sub[keep].to_csv(os.path.join(outdir, "te_divergence_per_ins.tsv"),
                     sep="\t", index=False)

    # ---- stats by family x group ----
    try:
        from scipy.stats import mannwhitneyu, ks_2samp
        HAVE_SCIPY = True
    except ImportError:
        HAVE_SCIPY = False
        print("[warn] scipy missing; p-values skipped")

    stat_rows = []
    for fam in fams:
        base = sub[(sub["family"] == fam) & sub["w_div"].notna()]
        for g in ["ctrl_specific", "exp_specific", "shared"]:
            b = base[base["group"] == g]["w_div"]
            if len(b) == 0:
                continue
            stat_rows.append({
                "family": fam, "group": GROUP_LABEL[g], "n": len(b),
                "median_div": round(b.median(), 3), "mean_div": round(b.mean(), 3),
                "frac_div<0.5": round((b < 0.5).mean(), 4),
                "frac_div<1.0": round((b < 1.0).mean(), 4),
                "frac_div<2.0": round((b < 2.0).mean(), 4)})
        if HAVE_SCIPY:
            a = base[base["group"] == "ctrl_specific"]["w_div"].dropna()
            b = base[base["group"] == "exp_specific"]["w_div"].dropna()
            s = base[base["group"] == "shared"]["w_div"].dropna()
            if len(a) > 10 and len(b) > 10:
                u = mannwhitneyu(a, b, alternative="two-sided")
                k = ks_2samp(a, b)
                stat_rows.append({"family": fam, "group": "youth-vs-aged MWU p",
                                  "n": len(a) + len(b), "median_div": round(u.pvalue, 4),
                                  "mean_div": round(k.pvalue, 4)})
            if len(s) > 10 and len(a) + len(b) > 10:
                pool = pd.concat([a, b])
                k2 = ks_2samp(pool, s)
                stat_rows.append({"family": fam, "group": "specific-vs-shared KS p",
                                  "n": len(pool) + len(s), "median_div": round(k2.pvalue, 4)})
    stat = pd.DataFrame(stat_rows)
    stat.to_csv(os.path.join(outdir, "te_age_by_group.tsv"), sep="\t", index=False)

    # ---- subfamily table ----
    sf = (sub[sub["w_div"].notna()]
          .groupby(["te_repeat", "family", "group"])
          .agg(n=("id", "size"), median_div=("w_div", "median"),
               mean_div=("w_div", "mean")).reset_index())
    sf = sf.sort_values(["family", "n"], ascending=[True, False])
    sf.to_csv(os.path.join(outdir, "te_age_subfamily.tsv"), sep="\t", index=False)

    # ---- full-length L1 ----
    fl_cols = ["chrom", "pos", "id", "group", "te_repeat", "ins_len", "w_div",
               "jc_div", "mu1.0e-9_age_myr"]
    fl = sub[(sub["family"] == "L1") & (sub["svlen"] >= 6000) &
             sub["w_div"].notna()][fl_cols + (["polyA", "tsd"]
                                              if "polyA" in sub.columns else [])]
    fl.to_csv(os.path.join(outdir, "te_full_length_L1_age.tsv"), sep="\t", index=False)

    # ================= figures =================
    colors = {"ctrl_specific": "#4C72B0", "exp_specific": "#C44E52",
              "shared": "#55A868"}

    # fig1: divergence histograms per family per group + young zoom
    fig, axes = plt.subplots(2, len(top3), figsize=(5 * len(top3), 8),
                             squeeze=False)
    for c, fam in enumerate(top3):
        base = sub[(sub["family"] == fam) & sub["w_div"].notna()]
        for g in ["ctrl_specific", "exp_specific", "shared"]:
            b = base[base["group"] == g]["w_div"]
            axes[0, c].hist(b, bins=60, range=(0, 30), alpha=0.55,
                            label=f"{GROUP_LABEL[g]} (n={len(b)})",
                            color=colors[g])
            axes[1, c].hist(b, bins=40, range=(0, 3), alpha=0.55,
                            color=colors[g])
        axes[0, c].set_title(f"{fam}: divergence from consensus")
        axes[1, c].set_title(f"{fam}: zoom 0-3%")
        axes[1, c].set_xlabel("Kimura div (%)")
        for r_ in (0, 1):
            axes[r_, c].set_ylabel("insertions")
        axes[0, c].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "te_age_divergence_hist.png"), dpi=150)
    plt.close(fig)

    # fig2: age landscape: div bins x family, stacked by group, fraction
    bins = [(0, 0.5), (0.5, 1), (1, 2), (2, 5), (5, 10), (10, 20), (20, 100)]
    labels = ["<0.5", "0.5-1", "1-2", "2-5", "5-10", "10-20", ">20"]
    fig, axes = plt.subplots(1, len(top3), figsize=(5 * len(top3), 4.5),
                             squeeze=False)
    for c, fam in enumerate(top3):
        base = sub[(sub["family"] == fam) & sub["w_div"].notna()]
        bottom = None
        for g in ["ctrl_specific", "exp_specific", "shared"]:
            b = base[base["group"] == g]["w_div"]
            fr = [((b >= lo) & (b < hi)).sum() / max(len(b), 1)
                  for lo, hi in bins]
            axes[0, c].bar(labels, fr, bottom=bottom, alpha=0.8,
                           label=GROUP_LABEL[g], color=colors[g])
            bottom = fr if bottom is None else [x + y for x, y in zip(bottom, fr)]
        axes[0, c].set_title(f"{fam} (n={len(base)})")
        axes[0, c].set_ylabel("fraction")
        axes[0, c].tick_params(axis="x", rotation=45)
    axes[0, 0].legend()
    fig.suptitle("Divergence-bin landscape by group (rate-free relative age)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "te_age_landscape.png"), dpi=150)
    plt.close(fig)

    # fig3: boxplots of divergence by group per family (all analyzed fams)
    fig, ax = plt.subplots(figsize=(9, 5))
    data, ticks = [], []
    for fam in fams:
        base = sub[(sub["family"] == fam) & sub["w_div"].notna()]
        for g in ["ctrl_specific", "exp_specific", "shared"]:
            b = base[base["group"] == g]["w_div"]
            if len(b):
                data.append(b)
                ticks.append(f"{fam}\n{GROUP_LABEL[g]}\nn={len(b)}")
    if data:
        bp = ax.boxplot(data, tick_labels=ticks, showfliers=False,
                        patch_artist=True)
        for patch, tick in zip(bp["boxes"], [t.split("\n")[1] for t in ticks]):
            patch.set_facecolor({"youth-specific": "#4C72B0",
                                 "aged-specific": "#C44E52",
                                 "shared": "#55A868"}[tick])
    ax.set_ylabel("divergence from consensus (%)")
    ax.set_title("TE insertion age proxy by family and group")
    plt.xticks(rotation=60, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "te_age_divergence_box.png"), dpi=150)
    plt.close(fig)

    # fig4: full-length L1 divergence
    fig, ax = plt.subplots(figsize=(8, 5))
    data, ticks = [], []
    for g in ["ctrl_specific", "exp_specific", "shared"]:
        b = fl[fl["group"] == g]["w_div"]
        if len(b):
            data.append(b)
            ticks.append(f"full-length L1\n{GROUP_LABEL[g]} n={len(b)}")
    if data:
        bp = ax.boxplot(data, tick_labels=ticks, patch_artist=True)
        for patch, i in zip(bp["boxes"], range(len(data))):
            patch.set_facecolor(list(colors.values())[i % 3])
        ax.set_ylabel("divergence (%)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "te_age_full_length_L1.png"), dpi=150)
    plt.close(fig)

    # ---- markdown digest report ----
    grp = stat[stat["group"].isin(["youth-specific", "aged-specific", "shared"])]
    lines = ["# TE insertion age (divergence-from-consensus) report", "",
             f"- mei dir: `{mei}`", f"- insertions dated: {len(sub)}",
             "- age model: JC-corrected divergence / (2*mu); mu assumptions "
             + ", ".join(f"{c}={v:.1e}" for c, v in MU_LIST), "",
             "## median divergence by family x group", ""]
    if len(grp):
        hdr = list(grp.columns)
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("|" + "---|" * len(hdr))
        for _, rr in grp.iterrows():
            lines.append("| " + " | ".join(str(v) for v in rr) + " |")
    lines += ["", "## notes", "",
              "- group comparisons (MWU/KS) are rate-independent; absolute "
              "ages in Myr depend on the assumed mu.",
              "- SVA typically has too few events for its own fig1/fig2 panel "
              "but is included in fig3 boxplots and the stats tables."]
    with open(os.path.join(outdir, "te_age_report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # ---- console digest ----
    print("\n=== median divergence (age proxy) by family x group ===")
    print(grp.to_string(index=False))
    print(f"\nfull-length L1 with div: n={len(fl)}, "
          f"median={fl['w_div'].median() if len(fl) else 'NA'}")
    print("\n[done] outputs in", outdir)


if __name__ == "__main__":
    main()
