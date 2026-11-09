#!/usr/bin/env python3
"""Group comparison of per-region Fiber-seq modification metrics.

Joins per-sample region quant tables (quant_matrix.py) and per-region
accessibility/mCpG tables (censat_profile.py) across two sample groups, and
writes per-region and per-region-class difference tables plus deliverable
figures:

  region_diff.tsv    per region: both groups' values, diff (g2-g1), log2fc
  class_diff.tsv     per class x metric: mean/sd per group, diff, log2fc,
                     paired Wilcoxon p over regions (region-level consistency,
                     NOT biological replication)
  fig_class_metrics  grouped bars per class with per-region points
  fig_region_scatter g1 vs g2 per region, colored by class
  fig_profiles_*     class-averaged accessibility / mCpG profiles per group

Usage:
    python region_mod_compare.py --meta meta.tsv \
        --quant S1=path/S1.region_quant.tsv --quant S2=path/S2.region_quant.tsv \
        --regions-tsv S1=path/S1.censat_regions.tsv ... \
        --profile-tsv S1=path/S1.censat_profile.tsv ... \
        --g1 youth --g2 aged -o compare/
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_BIN_DIR = Path(__file__).resolve().parent
_ROOT_DIR = _BIN_DIR.parent.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

METRICS = [
    ("m6a_per_fiber", "m6A per fiber"),
    ("percent_accessible", "accessibility (MSP fiber fraction)"),
    ("nuc_occupancy", "nucleosome occupancy (fiber fraction)"),
    ("fire_frac", "FIRE fiber fraction"),
    ("mcpg", "mCpG marks per kb"),
    ("msp_pooled_frac", "MSP pooled coverage fraction"),
]
EPS = 1e-6
GROUP_COLORS = {"g1": "#4C72B0", "g2": "#DD8452"}


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("RegionModCompare", log_file=str(log_file))
    else:
        lg = setup_logger("RegionModCompare")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def read_tsv(path: Path) -> Tuple[List[str], List[List[str]]]:
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return header, rows


def to_float(s: str) -> float:
    try:
        return float(s)
    except ValueError:
        return float("nan")


def parse_pairs(specs: Sequence[str]) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    for spec in specs:
        sample, _, path = spec.partition("=")
        out[sample] = Path(path)
    return out


def load_meta(path: Path) -> Dict[str, str]:
    meta: Dict[str, str] = {}
    header, rows = read_tsv(path)
    isample = header.index("sample_id")
    igroup = header.index("group")
    for r in rows:
        meta[r[isample]] = r[igroup]
    return meta


def region_key(f: Sequence[str]) -> Tuple[str, str, str, str]:
    return f[0], f[1], f[2], f[3]


def load_quant(path: Path) -> Dict[Tuple, Dict[str, float]]:
    header, rows = read_tsv(path)
    out = {}
    for r in rows:
        out[region_key(r)] = {
            "m6a_per_fiber": to_float(r[header.index("m6a_per_fiber")]),
            "percent_accessible": to_float(r[header.index("percent_accessible")]),
            "nuc_occupancy": to_float(r[header.index("nuc_occupancy")]),
            "fire_frac": to_float(r[header.index("fire_frac")]),
        }
    return out


def load_censat_regions(path: Path) -> Dict[Tuple, Dict[str, float]]:
    header, rows = read_tsv(path)
    out = {}
    for r in rows:
        out[region_key(r)] = {
            "mcpg": to_float(r[header.index("mcpg")]),
            "msp_pooled_frac": to_float(r[header.index("msp_pooled_frac")]),
        }
    return out


def load_profiles(path: Path
                  ) -> Dict[str, Dict[int, List[Tuple[float, float]]]]:
    """class -> position -> [(accessibility, mcpg), ...] per region bin."""
    header, rows = read_tsv(path)
    out: Dict[str, Dict[int, List[Tuple[float, float]]]] = defaultdict(
        lambda: defaultdict(list))
    for r in rows:
        out[r[0]][int(r[3])].append((to_float(r[4]), to_float(r[5])))
    return out


def wilcoxon_paired(vals1: "np.ndarray", vals2: "np.ndarray") -> float:
    from scipy.stats import wilcoxon
    diffs = np.array(vals2) - np.array(vals1)
    diffs = diffs[np.isfinite(diffs)]
    if len(diffs) < 2 or np.all(diffs == 0):
        return float("nan")
    try:
        return float(wilcoxon(diffs).pvalue)
    except ValueError:
        return float("nan")


def build_tables(g1: str, g2: str,
                 samples_by_group: Dict[str, List[str]],
                 quant: Dict[str, Dict], censat: Dict[str, Dict],
                 outdir: Path) -> Tuple[List, List[Tuple], Dict]:
    """Write region_diff.tsv + class_diff.tsv; return class rows for figures."""
    # region universe: union of quant keys (same regions across samples)
    keys: List[Tuple] = []
    seen = set()
    for table in list(quant.values()) + list(censat.values()):
        for k in table:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    keys.sort(key=lambda k: (k[0], int(k[1]), int(k[2]), k[3]))

    def group_value(tables: Dict[str, Dict], group: str, k: Tuple,
                    metric: str) -> float:
        vals = [tables[s][k][metric] for s in samples_by_group[group]
                if s in tables and k in tables[s]]
        vals = [v for v in vals if np.isfinite(v)]
        return float(np.mean(vals)) if vals else float("nan")

    region_rows = []
    per_class_metric: Dict[Tuple[str, str], Dict[str, List]] = defaultdict(
        lambda: {"g1": [], "g2": []})
    for k in keys:
        chrom, start, end, cls = k
        row = [chrom, start, end, cls]
        for metric, _ in METRICS:
            # metric lives in exactly one of the two table families
            fam = quant if metric in ("m6a_per_fiber", "percent_accessible",
                                      "nuc_occupancy", "fire_frac") else censat
            v1 = group_value(fam, g1, k, metric)
            v2 = group_value(fam, g2, k, metric)
            row += [v1, v2, v2 - v1,
                    float(np.log2((v2 + EPS) / (v1 + EPS)))
                    if np.isfinite(v1) and np.isfinite(v2) else float("nan")]
            per_class_metric[(cls, metric)]["g1"].append(v1)
            per_class_metric[(cls, metric)]["g2"].append(v2)
        region_rows.append(row)

    header = ["chrom", "start", "end", "region_class"]
    for metric, _ in METRICS:
        header += [f"{metric}_{g1}", f"{metric}_{g2}",
                   f"{metric}_diff", f"{metric}_log2fc"]
    out1 = outdir / "region_diff.tsv"
    with open(out1, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in region_rows:
            fh.write("\t".join(
                x if isinstance(x, str) else f"{x:.6g}" for x in row) + "\n")
    logger.info(f"Wrote {out1} ({len(region_rows)} regions)")

    classes = sorted({k[3] for k in keys})
    class_rows = []
    for cls in classes:
        for metric, _ in METRICS:
            a = np.array(per_class_metric[(cls, metric)]["g1"], dtype=float)
            b = np.array(per_class_metric[(cls, metric)]["g2"], dtype=float)
            m1, m2 = np.nanmean(a), np.nanmean(b)
            s1 = np.nanstd(a, ddof=1) if len(a) > 1 else float("nan")
            s2 = np.nanstd(b, ddof=1) if len(b) > 1 else float("nan")
            p = wilcoxon_paired(a, b)
            class_rows.append([
                cls, metric, int(np.isfinite(a).sum()),
                m1, s1, m2, s2, m2 - m1,
                float(np.log2((m2 + EPS) / (m1 + EPS)))
                if np.isfinite(m1) and np.isfinite(m2) else float("nan"),
                p,
            ])
    out2 = outdir / "class_diff.tsv"
    with open(out2, "w") as fh:
        fh.write("\t".join([
            "region_class", "metric", "n_regions",
            f"mean_{g1}", f"sd_{g1}", f"mean_{g2}", f"sd_{g2}",
            "diff", "log2fc", "wilcoxon_p_paired_regions"]) + "\n")
        for row in class_rows:
            fh.write("\t".join(
                x if isinstance(x, str) else f"{x:.6g}" for x in row) + "\n")
    logger.info(f"Wrote {out2} ({len(class_rows)} class x metric rows)")
    return class_rows, keys, per_class_metric


def savefig(fig, outbase: Path) -> None:
    for ext in ("png", "svg"):
        fig.savefig(outbase.with_suffix(f".{ext}"), dpi=300, facecolor="white",
                    bbox_inches="tight", pad_inches=0.3)


def fig_class_metrics(per_class_metric: Dict, keys: List[Tuple], g1: str,
                      g2: str, outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = sorted({k[3] for k in keys})
    x = np.arange(len(classes))
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, (metric, label) in zip(axes.flat, METRICS):
        a = [per_class_metric[(c, metric)]["g1"] for c in classes]
        b = [per_class_metric[(c, metric)]["g2"] for c in classes]
        m1 = [np.nanmean(v) if v else np.nan for v in a]
        m2 = [np.nanmean(v) if v else np.nan for v in b]
        w = 0.36
        ax.bar(x - w / 2, m1, w, color=GROUP_COLORS["g1"],
               edgecolor="black", linewidth=0.6, label=g1)
        ax.bar(x + w / 2, m2, w, color=GROUP_COLORS["g2"],
               edgecolor="black", linewidth=0.6, label=g2)
        # per-region points
        rng = np.random.default_rng(0)
        for xi, vals, off in ((x - w / 2, a, -w / 2), (x + w / 2, b, w / 2)):
            for xi_, v in zip(xi, vals):
                arr = np.array(v, dtype=float)
                arr = arr[np.isfinite(arr)]
                if len(arr):
                    jitter = rng.uniform(-w / 4, w / 4, len(arr))
                    ax.plot(np.full(len(arr), xi_) + jitter, arr, ".",
                            color="black", alpha=0.35, markersize=3,
                            zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(classes, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel(label, fontsize=9)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if metric == METRICS[0][0]:
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Telomere/centromere region modifications by class\n"
                 f"({g1} vs {g2})", fontsize=12, y=1.02)
    fig.tight_layout()
    savefig(fig, outdir / "fig_class_metrics")
    plt.close(fig)


def fig_region_scatter(quant: Dict, censat: Dict, samples_by_group: Dict,
                       keys: List[Tuple], g1: str, g2: str, outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = sorted({k[3] for k in keys})
    cmap = plt.get_cmap("tab10")
    colors = {c: cmap(i % 10) for i, c in enumerate(classes)}
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.5))
    for ax, (metric, label) in zip(axes.flat, METRICS):
        fam = quant if metric in ("m6a_per_fiber", "percent_accessible",
                                  "nuc_occupancy", "fire_frac") else censat

        def gval(k, g):
            vals = [fam[s][k][metric] for s in samples_by_group[g]
                    if s in fam and k in fam[s]]
            return np.nanmean(vals) if vals else np.nan

        for cls in classes:
            xs, ys = [], []
            for k in keys:
                if k[3] != cls:
                    continue
                v1, v2 = gval(k, g1), gval(k, g2)
                if np.isfinite(v1) and np.isfinite(v2):
                    xs.append(v1)
                    ys.append(v2)
            if xs:
                ax.plot(xs, ys, ".", color=colors[cls], label=cls, markersize=5,
                        alpha=0.8)
        lims = [np.nanmin(ax.get_xlim() + ax.get_ylim()),
                np.nanmax(ax.get_xlim() + ax.get_ylim())]
        if np.all(np.isfinite(lims)):
            ax.plot(lims, lims, "--", color="grey", linewidth=0.8)
            ax.set_xlim(lims)
            ax.set_ylim(lims)
        ax.set_xlabel(f"{g1}: {label}", fontsize=8)
        ax.set_ylabel(f"{g2}: {label}", fontsize=8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if metric == METRICS[0][0]:
            ax.legend(frameon=False, fontsize=7, markerscale=2)
    fig.suptitle(f"Per-region {g1} vs {g2} concordance by class", fontsize=12,
                 y=1.01)
    fig.tight_layout()
    savefig(fig, outdir / "fig_region_scatter")
    plt.close(fig)


def fig_profiles(profiles: Dict[str, Dict], samples_by_group: Dict, keys,
                 g1: str, g2: str, outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = sorted({k[3] for k in keys})
    for metric_idx, metric_name, fname in (
            (0, "MSP fiber coverage (unique fibers per base)",
             "fig_profiles_accessibility"),
            (1, "mCpG marks per kb", "fig_profiles_mcpg")):
        fig, axes = plt.subplots(2, int(np.ceil(len(classes) / 2)),
                                 figsize=(14, 7), squeeze=False)
        for ax, cls in zip(axes.flat, classes):
            for g in (g1, g2):
                # per-position mean across sample curves (linear aggregation)
                agg: Dict[int, List[float]] = defaultdict(list)
                for s in samples_by_group[g]:
                    posmap = profiles.get(s, {}).get(cls, {})
                    for pos, vals in posmap.items():
                        v = float(np.nanmean([x[metric_idx] for x in vals]))
                        if np.isfinite(v):
                            agg[pos].append(v)
                if agg:
                    xs = sorted(agg)
                    ys = [float(np.mean(agg[p])) for p in xs]
                    ax.plot(xs, ys,
                            color=GROUP_COLORS["g1" if g == g1 else "g2"],
                            label=g)
            ax.axvline(0, ls="--", lw=0.7, color="grey")
            ax.set_title(cls, fontsize=9)
            ax.set_xlabel("position vs region center (bp)", fontsize=8)
            ax.set_ylabel(metric_name, fontsize=8)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
        for ax in axes.flat[len(classes):]:
            ax.axis("off")
        axes.flat[0].legend(frameon=False, fontsize=8)
        fig.suptitle(f"Region-class {metric_name} profiles\n({g1} vs {g2})",
                     fontsize=12, y=1.02)
        fig.tight_layout()
        savefig(fig, outdir / fname)
        plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="region_mod_compare.py",
        description="Two-group comparison of region modification metrics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--meta", type=Path, required=True)
    p.add_argument("--quant", action="append", default=[],
                   help="sample_id=region_quant.tsv (repeatable)")
    p.add_argument("--regions-tsv", action="append", default=[],
                   help="sample_id=censat_regions.tsv (repeatable)")
    p.add_argument("--profile-tsv", action="append", default=[],
                   help="sample_id=censat_profile.tsv (repeatable)")
    p.add_argument("--g1", required=True, help="first group label")
    p.add_argument("--g2", required=True, help="second group label")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("-l", "--log", type=Path, default=None)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")

    meta = load_meta(args.meta)
    samples_by_group = {args.g1: [], args.g2: []}
    for s in list(parse_pairs(args.quant)) + list(parse_pairs(args.regions_tsv)):
        g = meta.get(s)
        if g in samples_by_group and s not in samples_by_group[g]:
            samples_by_group[g].append(s)
    for g, ss in samples_by_group.items():
        if not ss:
            logger.error(f"no samples for group '{g}' in meta")
            return 1
    logger.info(f"groups: {samples_by_group}")

    quant = {s: load_quant(p) for s, p in parse_pairs(args.quant).items()}
    censat = {s: load_censat_regions(p)
              for s, p in parse_pairs(args.regions_tsv).items()}
    profiles = {s: load_profiles(p)
                for s, p in parse_pairs(args.profile_tsv).items()}

    args.outdir.mkdir(parents=True, exist_ok=True)
    _, keys, per_class_metric = build_tables(
        args.g1, args.g2, samples_by_group, quant, censat, args.outdir)
    fig_class_metrics(per_class_metric, keys, args.g1, args.g2, args.outdir)
    fig_region_scatter(quant, censat, samples_by_group, keys, args.g1,
                       args.g2, args.outdir)
    if profiles:
        fig_profiles(profiles, samples_by_group, keys, args.g1, args.g2,
                     args.outdir)
    logger.info(f"All comparison outputs written to {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
