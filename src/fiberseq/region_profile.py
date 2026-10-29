#!/usr/bin/env python3
"""Region-class accessibility and methylation profiles for Fiber-seq samples.

For each region class in a region BED (col4 = class, e.g. majorSat_array,
tel_array_5p, MT_genome) and each sample (ft extract element BEDs), computes:

  - accessibility: fraction of region bp covered by MSPs >= --msp-min-len
    (pooled, and per-fiber mean), the >50 bp MSP proportion metric of
    Nat Genet 2025
  - nucleosome occupancy (pooled fraction)
  - mCpG / m6A mark density per kb of fiber read overlap
  - per-region table + class-level summary with two-sample Wilcoxon rank-sum
  - binned metaprofiles (region normalized to 0..1) per class

Element BED semantics follow ft extract BED12: name column = fiber (read)
name, blocks are the calls. For fibertools m6a/cpg BED12 the first and last
blocks are sentinels and are dropped; msp/nuc blocks are all real.

Usage:
  python region_profile.py -r regions.censat.bed \
      -s DMSO=D.msp.bed.gz,D.nuc.bed.gz,D.cpg.bed.gz,D.m6a.bed.gz \
      -s PlaB=P.msp.bed.gz,P.nuc.bed.gz,P.cpg.bed.gz,P.m6a.bed.gz \
      -o profile_out/
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402
from src.common.util.FiberlibUtil import (  # noqa: E402
    ChromSweep,
    read_bed,
    write_tsv,
    safe_div,
)


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("RegionProfile", log_file=str(log_file))
    else:
        lg = setup_logger("RegionProfile")
    for h in lg.handlers:
        if type(h) is logging.StreamHandler:
            h.setStream(sys.stdout)
    return lg


logger = _get_logger()


def iter_calls(rec, sentinels: bool) -> List[Tuple[int, int]]:
    """Interior blocks of a BED12 record (drop first/last when sentinels)."""
    blocks = list(rec.iter_blocks())
    if sentinels and len(blocks) > 2:
        return blocks[1:-1]
    if sentinels:
        return []
    return blocks


def _sweep(paths: Sequence[str], anchors, visitor) -> None:
    """Chrom-sweep *paths* against anchor regions, calling visitor(i, recs)."""
    sweep = ChromSweep(list(paths), list(anchors))
    for idx, _anchor, per_path in sweep.advance():
        visitor(idx, [r for r in per_path[0]])


# --------------------------------------------------------------------------- #
# per-region metric collection for one sample
# --------------------------------------------------------------------------- #
def collect_sample(regions, paths: Dict[str, str], msp_min_len: int) -> List[dict]:
    """One metrics row per region for one sample."""
    n = len(regions)
    anchors = [(c, s, e) for c, s, e, _ in regions]
    region_len = np.array([e - s for _, s, e, _ in regions], dtype=np.int64)

    msp_cov = np.zeros(n)
    nuc_cov = np.zeros(n)
    fiber_sets: List[Dict[str, List[Tuple[int, int]]]] = [dict() for _ in range(n)]
    call_counts = {"cpg": np.zeros(n, dtype=np.int64),
                   "m6a": np.zeros(n, dtype=np.int64)}
    overlap_bp = np.zeros(n, dtype=np.int64)

    # MSP: pooled coverage (len-filtered) + per-fiber covered bp
    def on_msp(i, recs):
        c, s, e, _ = regions[i]
        cov = np.zeros(e - s, dtype=bool)
        for r in recs:
            bits = fiber_sets[i].get(r.name)
            if bits is None:
                bits = []
                fiber_sets[i][r.name] = bits
            for bs, be in r.iter_blocks():
                if be - bs < msp_min_len:
                    continue
                s0, e0 = max(bs, s), min(be, e)
                if e0 > s0:
                    cov[s0 - s:e0 - s] = True
                    bits.append((s0, e0))
        msp_cov[i] = cov.sum()

    _sweep([paths["msp"]], anchors, on_msp)

    def on_nuc(i, recs):
        c, s, e, _ = regions[i]
        cov = np.zeros(e - s, dtype=bool)
        for r in recs:
            for bs, be in r.iter_blocks():
                s0, e0 = max(bs, s), min(be, e)
                if e0 > s0:
                    cov[s0 - s:e0 - s] = True
        nuc_cov[i] = cov.sum()

    _sweep([paths["nuc"]], anchors, on_nuc)

    for key in ("cpg", "m6a"):
        def on_meth(i, recs, key=key):
            cnt = 0
            for r in recs:
                cnt += len(iter_calls(r, sentinels=True))
                s0 = max(r.start, regions[i][1])
                e0 = min(r.end, regions[i][2])
                if e0 > s0:
                    overlap_bp[i] += e0 - s0
            call_counts[key][i] += cnt
        _sweep([paths[key]], anchors, on_meth)

    rows = []
    for i, (c, s, e, cls) in enumerate(regions):
        fracs = []
        rlen = e - s
        for spans in fiber_sets[i].values():
            # fast per-fiber union coverage via difference array
            diff = np.zeros(rlen + 1, dtype=np.int32)
            for s0, e0 in spans:
                diff[s0 - s] += 1
                diff[e0 - s] -= 1
            fracs.append(float((np.cumsum(diff[:rlen]) > 0).mean()) if rlen else 0.0)
        rows.append({
            "sample": paths["label"],
            "chrom": c, "start": s, "end": e, "class": cls,
            "length_bp": int(region_len[i]),
            "msp_frac_pooled": safe_div(msp_cov[i], region_len[i]),
            "msp_frac_per_fiber_mean": float(np.mean(fracs)) if fracs else 0.0,
            "n_fibers": len(fiber_sets[i]),
            "nuc_frac_pooled": safe_div(nuc_cov[i], region_len[i]),
            "m6a_calls_per_kb_overlap": safe_div(call_counts["m6a"][i] * 1000.0, overlap_bp[i]),
            "mcpg_calls_per_kb_overlap": safe_div(call_counts["cpg"][i] * 1000.0, overlap_bp[i]),
            "read_overlap_bp": int(overlap_bp[i]),
        })
    logger.info(f"region metrics done for sample {paths['label']}: {n} regions")
    return rows


# --------------------------------------------------------------------------- #
# binned metaprofile
# --------------------------------------------------------------------------- #
def binned_metaprofile(regions, msp_path: str, n_bins: int, msp_min_len: int
                       ) -> Dict[str, np.ndarray]:
    """Mean MSP-accessibility metaprofile per class (region normalized 0..1)."""
    anchors = [(c, s, e) for c, s, e, _ in regions]
    prof: Dict[str, List[np.ndarray]] = defaultdict(list)
    drop = n_bins * 20  # skip absurdly long regions for the metaprofile

    def on_msp(i, recs):
        c, s, e, cls = regions[i]
        if e - s > drop:
            return
        cov = np.zeros(e - s, dtype=bool)
        for r in recs:
            for bs, be in r.iter_blocks():
                if be - bs < msp_min_len:
                    continue
                s0, e0 = max(bs, s), min(be, e)
                if e0 > s0:
                    cov[s0 - s:e0 - s] = True
        usable = (e - s) // n_bins * n_bins
        if usable <= 0:
            return
        binned = cov[:usable].reshape(n_bins, -1).mean(axis=1)
        prof[cls].append(binned)

    _sweep([msp_path], anchors, on_msp)
    return {k: np.vstack(v).mean(axis=0) for k, v in prof.items() if v}


# --------------------------------------------------------------------------- #
# summary + plots
# --------------------------------------------------------------------------- #
METRICS = ["msp_frac_pooled", "msp_frac_per_fiber_mean", "nuc_frac_pooled",
           "m6a_calls_per_kb_overlap", "mcpg_calls_per_kb_overlap"]


def summarize(rows: List[dict]) -> List[dict]:
    """Class-level summary; rank-sum p between the first two samples."""
    from scipy.stats import ranksums

    groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in rows:
        groups[(r["sample"], r["class"])].append(r)
    samples = sorted({r["sample"] for r in rows})
    classes = sorted({r["class"] for r in rows})
    out: List[dict] = []
    for cls in classes:
        base: dict = {"class": cls}
        for smp in samples:
            g = groups.get((smp, cls), [])
            base[f"{smp}_n_regions"] = len(g)
            for m in METRICS:
                vals = [x[m] for x in g]
                base[f"{smp}_{m}_mean"] = float(np.mean(vals)) if vals else 0.0
                base[f"{smp}_{m}_median"] = float(np.median(vals)) if vals else 0.0
        if len(samples) == 2 and all(groups.get((s, cls)) for s in samples):
            for m in METRICS:
                a = [x[m] for x in groups[(samples[0], cls)]]
                b = [x[m] for x in groups[(samples[1], cls)]]
                try:
                    _stat, p = ranksums(a, b)
                except ValueError:
                    p = float("nan")
                base[f"{m}_ranksum_p"] = p
        out.append(base)
    return out


def plot_profiles(summary: List[dict], outdir: Path, samples: Sequence[str],
                  formats: Sequence[str]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = [r["class"] for r in summary]
    titles = [("msp_frac_per_fiber_mean", "MSP accessibility (per-fiber mean)"),
              ("nuc_frac_pooled", "Nucleosome occupancy"),
              ("mcpg_calls_per_kb_overlap", "mCpG / kb overlap"),
              ("m6a_calls_per_kb_overlap", "m6A / kb overlap")]
    width = 0.8 / max(1, len(samples))
    x = np.arange(len(classes))
    fig, axes = plt.subplots(4, 1, figsize=(max(7, 1.1 * len(classes)), 11))
    for ax, (m, title) in zip(axes, titles):
        for i, smp in enumerate(samples):
            vals = [r.get(f"{smp}_{m}_mean", 0.0) or 0.0 for r in summary]
            ax.bar(x + i * width - 0.4 + width / 2, vals, width, label=smp)
        ax.set_xticks(x)
        ax.set_xticklabels(classes, rotation=30, ha="right", fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8)
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(outdir / f"region_class_profile.{fmt}", dpi=300)
    plt.close(fig)


def plot_metaprofiles(profiles: Dict[str, Dict[str, np.ndarray]], outdir: Path,
                      n_bins: int, formats: Sequence[str]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = sorted({c for prof in profiles.values() for c in prof})
    if not classes:
        return
    fig, axes = plt.subplots(1, len(classes), figsize=(3.0 * len(classes), 3),
                             squeeze=False)
    xs = np.arange(n_bins) / max(1, n_bins - 1)
    for ax, cls in zip(axes[0], classes):
        for smp, prof in sorted(profiles.items()):
            if cls in prof:
                ax.plot(xs, prof[cls], marker=".", ms=3, label=smp)
        ax.set_title(cls, fontsize=9)
        ax.set_xlabel("normalized region position")
        ax.set_ylabel("MSP accessibility")
        ax.legend(fontsize=7)
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(outdir / f"region_metaprofile.{fmt}", dpi=300)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_sample(spec: str) -> Tuple[str, Dict[str, str]]:
    """Parse LABEL=MSP,NUC,CPG,M6A."""
    label, _, rest = spec.partition("=")
    parts = [p for p in rest.split(",") if p]
    if len(parts) != 4:
        raise SystemExit(f"--sample {label}: need 4 comma-separated paths: msp,nuc,cpg,m6a")
    d = dict(zip(["msp", "nuc", "cpg", "m6a"], parts))
    d["label"] = label
    return label, d


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="region_profile.py",
        description="Region-class accessibility/methylation profiles for Fiber-seq.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-r", "--regions", type=Path, required=True,
                   help="region BED with class in col4")
    p.add_argument("-s", "--sample", action="append", required=True,
                   metavar="LABEL=MSP,NUC,CPG,M6A",
                   help="sample label + 4 element BED paths (ft extract output)")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("-m", "--msp-min-len", type=int, default=50,
                   help="min MSP length for accessibility metric")
    p.add_argument("-n", "--n-bins", type=int, default=20,
                   help="bins for metaprofile")
    p.add_argument("--no-plot", action="store_true", help="skip plots")
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

    regions = [(r.chrom, r.start, r.end, r.name) for r in read_bed(str(args.regions))]
    logger.info(f"regions: {len(regions)} from {args.regions}")
    samples = dict(parse_sample(s) for s in args.sample)
    logger.info(f"samples: {list(samples)}")

    rows: List[dict] = []
    for label, paths in samples.items():
        rows.extend(collect_sample(regions, paths, args.msp_min_len))
    header = ["sample", "chrom", "start", "end", "class", "length_bp",
              "msp_frac_pooled", "msp_frac_per_fiber_mean", "n_fibers",
              "nuc_frac_pooled", "m6a_calls_per_kb_overlap",
              "mcpg_calls_per_kb_overlap", "read_overlap_bp"]
    write_tsv(str(args.outdir / "region_metrics.tsv"), header,
              [[r[c] for c in header] for r in rows])

    summary = summarize(rows)
    if summary:
        sheader = ["class"] + sorted({k for r in summary for k in r if k != "class"})
        write_tsv(str(args.outdir / "region_class_summary.tsv"), sheader,
                  [[r.get(c, "") for c in sheader] for r in summary])

    profiles: Dict[str, Dict[str, np.ndarray]] = {}
    for label, paths in samples.items():
        profiles[label] = binned_metaprofile(regions, paths["msp"],
                                             args.n_bins, args.msp_min_len)
        prof_rows = [(label, cls, i, float(v))
                     for cls, arr in sorted(profiles[label].items())
                     for i, v in enumerate(arr)]
        write_tsv(str(args.outdir / f"metaprofile_{label}.tsv"),
                  ["sample", "class", "bin", "msp_accessibility"], prof_rows)

    if not args.no_plot:
        if summary:
            plot_profiles(summary, args.outdir, list(samples), args.formats)
        plot_metaprofiles(profiles, args.outdir, args.n_bins, args.formats)
    logger.info(f"done -> {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
