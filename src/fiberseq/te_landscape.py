#!/usr/bin/env python3
"""TE-family chromatin landscape from Fiber-seq element BEDs (multi-mark).

For each TE copy in a TE BED (chrom start end name . strand family class,
the repo's rmsk projection layout) and each sample, computes per-copy:

  - accessibility: fraction of TE bp covered by MSPs >= --msp-min-len
  - nucleosome occupancy (pooled fraction)
  - mCpG and m6A mark density per kb of fiber read overlap

Then aggregates to class / family / subfamily level with two-sample
rank-sum p-values and a per-fiber-mean accessibility variant. Answers
"chromatin openness over repeats" and "5mC / 6mA distribution over TEs"
for DMSO vs PlaB style contrasts in one pass.

BED12 sentinel semantics: m6a/cpg first+last blocks are sentinels (dropped);
msp/nuc blocks are all real calls.

Usage:
  python te_landscape.py --te-bed GRCm39_GENCODE_rmsk_TE.lifted.bed \
      -s DMSO=D.msp,D.nuc,D.cpg,D.m6a -s PlaB=P.msp,P.nuc,P.cpg,P.m6a \
      -o te_out/ --level family
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
    ChromSweep, write_tsv, safe_div)


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("TeLandscape", log_file=str(log_file))
    else:
        lg = setup_logger("TeLandscape")
    for h in lg.handlers:
        if type(h) is logging.StreamHandler:
            h.setStream(sys.stdout)
    return lg


logger = _get_logger()


def load_te_bed(path: str) -> List[Tuple[str, int, int, str, str, str]]:
    """TE BED -> list of (chrom, start, end, subfamily, family, class)."""
    out: List[Tuple[str, int, int, str, str, str]] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            name = f[3] if len(f) > 3 else "."
            family = f[6] if len(f) > 6 else "."
            cls = f[7] if len(f) > 7 else "."
            out.append((f[0], s, e, name, family, cls))
    out.sort(key=lambda r: (r[0], r[1]))
    return out


def iter_calls(rec, sentinels: bool) -> List[Tuple[int, int]]:
    blocks = list(rec.iter_blocks())
    if sentinels and len(blocks) > 2:
        return blocks[1:-1]
    if sentinels:
        return []
    return blocks


def _sweep(paths: Sequence[str], anchors, visitor) -> None:
    sweep = ChromSweep(list(paths), list(anchors))
    for idx, _anchor, per_path in sweep.advance():
        visitor(idx, [r for r in per_path[0]])


def collect_sample(tes, paths: Dict[str, str], msp_min_len: int,
                   min_overlap_reads: int) -> List[dict]:
    """Per-TE-copy metrics for one sample (TEs with >= min reads kept)."""
    n = len(tes)
    anchors = [(c, s, e) for c, s, e, *_ in tes]

    msp_cov = np.zeros(n)
    nuc_cov = np.zeros(n)
    n_reads = np.zeros(n, dtype=np.int64)
    overlap_bp = np.zeros(n, dtype=np.int64)
    cpg = np.zeros(n, dtype=np.int64)
    m6a = np.zeros(n, dtype=np.int64)

    def on_msp(i, recs):
        c, s, e = tes[i][0], tes[i][1], tes[i][2]
        cov = np.zeros(e - s, dtype=bool)
        names = set()
        for r in recs:
            names.add(r.name)
            for bs, be in r.iter_blocks():
                if be - bs < msp_min_len:
                    continue
                s0, e0 = max(bs, s), min(be, e)
                if e0 > s0:
                    cov[s0 - s:e0 - s] = True
        msp_cov[i] = cov.sum()
        n_reads[i] = len(names)

    _sweep([paths["msp"]], anchors, on_msp)

    def on_nuc(i, recs):
        c, s, e = tes[i][0], tes[i][1], tes[i][2]
        cov = np.zeros(e - s, dtype=bool)
        for r in recs:
            for bs, be in r.iter_blocks():
                s0, e0 = max(bs, s), min(be, e)
                if e0 > s0:
                    cov[s0 - s:e0 - s] = True
        nuc_cov[i] = cov.sum()

    _sweep([paths["nuc"]], anchors, on_nuc)

    for key, acc in (("cpg", cpg), ("m6a", m6a)):
        def on_meth(i, recs, key=key, acc=acc):
            cnt = 0
            for r in recs:
                cnt += len(iter_calls(r, sentinels=True))
                s0 = max(r.start, tes[i][1])
                e0 = min(r.end, tes[i][2])
                if e0 > s0:
                    overlap_bp[i] += e0 - s0
            acc[i] += cnt
        _sweep([paths[key]], anchors, on_meth)

    rows = []
    label = paths["label"]
    for i, (c, s, e, sub, fam, cls) in enumerate(tes):
        if n_reads[i] < min_overlap_reads:
            continue
        rows.append({
            "sample": label, "chrom": c, "start": s, "end": e,
            "subfamily": sub, "family": fam, "class": cls,
            "length_bp": e - s,
            "n_read_overlaps": int(n_reads[i]),
            "msp_frac": safe_div(msp_cov[i], e - s),
            "nuc_frac": safe_div(nuc_cov[i], e - s),
            "mcpg_per_kb_overlap": safe_div(cpg[i] * 1000.0, overlap_bp[i]),
            "m6a_per_kb_overlap": safe_div(m6a[i] * 1000.0, overlap_bp[i]),
            "read_overlap_bp": int(overlap_bp[i]),
        })
    logger.info(f"{label}: kept {len(rows)}/{n} TE copies (>= {min_overlap_reads} reads)")
    return rows


METRICS = ["msp_frac", "nuc_frac", "mcpg_per_kb_overlap", "m6a_per_kb_overlap"]


def aggregate(rows: List[dict], level: str) -> List[dict]:
    """Aggregate per-copy rows to class/family/subfamily with rank-sum p."""
    from scipy.stats import ranksums

    key = {"class": "class", "family": "family", "subfamily": "subfamily"}[level]
    groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in rows:
        groups[(r["sample"], r[key])].append(r)
    samples = sorted({r["sample"] for r in rows})
    keys = sorted({r[key] for r in rows})
    out: List[dict] = []
    for k in keys:
        base: dict = {key: k}
        for smp in samples:
            g = groups.get((smp, k), [])
            base[f"{smp}_n_copies"] = len(g)
            base[f"{smp}_read_overlap_Mb"] = round(sum(x["read_overlap_bp"] for x in g) / 1e6, 3)
            for m in METRICS:
                vals = [x[m] for x in g]
                base[f"{smp}_{m}_mean"] = float(np.mean(vals)) if vals else 0.0
                base[f"{smp}_{m}_median"] = float(np.median(vals)) if vals else 0.0
        if len(samples) == 2 and all(groups.get((s, k)) for s in samples):
            for m in METRICS:
                a = [x[m] for x in groups[(samples[0], k)]]
                b = [x[m] for x in groups[(samples[1], k)]]
                try:
                    _st, p = ranksums(a, b)
                except ValueError:
                    p = float("nan")
                base[f"{m}_ranksum_p"] = p
                base[f"{m}_delta_median"] = base[f"{samples[1]}_{m}_median"] - base[f"{samples[0]}_{m}_median"]
        out.append(base)
    out.sort(key=lambda r: -r.get(f"{samples[0]}_read_overlap_Mb", 0))
    return out


def plot_heatmap(summary: List[dict], level: str, outdir: Path, samples: Sequence[str],
                 formats: Sequence[str], top_n: int) -> None:
    """Metric x TE-key heatmap of sample medians (top-N by overlap)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top = summary[:top_n]
    keys = [r[level] for r in top]
    cols = []
    labels = []
    for smp in samples:
        for m in METRICS:
            cols.append([r.get(f"{smp}_{m}_median", 0.0) or 0.0 for r in top])
            labels.append(f"{smp}\n{m}")
    mat = np.array(cols)
    fig, ax = plt.subplots(figsize=(0.35 * len(keys) + 3, 0.5 * len(labels) + 2))
    im = ax.imshow(mat, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(keys, rotation=90, fontsize=6)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.7)
    ax.set_title(f"TE {level} chromatin landscape")
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(outdir / f"te_{level}_heatmap.{fmt}", dpi=300)
    plt.close(fig)


def parse_sample(spec: str) -> Tuple[str, Dict[str, str]]:
    label, _, rest = spec.partition("=")
    parts = [p for p in rest.split(",") if p]
    if len(parts) != 4:
        raise SystemExit(f"--sample {label}: need 4 comma-separated paths: msp,nuc,cpg,m6a")
    d = dict(zip(["msp", "nuc", "cpg", "m6a"], parts))
    d["label"] = label
    return label, d


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="te_landscape.py",
        description="TE-family accessibility/5mC/6mA landscape from Fiber-seq.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--te-bed", required=True,
                   help="TE BED: chrom start end name . strand family class")
    p.add_argument("-s", "--sample", action="append", required=True,
                   metavar="LABEL=MSP,NUC,CPG,M6A")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("--level", default="family", choices=["class", "family", "subfamily"],
                   help="aggregation level for summary table + heatmap")
    p.add_argument("-m", "--msp-min-len", type=int, default=50)
    p.add_argument("--min-overlap-reads", type=int, default=5,
                   help="min overlapping fibers to keep a TE copy")
    p.add_argument("--top-n", type=int, default=50, help="heatmap top-N rows")
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

    tes = load_te_bed(args.te_bed)
    logger.info(f"TE copies: {len(tes)}")
    samples = dict(parse_sample(s) for s in args.sample)

    rows: List[dict] = []
    for label, paths in samples.items():
        rows.extend(collect_sample(tes, paths, args.msp_min_len, args.min_overlap_reads))
    header = ["sample", "chrom", "start", "end", "subfamily", "family", "class",
              "length_bp", "n_read_overlaps", "msp_frac", "nuc_frac",
              "mcpg_per_kb_overlap", "m6a_per_kb_overlap", "read_overlap_bp"]
    write_tsv(str(args.outdir / "te_copy_metrics.tsv"), header,
              [[r[c] for c in header] for r in rows])

    samples_sorted = sorted(samples)
    for level in (["class", "family", "subfamily"]
                  if args.level == "family" else [args.level]):
        summary = aggregate(rows, level)
        sheader = [level] + sorted({k for r in summary for k in r if k != level})
        write_tsv(str(args.outdir / f"te_{level}_summary.tsv"), sheader,
                  [[r.get(c, "") for c in sheader] for r in summary])
        if not args.no_plot and summary:
            plot_heatmap(summary, level, args.outdir, samples_sorted,
                         args.formats, args.top_n)
    logger.info(f"done -> {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
