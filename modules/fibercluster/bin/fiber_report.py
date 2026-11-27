#!/usr/bin/env python3
"""Curate fiber_cluster.py run outputs into a reader-facing report package.

Consumes one or more run directories (each produced by fiber_cluster.py with
one channel/k setting) and writes:

    {outdir}/figures/<region>.clusters.png     # canonical run, regions with >=2 clusters
    {outdir}/summary/loci_summary.tsv          # ALL regions incl. single-cluster ones
    {outdir}/summary/cluster_proportions.png   # cross-locus cluster-fraction overview
    {outdir}/sensitivity/sensitivity_table.tsv # cluster structure across settings

Regions whose fibers form a single cluster (no detectable heterogeneity in
that setting) get a summary-table row instead of a figure.

Usage:
    python fiber_report.py \
        --run fire_k2=.../fiber_cluster_pilot/fire_k2 \
        --run fire_k3=.../fiber_cluster_pilot/fire_k3 \
        --run k2=.../fiber_cluster_pilot/k2 \
        --run k3=.../fiber_cluster_pilot/k3 \
        --canonical fire_k2 \
        --outdir .../fiber_cluster_pilot
"""
from __future__ import annotations

import argparse
import csv
import logging
import shutil
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from typing import TextIO  # noqa: E402


def setup_logger(
    logger_name: str,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    stream: Optional[TextIO] = None,
) -> logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


logger = setup_logger("FiberReport", stream=sys.stdout)


class RegionStats(NamedTuple):
    """Per-region cluster structure read from one run directory."""
    region: str
    n_fibers: int
    n_per_sample: Dict[str, int]
    clusters: List[int]
    sizes: Dict[int, int]
    counts: Dict[Tuple[int, str], int]


def read_run(run_dir: Path) -> Dict[str, RegionStats]:
    """Parse every region's cluster_summary.tsv under one run directory.

    Parameters:
        run_dir: a fiber_cluster.py output directory (one region per subdir).

    Returns:
        region name -> RegionStats.

    Errors:
        ValueError when a summary file is malformed.
    """
    out: Dict[str, RegionStats] = {}
    for sub in sorted(p for p in run_dir.iterdir() if p.is_dir()):
        summary = sub / "cluster_summary.tsv"
        if not summary.is_file():
            continue
        counts: Dict[Tuple[int, str], int] = {}
        n_per_sample: Dict[str, int] = {}
        sizes: Dict[int, int] = {}
        region = sub.name
        with open(summary) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                region = row["region"]
                c, s = int(row["cluster"]), row["sample"]
                n = int(row["n_fibers"])
                counts[(c, s)] = n
                sizes[c] = sizes.get(c, 0) + n
                n_per_sample[s] = n_per_sample.get(s, 0) + n
        out[region] = RegionStats(region, sum(n_per_sample.values()), n_per_sample,
                                  sorted(sizes), sizes, counts)
    return out


def parse_run_entries(entries: Sequence[str]) -> Dict[str, Path]:
    """Parse repeated '<label>=<run_dir>' CLI entries.

    Parameters:
        entries: raw CLI strings in '<label>=<path>' form.

    Returns:
        label -> validated run directory.

    Errors:
        ValueError on malformed entry or missing directory.
    """
    out: Dict[str, Path] = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError(f"--run: expected '<label>=<dir>', got {entry!r}")
        label, path_str = entry.split("=", 1)
        path = Path(path_str)
        if not path.is_dir():
            raise ValueError(f"--run: not a directory for label {label!r}: {path}")
        out[label] = path
    return out


def write_loci_summary(stats: Dict[str, RegionStats], out_tsv: Path) -> None:
    """Write the all-loci summary table (including single-cluster loci).

    Parameters:
        stats: canonical run's per-region stats.
        out_tsv: output table path.
    """
    with open(out_tsv, "w") as fh:
        fh.write("region\tn_fibers\tn_clusters\tcluster_sizes\tsample_counts\tnote\n")
        for region, st in stats.items():
            sizes = ",".join(f"{c}:{st.sizes[c]}" for c in st.clusters)
            scounts = ",".join(f"{s}:{st.n_per_sample[s]}" for s in sorted(st.n_per_sample))
            note = ("single cluster: no detectable fiber-state heterogeneity in this setting"
                    if len(st.clusters) < 2 else "see figures/<region>.clusters.png")
            fh.write(f"{region}\t{st.n_fibers}\t{len(st.clusters)}\t{sizes}\t{scounts}\t{note}\n")


def write_sensitivity_table(runs: Dict[str, Path], per_run: Dict[str, Dict[str, RegionStats]],
                            out_tsv: Path) -> None:
    """Write the cross-setting sensitivity table.

    Parameters:
        runs: label -> run directory.
        per_run: label -> region -> stats.
        out_tsv: output table path.
    """
    with open(out_tsv, "w") as fh:
        fh.write("run\tregion\tn_fibers\tn_clusters\tcluster_sizes\tper_cluster_sample_counts\n")
        for label in runs:
            for region, st in per_run[label].items():
                sizes = ",".join(f"{c}:{st.sizes[c]}" for c in st.clusters)
                detail = ";".join(
                    f"c{c}:" + ",".join(f"{s}={st.counts.get((c, s), 0)}"
                                        for s in sorted(st.n_per_sample))
                    for c in st.clusters)
                fh.write(f"{label}\t{region}\t{st.n_fibers}\t{len(st.clusters)}\t{sizes}\t{detail}\n")


def plot_overview(stats: Dict[str, RegionStats], sample_order: List[str],
                  out_png: Path) -> None:
    """Cross-locus cluster-fraction overview (the conclusion figure).

    For every region with >=2 clusters: one grouped bar per sample showing
    the within-sample cluster fractions (stacked, sums to 1). Single-cluster
    regions are excluded (they carry no class structure to compare).

    Parameters:
        stats: canonical run's per-region stats.
        sample_order: display order of samples (reference first).
        out_png: output image path.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    regions = [r for r, st in stats.items() if len(st.clusters) >= 2]
    if not regions:
        logger.warning("no region with >=2 clusters; overview figure skipped")
        return
    n_reg, n_sam = len(regions), len(sample_order)
    fig_w = max(6.0, 2.4 * n_reg * n_sam)
    fig, ax = plt.subplots(figsize=(fig_w, 4.2))
    width = 0.8 / n_sam
    for i, region in enumerate(regions):
        st = stats[region]
        for j, sample in enumerate(sample_order):
            base = i + (j - (n_sam - 1) / 2) * width
            bottom = 0.0
            n_s = st.n_per_sample.get(sample, 0)
            for c in st.clusters:
                frac = st.counts.get((c, sample), 0) / n_s if n_s else 0.0
                ax.bar(base, frac, width=width * 0.92, bottom=bottom,
                       color=plt.cm.tab10((c - 1) % 10), edgecolor="white", lw=0.4)
                bottom += frac
    ax.set_xticks([i + (j - (n_sam - 1) / 2) * width
                   for i in range(n_reg) for j in range(n_sam)])
    ax.set_xticklabels([f"{region}\n{sample}" for region in regions for sample in sample_order],
                       fontsize=7)
    ax.set_ylabel("fraction within sample")
    ax.set_ylim(0, 1)
    handles = [Rectangle((0, 0), 1, 1, facecolor=plt.cm.tab10((c - 1) % 10))
               for c in sorted({c for st in stats.values() for c in st.clusters})]
    labels = [f"cluster {c}" for c in sorted({c for st in stats.values() for c in st.clusters})]
    ax.legend(handles, labels, fontsize=7, loc="lower right")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser.

    Returns:
        Configured argparse.ArgumentParser.
    """
    p = argparse.ArgumentParser(
        description="Curate fiber_cluster.py runs into figures/ + summary/ + sensitivity/.")
    p.add_argument("-r", "--run", action="append", required=True, metavar="L=DIR",
                   help="label=run_dir entry (repeatable; one per channel/k setting)")
    p.add_argument("-c", "--canonical", type=str, default=None,
                   help="label of the run whose figures/summary become the deliverable "
                        "(default: the first --run label)")
    p.add_argument("-o", "--outdir", type=Path, required=True,
                   help="report package directory (e.g. the pilot root)")
    p.add_argument("-S", "--sample-order", type=str, default="",
                   help="comma-separated sample display order (reference first)")
    p.add_argument("-l", "--log", type=Path, default=None, help="optional log file path")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Build the curated report package.

    Parameters:
        argv: CLI arguments (None = sys.argv[1:]).

    Returns:
        Process exit code: 0 on success, 2 on parameter/input errors.
    """
    global logger
    args = build_parser().parse_args(argv)
    if args.log is not None:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        logger = setup_logger("FiberReport", log_file=str(args.log), stream=sys.stdout)

    try:
        runs = parse_run_entries(args.run)
        canonical = args.canonical or next(iter(runs))
        if canonical not in runs:
            raise ValueError(f"--canonical {canonical!r} is not among --run labels "
                             f"{sorted(runs)}")
    except ValueError as exc:
        logger.error(str(exc))
        return 2

    per_run = {label: read_run(path) for label, path in runs.items()}
    canon = per_run[canonical]
    if not canon:
        logger.error("canonical run %s has no region outputs", canonical)
        return 2

    if args.sample_order:
        sample_order = [s for s in args.sample_order.split(",") if s]
    else:
        sample_order = sorted({s for st in canon.values() for s in st.n_per_sample})

    fig_dir = args.outdir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    copied = 0
    for region, st in canon.items():
        if len(st.clusters) < 2:
            logger.info("region %s: single cluster -> summary only, no figure", region)
            continue
        src = runs[canonical] / region / f"{region}.clusters.png"
        if not src.is_file():
            logger.warning("region %s: figure missing at %s", region, src)
            continue
        shutil.copy2(src, fig_dir / src.name)
        copied += 1

    sum_dir = args.outdir / "summary"
    sum_dir.mkdir(parents=True, exist_ok=True)
    write_loci_summary(canon, sum_dir / "loci_summary.tsv")
    plot_overview(canon, sample_order, sum_dir / "cluster_proportions.png")

    sens_dir = args.outdir / "sensitivity"
    sens_dir.mkdir(parents=True, exist_ok=True)
    write_sensitivity_table(runs, per_run, sens_dir / "sensitivity_table.tsv")

    logger.info("report package written: %d figures, %d loci (summary), "
                "%d settings (sensitivity)", copied, len(canon), len(runs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
