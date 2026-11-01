#!/usr/bin/env python3
"""scCODA differential abundance analysis (pertpy).

Bayesian Dirichlet-Multinomial compositional analysis of cell-type
abundance across two conditions (Büttner et al. 2021, Nat Comm), following
https://www.sc-best-practices.org/conditions/compositional/.

Why scCODA over Fisher's exact:
  - Jointly models all cell-type proportions (compositional constraint).
  - Reference cell type auto-selected (least dispersion across samples).
  - Bayesian credible-effect output, designed for low-replicate data.
  - Supports covariate formulas (batch correction, multi-factor designs).

Inputs per tissue: an h5ad with ``obs['sample_id']`` and ``obs['cell_type']``
(e.g. the scanpy_auto output). Each comparison names a reference sample
(LEFT, baseline) and a treatment sample (RIGHT, perturbed); effects are
reported as log2(treat / ref).

Outputs under ``--out-dir`` (written on every code path, including
zero-credible-effect runs):
  - sccoda_report.html                     summary report (this run's target)
  - da/{prefix}_da.csv                     credible effects + final_prob
  - da/{prefix}_effects_barplot.png        per-cell-type log2FC bar
  - da/{prefix}_stacked_barplot.png        cell-type proportions per sample

where ``prefix`` = ``{tissue}_{label}`` and ``label`` is the ``=label`` part
of the compare entry (default ``{ref}_vs_{treat}``).

Usage:
  python sccoda.py \\
    --h5ad ovaries=/path/ovaries.h5ad \\
    --compare ovaries:luanchao-21310:luanchao-11238=Aging_vs_Youth \\
    --out-dir ./da_sccoda

  # explicit reference cell type and stricter credible threshold:
  python sccoda.py \\
    --h5ad uterus=/path/uterus.h5ad \\
    --compare uterus:zigong-21310:ZIGONG-21224=Sham_vs_OV-POI \\
    --out-dir ./da_sccoda --ref-celltype Fibroblast --credible-threshold 0.95

--h5ad   tissue=/path/to.h5ad         (repeatable)
--compare tissue:ref:treat[=label]    (repeatable, label optional)
"""
import argparse
import logging
import os
import re
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

warnings.filterwarnings("ignore")
os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import pertpy as pt
import pertpy.tools  # noqa: F401  -- ensure jax deps register


def setup_logger(logger_name: str, level: int = logging.INFO,
                 log_file: Optional[str] = None) -> logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.

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
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


logger = setup_logger(__name__)


# ---------------------------------------------------------------------------
# CLI parsing / validation
# ---------------------------------------------------------------------------

def parse_compare_items(items: List[str],
                        tissues: Dict[str, Path]) -> List[dict]:
    """Parse ``tissue:ref:treat[=label]`` compare entries.

    Parameters
    ----------
    items : list of str
        Raw ``--compare`` values.
    tissues : dict
        Mapping of tissue name to h5ad path from ``--h5ad``.

    Returns
    -------
    list of dict
        One dict per comparison: ``{tissue, ref, treat, label, key}``.

    Raises
    ------
    SystemExit
        On malformed entries, unknown tissues, or ref == treat.
    """
    comparisons = []
    seen = set()
    for item in items:
        main, sep, label = item.partition("=")
        label = label.strip() if sep else ""
        parts = main.split(":")
        if len(parts) != 3 or not all(parts):
            raise SystemExit(
                f"--compare format: tissue:ref:treat[=label], got: {item}")
        tissue, ref, treat = (p.strip() for p in parts)
        if tissue not in tissues:
            raise SystemExit(
                f"Unknown tissue '{tissue}' in --compare '{item}', "
                f"define --h5ad first")
        if ref == treat:
            raise SystemExit(
                f"--compare '{item}': ref and treat are the same sample")
        key = f"{ref}_vs_{treat}"
        if not label:
            label = key
        dup = (tissue, ref, treat)
        if dup in seen:
            raise SystemExit(f"Duplicate --compare entry: {item}")
        seen.add(dup)
        comparisons.append(dict(tissue=tissue, ref=ref, treat=treat,
                                label=label, key=key))
    return comparisons


# ---------------------------------------------------------------------------
# scCODA per comparison
# ---------------------------------------------------------------------------

def run_sccoda(adata: sc.AnnData, comp: dict, formula: str,
               ref_celltype: str, credible_threshold: float,
               num_samples: int, num_warmup: int,
               da_dir: Path) -> Tuple[pd.DataFrame, List[str]]:
    """Run scCODA for one ref-vs-treat comparison and write its artifacts.

    Parameters
    ----------
    adata : AnnData
        Object with ``obs['sample_id']`` and ``obs['cell_type']``.
    comp : dict
        ``{tissue, ref, treat, label, key}`` from :func:`parse_compare_items`.
    formula : str
        Design formula passed to ``sccoda.prepare`` (e.g. ``condition`` or
        ``condition + batch``). Terms other than ``condition`` must exist in
        ``obs``.
    ref_celltype : str
        ``automatic`` or an explicit cell-type name.
    credible_threshold : float
        Posterior probability cutoff for ``is_credible``.
    num_samples, num_warmup : int
        NUTS sampler budgets.
    da_dir : Path
        Directory for the per-comparison CSV and plots (created if missing).

    Returns
    -------
    (DataFrame, list of str)
        Per-cell-type effect table and human-readable notes (skipped plots
        etc.) for the report.
    """
    ref, treat = comp["ref"], comp["treat"]
    tissue, label = comp["tissue"], comp["label"]
    notes: List[str] = []
    prefix = f"{tissue}_{label}"
    da_dir.mkdir(parents=True, exist_ok=True)

    sub = adata[adata.obs["sample_id"].isin([ref, treat])].copy()
    for sid in (ref, treat):
        n = int((sub.obs["sample_id"] == sid).sum())
        if n == 0:
            raise RuntimeError(
                f"[{tissue}] comparison {comp['key']}: sample '{sid}' has "
                f"0 cells in the h5ad — check the sample_id annotation")
        logger.info("  %s: %d cells", sid, n)
    # Order matters for scCODA: the FIRST level of `condition` becomes the
    # base category (no [T.] marker in covariate names); the SECOND level
    # becomes the treatment ([T.<level>]). We force `categories=[ref, treat]`
    # so effects are log2(treat / ref) — "reference first, treatment second".
    sub.obs["condition"] = pd.Categorical(
        sub.obs["sample_id"].values, categories=[ref, treat])

    covariate_terms = [t for t in re.split(r"\W+", formula) if t and t != "1"]
    if "condition" not in covariate_terms:
        raise RuntimeError(
            f"--formula '{formula}' must include 'condition' (the modeled "
            f"ref/treat column this script builds)")
    for term in covariate_terms:
        if term != "condition" and term not in sub.obs.columns:
            raise RuntimeError(
                f"--formula term '{term}' not found in h5ad obs columns "
                f"({sorted(sub.obs.columns)})")

    sccoda = pt.tl.Sccoda()
    logger.info("Aggregating cell-level to sample-level (cell_type x sample)")
    mdata = sccoda.load(
        sub,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier="cell_type",
        sample_identifier="sample_id",
        covariate_obs=covariate_terms,
    )
    logger.info("  sample-level shape: %s x %s",
                mdata["coda"].shape[0], mdata["coda"].shape[1])

    sccoda.prepare(mdata, formula=formula,
                   reference_cell_type=ref_celltype)
    logger.info("Running scCODA NUTS sampler (num_samples=%d, num_warmup=%d)",
                num_samples, num_warmup)
    sccoda.run_nuts(mdata, num_samples=num_samples, num_warmup=num_warmup,
                    rng_key=0)

    summary = sccoda.get_effect_df(mdata)
    summary = summary.reset_index().rename(columns={
        "Cell Type": "cell_type",
        "log2-fold change": "log2FC",
        "Inclusion probability": "final_prob",
    })
    summary["final_prob"] = summary["final_prob"].astype(float)
    summary["log2FC"] = summary["log2FC"].astype(float)
    summary["is_credible"] = summary["final_prob"] >= credible_threshold

    n_up = int(((summary["log2FC"] > 0) & summary["is_credible"]).sum())
    n_dn = int(((summary["log2FC"] < 0) & summary["is_credible"]).sum())
    n_total = int(summary["is_credible"].sum())
    ref_used = mdata["coda"].uns["scCODA_params"]["reference_cell_type"]
    logger.info("  credible: %d (%d up / %d dn), reference cell type=%s",
                n_total, n_up, n_dn, ref_used)
    notes.append(f"scCODA reference cell type: {ref_used}")
    if n_total:
        logger.info("  credible cell types:\n%s",
                    summary[summary["is_credible"]]
                    [[ "cell_type", "log2FC", "final_prob"]]
                    .to_string(index=False))

    # CSV is the primary result — always written, empty results keep the
    # real header row.
    csv_path = da_dir / f"{prefix}_da.csv"
    summary.to_csv(csv_path, index=False)
    logger.info("Saved: %s", csv_path)

    title = f"{tissue.upper()} DA: {label} ({ref} vs {treat})"

    if n_total > 0:
        try:
            fg = sccoda.plot_effects_barplot(
                mdata,
                covariates="condition",
                parameter="log2-fold change",
                plot_facets=True,
                plot_zero_covariate=False,
                plot_zero_cell_type=False,
                figsize=(10, 4),
                return_fig=True,
            )
            if hasattr(fg, "fig"):
                # Force the actual canvas size — pertpy leaves a tall FacetGrid
                # even when there's only one facet, so resize manually.
                fg.fig.set_size_inches(10, 5.5)
                fg.fig.suptitle(title, fontsize=14, y=0.97)
                for ax in fg.axes.flat:
                    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                             fontsize=10)
                    ax.set_title("")
                    ax.set_ylabel(ax.get_ylabel(), fontsize=12)
                fg.fig.subplots_adjust(left=0.08, right=0.98, top=0.85,
                                       bottom=0.38)
                out = da_dir / f"{prefix}_effects_barplot.png"
                fg.fig.savefig(out, dpi=300, facecolor="white")
                plt.close(fg.fig)
                logger.info("Saved: %s", out)
        except Exception as e:
            plt.close("all")
            notes.append(f"effects_barplot skipped (plot error: {e})")
            logger.warning("plot_effects_barplot failed: %s", e)
    else:
        notes.append("effects_barplot skipped: no credible effects")
        logger.info("  skip effects_barplot: no credible effects to plot")

    try:
        fig = sccoda.plot_stacked_barplot(
            mdata, feature_name="samples", figsize=(10, 4.5),
            return_fig=True,
        )
        if fig is None:
            # pertpy sometimes still leaves the figure open — grab current
            fig = plt.gcf()
        fig.set_size_inches(10, 5.0)
        fig.suptitle(title, fontsize=14, y=0.97)
        fig.subplots_adjust(left=0.08, right=0.78, top=0.86, bottom=0.22)
        out = da_dir / f"{prefix}_stacked_barplot.png"
        fig.savefig(out, dpi=300, facecolor="white")
        plt.close(fig)
        logger.info("Saved: %s", out)
    except Exception as e:
        plt.close("all")
        notes.append(f"stacked_barplot skipped (plot error: {e})")
        logger.warning("plot_stacked_barplot failed: %s", e)

    return summary, notes


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(out_dir: Path, comparisons: List[dict],
                 results: Dict[str, dict], run_params: dict) -> None:
    """Write ``<out_dir>/sccoda_report.html``: the run's declared output.

    Emitted unconditionally at the end of EVERY code path (including
    zero-comparison and zero-credible-effect runs), built by walking
    *out_dir*, so the declared output can never go missing.

    Parameters
    ----------
    out_dir : Path
        Output root (the report is written directly into it).
    comparisons : list of dict
        Comparison dicts from :func:`parse_compare_items`.
    results : dict
        ``{key: {'summary': DataFrame, 'notes': [str]}}`` per comparison.
    run_params : dict
        CLI parameters echoed into the report header.
    """
    from html import escape as _escape

    report_path = out_dir / "sccoda_report.html"

    # File inventory: every artifact the run produced (report excluded).
    inventory = []
    for root, _dirs, files in os.walk(out_dir):
        for name in sorted(files):
            if name == report_path.name:
                continue
            path = Path(root) / name
            inventory.append((os.path.relpath(path, out_dir),
                              path.stat().st_size))
    inventory.sort()

    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>scCODA differential abundance report</title>",
        "<style>body{font-family:sans-serif;margin:2em;}"
        "table{border-collapse:collapse;}"
        "td,th{border:1px solid #ccc;padding:4px 8px;}"
        "img{max-width:640px;display:block;margin:6px 0;}"
        "code{background:#f4f4f4;padding:2px 4px;}</style>",
        "</head><body>",
        "<h1>scCODA differential abundance report</h1>",
        f"<p>Output root: <code>{_escape(str(out_dir))}</code></p>",
        "<h2>Run parameters</h2><table>",
    ]
    for k, v in run_params.items():
        parts.append(f"<tr><th>{_escape(str(k))}</th>"
                     f"<td>{_escape(str(v))}</td></tr>")
    parts.append("</table>")

    if not comparisons:
        parts.append("<p>No comparisons were configured for this run.</p>")

    for comp in comparisons:
        key = comp["key"]
        res = results.get(key)
        parts.append(f"<h2>{_escape(comp['tissue'])}: "
                     f"{_escape(comp['ref'])} vs {_escape(comp['treat'])} "
                     f"({_escape(comp['label'])})</h2>")
        if res is None:
            parts.append("<p>Comparison did not complete.</p>")
            continue
        summary: pd.DataFrame = res["summary"]
        n_cred = int(summary["is_credible"].sum()) if len(summary) else 0
        parts.append(f"<p>{len(summary)} cell types tested, "
                     f"{n_cred} credible effects "
                     f"(final_prob &ge; {run_params['credible_threshold']})"
                     f"</p>")
        if len(summary):
            parts.append("<table><tr><th>cell_type</th><th>log2FC</th>"
                         "<th>final_prob</th><th>is_credible</th></tr>")
            for _, row in summary.iterrows():
                parts.append(
                    f"<tr><td>{_escape(str(row['cell_type']))}</td>"
                    f"<td>{row['log2FC']:.4f}</td>"
                    f"<td>{row['final_prob']:.4f}</td>"
                    f"<td>{'yes' if row['is_credible'] else 'no'}</td></tr>")
            parts.append("</table>")
        for note in res.get("notes", []):
            parts.append(f"<p>Note: {_escape(note)}</p>")
        for fname in sorted((out_dir / "da").glob(
                f"{comp['tissue']}_{comp['label']}_*.png")):
            rel = os.path.relpath(fname, out_dir)
            parts.append(f"<img src='{_escape(rel)}' "
                         f"alt='{_escape(fname.name)}'>")

    parts.append("<h2>Output files</h2><table>"
                 "<tr><th>file</th><th>bytes</th></tr>")
    for rel, size in inventory:
        parts.append(f"<tr><td><code>{_escape(rel)}</code></td>"
                     f"<td>{size}</td></tr>")
    parts.append("</table>")
    parts.append("</div></body></html>")

    tmp = str(report_path) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write("\n".join(parts))
    os.replace(tmp, report_path)
    logger.info("Report: %s", report_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "scCODA-based differential abundance analysis.\n"
            "\n"
            "PARAMETER CONVENTION (READ FIRST):\n"
            "  --compare tissue:ref:treat   <- ref FIRST, treat SECOND\n"
            "  In every output:\n"
            "    - File names: {tissue}_{ref}_vs_{treat}_*  (label overrides)\n"
            "    - Plot titles: 'ref vs treat'  (ref on the LEFT)\n"
            "  In the CSV covariate column, treat is tagged with [T.]:\n"
            "  [T.treat] so log2FC = log2(treat / ref).  ref = baseline /\n"
            "  control / unperturbed; treat = experimental / treated /\n"
            "  perturbed.  Convention: 'reference first, treatment second'."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Ovaries - youth (ref) vs aged (treat):
  python sccoda.py \\
    --h5ad ovaries=/path/ovaries.h5ad \\
    --compare ovaries:luanchao-21310:luanchao-11238=Aging_vs_Youth \\
    --out-dir ./da_sccoda

  # Uterus - sham (ref) vs POI (treat):
  python sccoda.py \\
    --h5ad uterus=/path/uterus.h5ad \\
    --compare uterus:zigong-21310:ZIGONG-21224=Sham_vs_OV-POI \\
    --out-dir ./da_sccoda
""")
    p.add_argument("-a", "--h5ad", action="append", required=True,
                   metavar="tissue=path",
                   help="Tissue name = h5ad path (repeatable)")
    p.add_argument("-c", "--compare", action="append", required=True,
                   metavar="tissue:ref:treat[=label]",
                   help="Comparison ref vs treat (repeatable, label optional)")
    p.add_argument("-o", "--out-dir", required=True,
                   help="Output directory (report + da/ artifacts)")
    p.add_argument("-r", "--ref-celltype", default="automatic",
                   help="Reference cell type ('automatic' or explicit name)")
    p.add_argument("-f", "--formula", default="condition",
                   help="Design formula (default: condition)")
    p.add_argument("-T", "--credible-threshold", type=float, default=0.90,
                   help="Posterior probability threshold for credible "
                        "effect (default: 0.90)")
    p.add_argument("-n", "--num-samples", type=int, default=2000,
                   help="NUTS posterior samples (default: 2000)")
    p.add_argument("-w", "--num-warmup", type=int, default=500,
                   help="NUTS warmup samples (default: 500)")
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tissues: Dict[str, Path] = {}
    for item in args.h5ad:
        if "=" not in item:
            raise SystemExit(f"--h5ad format: tissue=path, got: {item}")
        tissue, path = item.split("=", 1)
        tissues[tissue.strip()] = Path(path.strip())

    comparisons = parse_compare_items(args.compare, tissues)
    run_params = {
        "formula": args.formula,
        "ref_celltype": args.ref_celltype,
        "credible_threshold": args.credible_threshold,
        "num_samples": args.num_samples,
        "num_warmup": args.num_warmup,
    }

    results: Dict[str, dict] = {}
    loaded: Dict[str, sc.AnnData] = {}
    for comp in comparisons:
        tissue = comp["tissue"]
        logger.info("=" * 60)
        logger.info("[%s] %s vs %s (%s)", tissue, comp["ref"],
                    comp["treat"], comp["label"])
        logger.info("=" * 60)
        if tissue not in loaded:
            adata = sc.read_h5ad(tissues[tissue])
            missing = [c for c in ("sample_id", "cell_type")
                       if c not in adata.obs.columns]
            if missing:
                raise RuntimeError(
                    f"{tissues[tissue]}: obs missing required column(s) "
                    f"{missing} — expected scanpy_auto-style annotations")
            logger.info("Loaded %s: shape=%s", tissue, adata.shape)
            loaded[tissue] = adata
        summary, notes = run_sccoda(
            loaded[tissue], comp,
            formula=args.formula,
            ref_celltype=args.ref_celltype,
            credible_threshold=args.credible_threshold,
            num_samples=args.num_samples,
            num_warmup=args.num_warmup,
            da_dir=out_dir / "da")
        results[comp["key"]] = {"summary": summary, "notes": notes}

    write_report(out_dir, comparisons, results, run_params)
    logger.info("DONE")


if __name__ == "__main__":
    main()
