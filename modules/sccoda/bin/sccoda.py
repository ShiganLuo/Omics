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
from matplotlib.patches import Rectangle
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

def _plot_stacked_ribbon(counts_df: pd.DataFrame, out_path: Path,
                         title: str) -> None:
    """Stacked proportion bars per sample with light same-color ribbons.

    Reference da_sccoda style with an alluvial touch: bars are stacked cell-
    type proportions per sample, and same-cell-type segments of adjacent bars
    are connected by a lighter fill of that color. No x-axis label (sample
    names on the ticks are self-explanatory), two-line title, tight bbox.

    Parameters
    ----------
    counts_df : DataFrame
        Samples x cell types raw counts (index = sample ids in plot order).
    out_path : Path
        Target PNG.
    title : str
        Plot title; the `` (`` before the comparison is turned into a line
        break so long titles cannot clip at the canvas edge.
    """
    counts = counts_df.to_numpy(dtype=np.float64)
    samples = [str(s) for s in counts_df.index]
    cell_types = [str(c) for c in counts_df.columns]
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    props = counts / totals * 100.0

    # Stable stacking order: most abundant cell type at the bottom.
    order = np.argsort(-counts.sum(axis=0))
    props = props[:, order]
    cell_types = [cell_types[i] for i in order]

    n_s, n_c = props.shape
    cmap = plt.get_cmap("tab20" if n_c <= 20 else "tab20b")
    colors = [cmap(i / max(n_c - 1, 1)) for i in range(n_c)]

    fig, ax = plt.subplots(figsize=(max(6.0, 1.2 * n_s + 4.5), 6))
    width = 0.55
    xs = np.arange(n_s, dtype=float)

    # Segment bounds per sample (n_s x n_c), bottom to top.
    y1 = np.cumsum(props, axis=1)
    y0 = y1 - props

    # Ribbons first (behind the bars): lighter fill of the same color
    # connecting same-cell-type segments between adjacent bars. Vertical
    # transition uses a smoothstep (S-curve) so ribbons read as soft flows
    # instead of angular trapezoids.
    t = np.linspace(0.0, 1.0, 40)
    s = t * t * (3.0 - 2.0 * t)  # smoothstep
    for i in range(n_c):
        for j in range(n_s - 1):
            xl0, xl1 = xs[j] + width / 2, xs[j + 1] - width / 2
            xline = xl0 + (xl1 - xl0) * t
            y0line = y0[j, i] + (y0[j + 1, i] - y0[j, i]) * s
            y1line = y1[j, i] + (y1[j + 1, i] - y1[j, i]) * s
            ax.fill(np.concatenate([xline, xline[::-1]]),
                    np.concatenate([y0line, y1line[::-1]]),
                    color=colors[i], alpha=0.22, linewidth=0, zorder=1)

    # Bars on top.
    for i in range(n_c):
        ax.bar(xs, props[:, i], bottom=y0[:, i], width=width,
               color=colors[i], edgecolor="white", linewidth=0.4, zorder=2)

    ax.set_xticks(xs)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Proportion (%)", fontsize=12)
    ax.set_ylim(0, 100)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = [Rectangle((0, 0), 1, 1, color=colors[i])
               for i in range(n_c)]
    ax.legend(handles, cell_types, bbox_to_anchor=(1.02, 1.0),
              loc="upper left", fontsize=9, frameon=False, title="cell_type")
    # Title OUTSIDE the axes (pad-based): a y=0.97 multi-line title dips
    # into the plot and overlaps the top y tick label.
    ax.set_title(title.replace(" (", "\n("), fontsize=14, pad=16)
    fig.subplots_adjust(top=0.88)
    fig.savefig(out_path, dpi=300, facecolor="white",
                bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)


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
    # Reference naming (analysis/results/da_sccoda convention):
    # da_sccoda_<ref>_vs_<treat>_{da.csv,effects_barplot.png,stacked_barplot.png}
    prefix = f"da_sccoda_{comp['ref']}_vs_{comp['treat']}"
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

    # --- numpyro/anndata compatibility shim (pertpy 1.0.x quirks) ----------
    # 1) anndata>=0.12 coerces 1-D obsm entries into (N,1) columns, so
    #    prepare()'s sample_counts (= np.sum(X, 1), meant to be (N,)) is
    #    stored and passed on as (N,1). numpyro's DirichletMultinomial
    #    broadcasts that into an (N,N) batch whose multinomial support check
    #    marks every obs row out-of-support -> -inf log-probs -> "Cannot find
    #    valid initial parameters" (or "invalid total_count parameter" when
    #    fractional). obsm assignment cannot hold a 1-D array, so flatten at
    #    the model boundary instead.
    # 2) prepare() pseudocounts zeros with 0.5, yielding fractional counts
    #    that numpyro's nonnegative_integer total_count check rejects. Pre-
    #    fill zeros with 1 (integer) before prepare so its 0.5 path never
    #    fires and counts/totals stay integer-valued.
    _orig_model = sccoda.model

    def _model(counts, covariates, n_total, ref_index, sample_adata):
        return _orig_model(counts, covariates,
                           np.asarray(n_total, dtype=np.float64).reshape(-1),
                           ref_index, sample_adata)

    sccoda.model = _model

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

    # Zero pre-fill for the model shim above (integer pseudocount).
    coda = mdata["coda"]
    X = np.asarray(coda.X.todense() if hasattr(coda.X, "todense") else coda.X,
                   dtype=np.float64)
    n_zero = int((X == 0).sum())
    if n_zero:
        logger.info("  %d zero sample-level counts; pre-filling with 1 "
                    "(integer pseudocount)", n_zero)
        X[X == 0] = 1.0
    coda.X = X

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

    # DA effect plot — reference style (analysis/results/da_sccoda): one
    # horizontal bar per cell type, sorted by log2FC, credible effects in
    # distinct colors with black edges, non-credible grey.
    if len(summary):
        try:
            dfp = summary.sort_values("log2FC").reset_index(drop=True)
            n = len(dfp)
            n_credible = int(dfp["is_credible"].sum())
            cmap = plt.get_cmap("tab20" if n_credible <= 20 else "tab20b")
            cred_colors = [cmap(i / max(n_credible - 1, 1))
                           for i in range(n_credible)]
            ci = 0
            colors = []
            for _, row in dfp.iterrows():
                if row["is_credible"]:
                    colors.append(cred_colors[ci])
                    ci += 1
                else:
                    colors.append("#cccccc")
            fig, ax = plt.subplots(figsize=(10, max(4, n * 0.45 + 1.5)))
            for i, (_, row) in enumerate(dfp.iterrows()):
                edgecolor = "black" if row["is_credible"] else "none"
                lw = 0.8 if row["is_credible"] else 0
                ax.barh(i, row["log2FC"], color=colors[i],
                        edgecolor=edgecolor, linewidth=lw, height=0.7,
                        zorder=2)
            ax.set_yticks(range(n))
            ax.set_yticklabels(dfp["cell_type"], fontsize=10)
            ax.set_xlabel("log2-fold change", fontsize=12)
            # Two-line title + tight bbox: long "<label> (<ref> vs <treat>)"
            # titles clipped at the figure's right edge otherwise.
            ax.set_title(title.replace(" (", "\n("), fontsize=14, y=0.97)
            ax.axvline(0, color="black", linewidth=0.5, alpha=0.5, zorder=1)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(axis="x", labelsize=10)
            fig.subplots_adjust(left=0.25, right=0.98, top=0.88, bottom=0.10)
            out = da_dir / f"{prefix}_effects_barplot.png"
            fig.savefig(out, dpi=300, facecolor="white",
                        bbox_inches="tight", pad_inches=0.3)
            plt.close(fig)
            logger.info("Saved: %s", out)
        except Exception as e:
            plt.close("all")
            notes.append(f"effects_barplot skipped (plot error: {e})")
            logger.warning("effects_barplot failed: %s", e)
    else:
        notes.append("effects_barplot skipped: empty result table")
        logger.info("  skip effects_barplot: empty result table")

    try:
        # Stacked proportions per sample with same-color ribbon links.
        # Built from the raw cell-level crosstab (the coda counts carry the
        # integer zero pre-fill, which would distort proportions).
        counts_df = pd.crosstab(sub.obs["sample_id"], sub.obs["cell_type"])
        counts_df = counts_df.reindex(
            index=[ref, treat], columns=mdata["coda"].var.index, fill_value=0)
        out = da_dir / f"{prefix}_stacked_barplot.png"
        _plot_stacked_ribbon(counts_df, out, title)
        logger.info("Saved: %s", out)
    except Exception as e:
        plt.close("all")
        notes.append(f"stacked_barplot skipped (plot error: {e})")
        logger.warning("stacked_barplot failed: %s", e)

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
                f"da_sccoda_{comp['ref']}_vs_{comp['treat']}_*.png")):
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
