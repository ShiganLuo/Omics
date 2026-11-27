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
(e.g. the scanpy_auto output). The comparison design is JSON
(``--compare-json``), mirroring the DESeq2 ``group_pairs`` shape::

    {"<tissue>": {"<contrast>": {
        "control_group_name": "youth",
        "experimental_group_name": "aged",
        "control_samples": ["luanchao-21310-10XSC3"],
        "experimental_samples": ["luanchao-11238-10XSC3"]}}}

A multi-group contrast (several consecutive bars in ONE stacked
barplot) uses an ordered group list instead; the FIRST group is the
scCODA base level and every later group is compared against it::

    {"<tissue>": {"<contrast>": {"groups": [
        {"name": "youth", "samples": ["s1", "s2"]},
        {"name": "aged", "samples": ["s3", "s4"]}]}}}

Effects are log2FC = log2(group / control). Multi-sample groups act as
model replicates and are summed into ONE bar of the stacked barplot.
The stacked barplot renders DA status with ↑/↓ arrows directly after
each cell-type legend name — one glyph per comparison, in that
comparison's color (the color, not the position, identifies the
comparison), centered on the name's row. Below, a DA key legend maps
each color to its meaning: the ``↑↓`` glyph pair as the handle
(leftmost, left-aligned with the cell-type swatches) and the reading as
the label (``ctrl vs exp`` + ↑/↓ semantics — direction phrased relative
to the EXPERIMENTAL group: ↑ = higher in exp, ↓ = lower in exp).

Outputs under ``--out-dir`` (written on every code path, including
zero-credible-effect runs):
  - sccoda_report.html                     summary report (this run's target)
  - da/{prefix}_da.csv                     credible effects + final_prob
  - da/{prefix}_effects_barplot.png        per-cell-type log2FC bars
  - da/{prefix}_stacked_barplot.png        cell-type proportions per group

where ``prefix`` = ``da_sccoda_{contrast}`` (the design's contrast key).

Usage:
  python sccoda.py \\
    --h5ad ovaries=/path/ovaries.h5ad \\
    --compare-json design.json \\
    --out-dir ./da_sccoda

  # explicit reference cell type and stricter credible threshold:
  python sccoda.py \\
    --h5ad uterus=/path/uterus.h5ad \\
    --compare-json design.json \\
    --out-dir ./da_sccoda --ref-celltype Fibroblast --credible-threshold 0.95

--h5ad   tissue=/path/to.h5ad         (repeatable)
--compare-json FILE                   (repeatable; design shape above)
"""
import argparse
import json
import logging
import os
import re
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

warnings.filterwarnings("ignore")
os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Rectangle
from matplotlib.text import Text
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

PlotFormat = str

VALID_FORMATS = ("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg",
                 "jpeg", "pgf", "raw", "rgba")


def save_fig(fig, out_base: Union[str, Path],
             image_formats: Optional[List[PlotFormat]] = None,
             dpi: int = 300) -> None:
    """Save *fig* as ``<out_base>.<fmt>`` for every requested format.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save (closed afterwards).
    out_base : str or Path
        Output path without extension.
    image_formats : list of str, optional
        Save formats. Default ``["png"]``.
    dpi : int, default=300
        Output resolution (raster formats; vector formats ignore it).
    """
    if image_formats is None:
        image_formats = ["png"]
    bad = [f for f in image_formats if f not in VALID_FORMATS]
    if bad:
        raise ValueError(
            f"Unsupported image format(s): {', '.join(bad)}; "
            f"valid: {', '.join(VALID_FORMATS)}")
    for fmt in image_formats:
        out_path = f"{out_base}.{fmt}"
        fig.savefig(out_path, dpi=dpi, facecolor="white",
                    bbox_inches="tight", pad_inches=0.3)
        logger.info("Saved: %s", out_path)
    plt.close(fig)



def _parse_group_entry(tissue: str, key: str, entry: dict
                       ) -> Tuple[List[Tuple[str, List[str]]], str]:
    """Normalize one design entry into ``(ordered groups, base group)``.

    Two accepted shapes:

    - pairwise (DESeq2 ``group_pairs`` keys)::

        {"control_group_name": "youth",
         "experimental_group_name": "aged",
         "control_samples": ["s1"], "experimental_samples": ["s2"]}

    - multi-group (ordered; bars follow the listed order)::

        {"groups": [{"name": "youth", "samples": ["s1", "s2"]}, ...],
         "control_group_name": "youth"}   # optional, default = first

    ``control_group_name`` names the scCODA base level (every other group
    is compared against it); it need not be the first group — the listed
    order only sets the stacked-barplot bar order.

    Parameters
    ----------
    tissue, key : str
        Location of the entry in the design (for error messages).
    entry : dict
        The raw design entry.

    Returns
    -------
    (list of (str, list of str), str)
        Group names and sample ids in plot order, and the base group name.

    Raises
    ------
    SystemExit
        On unknown shapes, missing/empty fields, duplicate group names,
        an unknown control group, or a sample claimed by two groups.
    """
    where = f"compare design entry '{tissue}:{key}'"
    if not isinstance(entry, dict):
        raise SystemExit(f"{where}: entry must be an object, got {type(entry).__name__}")
    if "groups" in entry:
        raw_groups = entry.get("groups") or []
        if not isinstance(raw_groups, list) or len(raw_groups) < 2:
            raise SystemExit(
                f"{where}: 'groups' must be a list of at least 2 "
                f"{{'name', 'samples'}} objects")
        groups: List[Tuple[str, List[str]]] = []
        for g in raw_groups:
            if not isinstance(g, dict) or not g.get("name") or not g.get("samples"):
                raise SystemExit(
                    f"{where}: each group needs a non-empty 'name' and "
                    f"'samples' list, got {g!r}")
            name = str(g["name"]).strip()
            samples = [str(s).strip() for s in g["samples"]]
            groups.append((name, samples))
        control = str(entry.get("control_group_name") or groups[0][0]).strip()
    else:
        required = ("control_group_name", "experimental_group_name",
                    "control_samples", "experimental_samples")
        missing = [k for k in required if k not in entry]
        if missing:
            raise SystemExit(
                f"{where}: missing {missing}; expected the DESeq2 "
                f"group_pairs keys {list(required)} (or an ordered "
                f"'groups' list for a multi-group contrast)")
        groups = [
            (str(entry["control_group_name"]).strip(),
             [str(s).strip() for s in entry["control_samples"]]),
            (str(entry["experimental_group_name"]).strip(),
             [str(s).strip() for s in entry["experimental_samples"]]),
        ]
        control = groups[0][0]
    names = [n for n, _ in groups]
    if any(not n for n in names):
        raise SystemExit(f"{where}: empty group name")
    if len(set(names)) != len(names):
        raise SystemExit(f"{where}: duplicate group name in {names}")
    if control not in names:
        raise SystemExit(
            f"{where}: control_group_name '{control}' is not one of {names}")
    for name, samples in groups:
        if not samples or any(not s for s in samples):
            raise SystemExit(f"{where}: group '{name}' needs non-empty sample ids")
        if len(set(samples)) != len(samples):
            raise SystemExit(f"{where}: duplicate sample in group '{name}': {samples}")
    flat = [s for _, samples in groups for s in samples]
    dupes = sorted({s for s in flat if flat.count(s) > 1})
    if dupes:
        raise SystemExit(f"{where}: samples claimed by two groups: {dupes}")
    return groups, control


def _covariate_group(cov: str, group_names: List[str]) -> str:
    """Map a scCODA covariate row label back to its compared group name.

    pertpy labels the dummy covariates ``conditionT.<group>`` (e.g.
    ``conditionT.aged``); every non-base group gets one covariate.
    Longest group names are matched first so one name cannot shadow a
    suffix of another. Unmatched labels (extra formula terms) map to
    themselves.
    """
    for name in sorted(group_names[1:], key=len, reverse=True):
        if cov == f"conditionT.{name}" or cov.endswith(f"T.{name}"):
            return name
    return cov


def parse_compare_json(items: List[str],
                       tissues: Dict[str, Path]) -> List[dict]:
    """Parse ``--compare-json`` design files into comparison dicts.

    Design shape (one file, repeatable and merged)::

        {"<tissue>": {"<contrast>": <entry>}}

    where ``<entry>`` is a pairwise ``group_pairs`` object or a
    multi-group ``groups`` list (see :func:`_parse_group_entry`). The
    contrast key names the outputs (``da_sccoda_<contrast>_*``).

    Parameters
    ----------
    items : list of str
        Paths from ``--compare-json``.
    tissues : dict
        Mapping of tissue name to h5ad path from ``--h5ad``.

    Returns
    -------
    list of dict
        One dict per comparison: ``{tissue, key, label, groups}`` where
        ``groups`` is an ordered ``[(name, [samples]), ...]`` (FIRST =
        base level).

    Raises
    ------
    SystemExit
        On unreadable files, wrong shapes, unknown tissues, duplicate
        contrasts, or invalid group/sample specs.
    """
    design: Dict[str, Dict[str, dict]] = {}
    for path in items:
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError) as e:
            raise SystemExit(f"--compare-json {path}: cannot read design: {e}")
        if not isinstance(data, dict):
            raise SystemExit(
                f"--compare-json {path}: top level must be an object of "
                f"{{tissue: {{contrast: entry}}}}")
        for tissue, contrasts in data.items():
            if not isinstance(contrasts, dict) or not contrasts:
                raise SystemExit(
                    f"--compare-json {path}: tissue '{tissue}' must map to a "
                    f"non-empty object of contrast entries")
            tgt = design.setdefault(str(tissue), {})
            for key, entry in contrasts.items():
                if str(key) in tgt:
                    raise SystemExit(
                        f"--compare-json {path}: duplicate contrast '{key}' "
                        f"for tissue '{tissue}'")
                tgt[str(key)] = entry
    if not design:
        raise SystemExit("--compare-json: design is empty")

    comparisons = []
    for tissue, contrasts in design.items():
        if tissue not in tissues:
            raise SystemExit(
                f"Unknown tissue '{tissue}' in --compare-json design, "
                f"define --h5ad first")
        for key, entry in contrasts.items():
            groups, control = _parse_group_entry(tissue, key, entry)
            comparisons.append(dict(tissue=tissue, key=key, label=key,
                                    groups=groups, control=control))
    return comparisons



def _plot_stacked_ribbon(counts_df: pd.DataFrame, out_base: Path,
                         image_formats: Optional[List[PlotFormat]] = None,
                         cell_arrows: Optional[Dict[str, List[Tuple[str, Any]]]]
                         = None,
                         da_key_rows: Optional[List[Tuple[str, Any]]] = None
                         ) -> None:
    """Stacked proportion bars per group with light same-color ribbons.

    Reference da_sccoda style with an alluvial touch: ONE bar per group
    (counts pre-summed across the group's samples), the groups displayed
    as consecutive bars along the x axis with the GROUP NAME as tick
    label, and same-cell-type segments of adjacent bars connected by a
    lighter fill of that color. No title.

    DA status rendering: the ↑/↓ arrows sit DIRECTLY AFTER each
    cell-type legend name, one glyph per comparison in that
    comparison's color (the color identifies the comparison — multi-group
    safe), vertically centered on the name's row so each arrow unambiguously
    belongs to its cell type. Below the cell-type legend a DA key legend
    maps each comparison's color to its reading: one row per comparison
    whose HANDLE is the ``↑↓`` glyph pair in that color (leftmost column,
    left-aligned with the cell-type swatches) and whose label is the
    meaning (``ctrl vs exp`` + ↑/↓ semantics). Placements are measured at
    the save dpi (legend layout is dpi-dependent).

    Parameters
    ----------
    counts_df : DataFrame
        Groups x cell types raw counts (index = group names in plot
        order; a group with several samples must already be summed).
    out_base : Path
        Output path without extension.
    image_formats : list of str, optional
        Save formats. Default ``["png"]``.
    cell_arrows : dict, optional
        ``{cell_type: [(arrow, color), ...]}`` — arrows appended after
        the legend name, in listed order.
    da_key_rows : list of (label, color), optional
        Arrow-legend key rows, one per comparison (bold).
    """
    counts = counts_df.to_numpy(dtype=np.float64)
    groups = [str(g) for g in counts_df.index]
    cell_types = [str(c) for c in counts_df.columns]
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    props = counts / totals * 100.0

    order = np.argsort(-counts.sum(axis=0))
    props = props[:, order]
    cell_types = [cell_types[i] for i in order]

    n_g, n_c = props.shape
    cmap = plt.get_cmap("tab20" if n_c <= 20 else "tab20b")
    colors = [cmap(i / max(n_c - 1, 1)) for i in range(n_c)]

    fig, ax = plt.subplots(figsize=(max(6.0, 1.2 * n_g + 4.5), 6))
    width = 0.55
    xs = np.arange(n_g, dtype=float)

    y1 = np.cumsum(props, axis=1)
    y0 = y1 - props

    t = np.linspace(0.0, 1.0, 40)
    s = t * t * (3.0 - 2.0 * t)
    for i in range(n_c):
        for j in range(n_g - 1):
            xl0, xl1 = xs[j] + width / 2, xs[j + 1] - width / 2
            xline = xl0 + (xl1 - xl0) * t
            y0line = y0[j, i] + (y0[j + 1, i] - y0[j, i]) * s
            y1line = y1[j, i] + (y1[j + 1, i] - y1[j, i]) * s
            ax.fill(np.concatenate([xline, xline[::-1]]),
                    np.concatenate([y0line, y1line[::-1]]),
                    color=colors[i], alpha=0.22, linewidth=0, zorder=1)

    for i in range(n_c):
        ax.bar(xs, props[:, i], bottom=y0[:, i], width=width,
               color=colors[i], edgecolor="white", linewidth=0.4, zorder=2)

    ax.set_xticks(xs)
    rotate = max((len(g) for g in groups), default=0) > 10
    ax.set_xticklabels(groups, rotation=45 if rotate else 0,
                       ha="right" if rotate else "center", fontsize=11)
    ax.set_ylabel("Proportion (%)", fontsize=12)
    ax.set_ylim(0, 100)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = [Rectangle((0, 0), 1, 1, color=colors[i])
               for i in range(n_c)]
    leg_cells = ax.legend(handles, cell_types, bbox_to_anchor=(1.02, 1.0),
                          loc="upper left", fontsize=9, frameon=False,
                          title="cell_type")
    if da_key_rows or cell_arrows:
        leg_da = None
        if da_key_rows:
            ax.add_artist(leg_cells)
            key_handles: Any = [matplotlib.colors.to_rgba(color)
                                for _label, color in da_key_rows]
            leg_da = ax.legend(
                key_handles, [label for label, _ in da_key_rows],
                bbox_to_anchor=(1.02, 1.0), loc="upper left",
                fontsize=9, frameon=False,
                handler_map={tuple: _ArrowPairHandler()})
            for text, (_label, color) in zip(leg_da.get_texts(), da_key_rows):
                text.set_color(color)
                text.set_fontweight("bold")
        fig.set_dpi(300)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        if leg_da is not None:
            bb = leg_cells.get_window_extent()
            (_, y0), _ = ax.transAxes.inverted().transform(bb.get_points())
            leg_da.set_bbox_to_anchor((1.02, y0 - 0.01),
                                      transform=ax.transAxes)
        if cell_arrows:
            inv_fig = fig.transFigure.inverted()
            gap = 0.004
            for name_text, ct in zip(leg_cells.get_texts(), cell_types):
                marks = cell_arrows.get(ct) or []
                if not marks:
                    continue
                bb_name = name_text.get_window_extent(renderer)
                x_end, y_mid = inv_fig.transform(
                    (bb_name.x1, (bb_name.y0 + bb_name.y1) / 2))
                x = x_end + gap
                for arrow, color in marks:
                    at = fig.text(x, y_mid, arrow, color=color,
                                  fontsize=9, va="center", ha="left")
                    bb_arrow = at.get_window_extent(renderer)
                    x = inv_fig.transform((bb_arrow.x1, 0))[0] + gap * 0.5
    save_fig(fig, out_base, image_formats, dpi=300)


class _ArrowPairHandler(HandlerBase):
    """Legend handle drawing the ``↑↓`` glyph pair in the handle color.

    Used by the stacked barplot's DA key legend so the arrow pair sits in
    the LEFTMOST handle column — left-aligned with the cell-type legend's
    color swatches — with the meaning as the row label.
    """

    def legend_artist(self, legend: Any, orig_handle: Any, fontsize: float,
                      handlebox: Any) -> Text:
        """Draw ``↑↓`` centered in the handle box, in ``orig_handle``'s color.

        Parameters
        ----------
        legend : Legend
            Owning legend.
        orig_handle : tuple
            RGBA color of the comparison whose arrows this row explains.
        fontsize : float
            Legend font size in points.
        handlebox : DrawingArea
            Reserved handle column of the legend row.

        Returns
        -------
        matplotlib.text.Text
            The drawn glyph-pair artist.
        """
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        txt = Text(x0 + width / 2, y0 + height / 2, "↑↓",
                   color=orig_handle, fontsize=fontsize,
                   ha="center", va="center")
        handlebox.add_artist(txt)
        return txt


def run_sccoda(adata: sc.AnnData, comp: dict, formula: str,
               ref_celltype: str, credible_threshold: float,
               num_samples: int, num_warmup: int,
               da_dir: Path,
               image_formats: Optional[List[PlotFormat]] = None
               ) -> Tuple[pd.DataFrame, List[str]]:
    """Run scCODA for one ref-vs-treat comparison and write its artifacts.

    Parameters
    ----------
    adata : AnnData
        Object with ``obs['sample_id']`` and ``obs['cell_type']``.
    comp : dict
        ``{tissue, key, label, groups}`` from :func:`parse_compare_json`
        (``groups`` is an ordered ``[(name, [samples]), ...]``; the FIRST
        group is the base level and every later group is compared
        against it; group samples are the model replicates).
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
    image_formats : list of str, optional
        Plot save formats. Default ``["png"]``.

    Returns
    -------
    (DataFrame, list of str)
        Per-cell-type effect table and human-readable notes (skipped plots
        etc.) for the report.
    """
    tissue, label = comp["tissue"], comp["label"]
    groups: List[Tuple[str, List[str]]] = comp["groups"]
    group_names = [name for name, _ in groups]
    base_name = comp["control"]
    compare_names = [n for n in group_names if n != base_name]
    model_names = [base_name] + compare_names
    notes: List[str] = []
    prefix = f"da_sccoda_{comp['key']}"
    da_dir.mkdir(parents=True, exist_ok=True)

    sample_group = {sid: name for name, samples in groups for sid in samples}
    sub = adata[adata.obs["sample_id"].isin(sample_group)].copy()
    for sid in sample_group:
        n = int((sub.obs["sample_id"] == sid).sum())
        if n == 0:
            raise RuntimeError(
                f"[{tissue}] comparison {comp['key']}: sample '{sid}' has "
                f"0 cells in the h5ad — check the sample_id annotation")
        logger.info("  %s (%s): %d cells", sid, sample_group[sid], n)
    sub.obs["condition"] = pd.Categorical(
        sub.obs["sample_id"].map(sample_group), categories=model_names)

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
    if "Covariate" in summary.columns:
        summary["group"] = summary["Covariate"].map(
            lambda c: _covariate_group(str(c), model_names))
    elif len(group_names) == 2:
        summary["group"] = group_names[1]
    else:
        raise RuntimeError(
            f"[{tissue}] comparison {comp['key']}: effect table has no "
            f"'Covariate' column, cannot attribute effects to the "
            f"compared groups {group_names[1:]}")

    ref_used = mdata["coda"].uns["scCODA_params"]["reference_cell_type"]
    logger.info("  reference cell type=%s", ref_used)
    notes.append(f"scCODA reference cell type: {ref_used}")
    for gname in compare_names:
        part = summary[summary["group"] == gname]
        n_up = int(((part["log2FC"] > 0) & part["is_credible"]).sum())
        n_dn = int(((part["log2FC"] < 0) & part["is_credible"]).sum())
        logger.info("  %s vs %s: credible %d (%d up / %d dn)",
                    base_name, gname, n_up + n_dn, n_up, n_dn)
    n_total = int(summary["is_credible"].sum())
    if n_total:
        logger.info("  credible cell types:\n%s",
                    summary[summary["is_credible"]]
                    [["group", "cell_type", "log2FC", "final_prob"]]
                    .to_string(index=False))

    csv_path = da_dir / f"{prefix}_da.csv"
    summary.to_csv(csv_path, index=False)
    logger.info("Saved: %s", csv_path)

    panel_names = list(compare_names)
    if len(summary):
        try:
            all_types = sorted(summary["cell_type"].astype(str).unique())
            cmap = plt.get_cmap("tab20" if len(all_types) <= 20 else "tab20b")
            type_colors = {ct: cmap(i / max(len(all_types) - 1, 1))
                           for i, ct in enumerate(all_types)}
            panel_h = [max(2.5, int((summary["group"] == g).sum()) * 0.45 + 1.6)
                       for g in panel_names]
            fig, axes = plt.subplots(len(panel_names), 1, squeeze=False,
                                     figsize=(10, sum(panel_h)))
            for pi, gname in enumerate(panel_names):
                ax = axes[pi][0]
                dfp = (summary[summary["group"] == gname]
                       .sort_values("log2FC").reset_index(drop=True))
                n = len(dfp)
                if not (dfp["log2FC"].abs() > 0).any():
                    ax.text(0.5, 0.5,
                            "all effects are zero (nothing selected)",
                            transform=ax.transAxes, ha="center", va="center",
                            fontsize=11, color="#999999")
                for i, (_, row) in enumerate(dfp.iterrows()):
                    credible = bool(row["is_credible"])
                    color = type_colors[str(row["cell_type"])] if credible else "#cccccc"
                    ax.barh(i, row["log2FC"], color=color,
                            edgecolor="black" if credible else "none",
                            linewidth=0.8 if credible else 0,
                            height=0.7, zorder=2)
                ylabels = [str(row["cell_type"]) for _, row in dfp.iterrows()]
                ax.set_yticks(range(n))
                ax.set_yticklabels(ylabels, fontsize=10)
                ax.set_xlabel(
                    f"log2-fold change ({base_name} vs {gname})", fontsize=12)
                ax.text(0.0, 1.01, f"← lower in {gname}",
                        transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=9, color="#555555")
                ax.text(1.0, 1.01, f"higher in {gname} →",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=9, color="#555555")
                ax.axvline(0, color="black", linewidth=0.5, alpha=0.5, zorder=1)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.tick_params(axis="x", labelsize=10)
            fig.subplots_adjust(left=0.25, right=0.98)
            out = da_dir / f"{prefix}_effects_barplot"
            save_fig(fig, out, image_formats)
        except Exception as e:
            plt.close("all")
            notes.append(f"effects_barplot skipped (plot error: {e})")
            logger.warning("effects_barplot failed: %s", e)
    else:
        notes.append("effects_barplot skipped: empty result table")
        logger.info("  skip effects_barplot: empty result table")

    try:
        group_of = sub.obs["sample_id"].map(sample_group)
        counts_df = pd.crosstab(group_of, sub.obs["cell_type"])
        counts_df = counts_df.reindex(index=group_names,
                                      columns=mdata["coda"].var.index,
                                      fill_value=0)
        comp_palette = plt.get_cmap("tab10")
        da_key_rows: List[Tuple[str, Any]] = []
        cell_arrows: Dict[str, List[Tuple[str, Any]]] = {}
        for gi, gname in enumerate(compare_names):
            color = comp_palette(gi)
            da_key_rows.append(
                (f"{base_name} vs {gname} (↑/↓ = higher/lower in {gname})",
                 color))
            part = summary[summary["group"] == gname]
            for _, row in part.sort_values("cell_type").iterrows():
                if row["is_credible"] and row["log2FC"] != 0:
                    arrow = "↑" if row["log2FC"] > 0 else "↓"
                    cell_arrows.setdefault(str(row["cell_type"]), []).append(
                        (arrow, color))
        out = da_dir / f"{prefix}_stacked_barplot"
        _plot_stacked_ribbon(
            counts_df, out, image_formats=image_formats,
            cell_arrows=cell_arrows, da_key_rows=da_key_rows)
    except Exception as e:
        plt.close("all")
        notes.append(f"stacked_barplot skipped (plot error: {e})")
        logger.warning("stacked_barplot failed: %s", e)

    return summary, notes



def _figure_variants(out_dir: Path, comp: dict) -> List[Path]:
    """One embeddable figure file per plot base (raster-first).

    Multi-format runs write each da figure once per format; the report
    gallery embeds one variant per base (png preferred, then jpg/jpeg/tif/
    tiff, else the first produced) so every figure appears exactly once.

    Parameters
    ----------
    out_dir : Path
        Run output root (figures live in ``da/``).
    comp : dict
        Comparison dict with ``ref``/``treat`` keys.

    Returns
    -------
    list of Path
        Selected figure files, sorted by base name.
    """
    preferred = ("png", "jpg", "jpeg", "tif", "tiff")
    best: Dict[Path, Tuple[int, Path]] = {}
    for fname in sorted((out_dir / "da").glob(
            f"da_sccoda_{comp['key']}_*")):
        fmt = fname.suffix.lstrip(".").lower()
        if fmt not in VALID_FORMATS:
            continue
        rank = preferred.index(fmt) if fmt in preferred else len(preferred)
        base = fname.with_suffix("")
        cur = best.get(base)
        if cur is None or rank < cur[0]:
            best[base] = (rank, fname)
    return [item[1] for item in sorted(best.values(), key=lambda x: str(x[1]))]


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
        Comparison dicts from :func:`parse_compare_json`.
    results : dict
        ``{key: {'summary': DataFrame, 'notes': [str]}}`` per comparison.
    run_params : dict
        CLI parameters echoed into the report header.
    """
    from html import escape as _escape

    report_path = out_dir / "sccoda_report.html"

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
                     f"{_escape(comp['label'])} "
                     f"({' vs '.join(_escape(n) for n, _ in comp['groups'])})"
                     f"</h2>")
        parts.append("<p>" + "; ".join(
            f"{_escape(name)} ({_escape(' + '.join(samples))})"
            for name, samples in comp["groups"]) + "</p>")
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
        for fname in _figure_variants(out_dir, comp):
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



def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "scCODA-based differential abundance analysis.\n"
            "\n"
            "PARAMETER CONVENTION (READ FIRST):\n"
            "  --compare-json FILE  <- JSON design, DESeq2 group_pairs shape:\n"
            "    {\"<tissue>\": {\"<contrast>\": {\n"
            "        \"control_group_name\": \"youth\",\n"
            "        \"experimental_group_name\": \"aged\",\n"
            "        \"control_samples\": [\"s1\", \"s2\"],\n"
            "        \"experimental_samples\": [\"s3\"]}}}\n"
            "  A multi-group contrast (several consecutive bars in ONE\n"
            "  stacked barplot) uses an ordered 'groups' list instead; the\n"
            "  FIRST group is the base level and every later group is\n"
            "  compared against it:\n"
            "    {\"<tissue>\": {\"<contrast>\": {\"groups\": [\n"
            "        {\"name\": \"youth\", \"samples\": [\"s1\", \"s2\"]}, ...]}}}\n"
            "  Multi-sample groups are model replicates and are summed\n"
            "  into ONE bar of the stacked barplot (x tick = group name).\n"
            "  DA status: ↑/↓ arrows sit directly after each cell-type\n"
            "  legend name, one glyph per comparison in that comparison's\n"
            "  color (e.g. a blue ↓ = youth vs aged), centered on the\n"
            "  name's row. The DA key legend below shows the ↑↓ pair as\n"
            "  the handle (leftmost, aligned with the cell-type swatches)\n"
            "  and the meaning as the label (ctrl vs exp; ↑ = higher in\n"
            "  exp, ↓ = lower in exp).\n"
            "  In every output:\n"
            "    - File names: da_sccoda_{contrast}_*\n"
            "  In the CSV covariate column, non-base groups are tagged\n"
            "  [T.<group>] so log2FC = log2(group / control).  control =\n"
            "  baseline / unperturbed; the later groups = experimental /\n"
            "  treated.  Convention: 'reference first, treatment second'."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  python sccoda.py \\
    --h5ad ovaries=/path/ovaries.h5ad \\
    --compare-json design.json \\
    --out-dir ./da_sccoda

  # design.json (pairwise + multi-group contrasts in one design):
  # {
  #   "ovaries": {
  #     "youth_vs_aged": {
  #       "control_group_name": "youth",
  #       "experimental_group_name": "aged",
  #       "control_samples": ["luanchao-21310-10XSC3"],
  #       "experimental_samples": ["luanchao-11238-10XSC3"]},
  #     "aging_series": {
  #       "groups": [
  #         {"name": "youth", "samples": ["s1", "s2"]},
  #         {"name": "aged", "samples": ["s3", "s4"]}]}}}
  #
  # explicit reference cell type and stricter credible threshold:
  python sccoda.py \\
    --h5ad uterus=/path/uterus.h5ad \\
    --compare-json design.json \\
    --out-dir ./da_sccoda --ref-celltype Fibroblast --credible-threshold 0.95
""")
    p.add_argument("-a", "--h5ad", action="append", required=True,
                   metavar="tissue=path",
                   help="Tissue name = h5ad path (repeatable)")
    p.add_argument("-c", "--compare-json", action="append", required=True,
                   metavar="FILE",
                   help="JSON design file (repeatable, merged): "
                        "{\"tissue\": {\"contrast\": entry}} where entry is "
                        "a DESeq2 group_pairs object "
                        "(control_group_name/experimental_group_name/"
                        "control_samples/experimental_samples) or an "
                        "ordered 'groups' list for a multi-group contrast. "
                        "The contrast key names the outputs and the "
                        "stacked-barplot x ticks use the group names")
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
    p.add_argument("--format", action="append", dest="formats", metavar="FMT",
                   help="Image format (png/pdf/svg/...). Repeatable. Default: png.")
    args = p.parse_args()
    args.formats = list(args.formats) if args.formats else ["png"]
    bad = [f for f in args.formats if f not in VALID_FORMATS]
    if bad:
        p.error(f"Unsupported image format(s): {', '.join(bad)}; "
                f"valid: {', '.join(VALID_FORMATS)}")
    return args


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

    comparisons = parse_compare_json(args.compare_json, tissues)
    run_params = {
        "formula": args.formula,
        "ref_celltype": args.ref_celltype,
        "credible_threshold": args.credible_threshold,
        "num_samples": args.num_samples,
        "num_warmup": args.num_warmup,
        "formats": ",".join(args.formats),
    }

    results: Dict[str, dict] = {}
    loaded: Dict[str, sc.AnnData] = {}
    for comp in comparisons:
        tissue = comp["tissue"]
        logger.info("=" * 60)
        logger.info("[%s] %s (%s)", tissue,
                    " vs ".join(n for n, _ in comp["groups"]), comp["label"])
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
            da_dir=out_dir / "da",
            image_formats=args.formats)
        results[comp["key"]] = {"summary": summary, "notes": notes}

    write_report(out_dir, comparisons, results, run_params)
    logger.info("DONE")


if __name__ == "__main__":
    main()
