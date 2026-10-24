#!/usr/bin/env python3
"""Evidence-3 figure: per-fragment RNA/background ratio heatmap, grouped by family.

Renders a single-panel heatmap of circle RNA expression evidence:

    x axis  : samples (columns)
    y axis  : circle fragments, grouped and colored by breakpoint family
    cell    : log2(mean-of-reps of ratio_vs_flank), annotated with the linear
              ratio (e.g. ``6.2x``); ``n.e.`` marks fragments whose linear
              background window contained zero RNA reads (ratio undefined)
    left    : one strength bar per family, colored by the family-mean log2
              ratio, with the family tag (f1..f9) inside
    right   : colorbar of the log2 ratio scale (clipped to [-3, 3])

Input
-----
Per-sample-replicate TSVs written by ``rna_dna_ratio.py`` (masked-flank
version). Required columns::

    circ_id          circle id (decoil), e.g. 104
    fragment         genomic interval ``chr:start-end``
    ratio_vs_flank   fragment RNA rpk / circle-masked linear-flank rpk
                     (``inf`` = empty background, rendered as n.e.)

File naming: ``<sample_prefix>_<rep>.tsv``. The sample prefix must be mapped
to a display name via ``--sample-map`` (comma-separated ``PREFIX=Name``
items), e.g. ``DMS_PO6=DMSO_P06,DMSO_P20=DMSO_P20``. Only samples present in
this mapping are plotted; column order follows the mapping order.

Family assignment
-----------------
``(sample, circ_id) -> family`` comes from ``--family-tsv`` (a TSV with
columns ``sample``, ``circ_id``, ``family``; e.g. the output of
``breakpoint_families.py``). Rows unassigned to a family are dropped.

Examples
--------
    python plot_rna_dna_ratio.py \
        -d rna_dna_ratio/ \
        --sample-map DMS_PO6=DMSO_P06,DMSO_P20=DMSO_P20,PlaB_PO6=PlaB_P06,PlaB_P20=PlaB_P20 \
        --family-tsv breakpoint_families.tsv \
        -o rna_dna_ratio_evidence -f png -f pdf
"""
from __future__ import annotations

import argparse
import logging
import os
from collections import defaultdict
from typing import Dict, List, Literal, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("plot_rna_dna_ratio")

#: Export file formats supported by :func:`main`.
PlotFormat = Literal["png", "pdf", "svg"]

#: Fill color for cells whose ratio is undefined (empty background window).
NE_COLOR = "#e8e8e8"
#: Ratio color scale (log2 units); values outside are clipped to the ends.
VMIN, VMAX = -3.0, 3.0
#: Fixed color per family index (family 1..9 -> position 0..8).
FAMILY_COLORS = ["#4a708b", "#b2182b", "#2f8f4e", "#7a5195", "#c9a227",
                 "#d95f02", "#1b9e77", "#7570b3", "#8c8c8c"]


def family_color(fam: Optional[int]) -> str:
    """Return the palette color for a family number.

    Parameters
    ----------
    fam : int or None
        Family number (1-based). ``None`` (unassigned) maps to the last
        (gray) palette entry.

    Returns
    -------
    str
        Hex color string.
    """
    if fam is None:
        return FAMILY_COLORS[-1]
    return FAMILY_COLORS[(fam - 1) % len(FAMILY_COLORS)]


def load_inputs(datadir: str, sample_map: Dict[str, str]) -> Dict[Tuple[str, int], pd.DataFrame]:
    """Load per-sample-replicate ratio TSVs into a keyed dict.

    A file qualifies when its name is ``<prefix>_<rep>.tsv`` and ``prefix``
    is a key of ``sample_map``. ``inf`` ratios (empty background windows)
    are converted to ``NaN`` so downstream averaging treats them as missing.

    Parameters
    ----------
    datadir : str
        Directory containing the ``<prefix>_<rep>.tsv`` files.
    sample_map : dict of str -> str
        File-name prefix -> display sample name. Column order of the figure
        follows the insertion order of this dict.

    Returns
    -------
    dict
        ``(display_name, rep) -> pandas.DataFrame`` with the TSV columns plus
        a normalized numeric ``ratio_vs_flank`` column (inf -> NaN).

    Raises
    ------
    FileNotFoundError
        If ``datadir`` does not exist.
    ValueError
        If no qualifying TSVs were found.
    """
    if not os.path.isdir(datadir):
        raise FileNotFoundError(f"data directory not found: {datadir}")
    out: Dict[Tuple[str, int], pd.DataFrame] = {}
    for fn in sorted(os.listdir(datadir)):
        if not fn.endswith(".tsv"):
            continue
        stem = fn[:-4]
        prefix, _, rep = stem.rpartition("_")
        if prefix not in sample_map or not rep.isdigit():
            continue
        df = pd.read_csv(os.path.join(datadir, fn), sep="\t")
        df["ratio_vs_flank"] = df["ratio_vs_flank"].replace(float("inf"), np.nan)
        out[(sample_map[prefix], int(rep))] = df
    if not out:
        raise ValueError(
            f"no input TSVs matched in {datadir}; expected files named "
            f"<prefix>_<rep>.tsv with prefixes {sorted(sample_map)}")
    return out


def load_family_map(path: str, sample_map: Dict[str, str]) -> Dict[Tuple[str, str], int]:
    """Read the (sample, circ_id) -> family assignment table.

    Parameters
    ----------
    path : str
        TSV with header columns ``sample``, ``circ_id``, ``family``.
        ``sample`` must use display names (values of ``sample_map``).
    sample_map : dict of str -> str
        File prefix -> display name mapping (for the error message only).

    Returns
    -------
    dict
        ``(display_name, circ_id_as_str) -> family number``.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If required columns are missing.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"family table not found: {path}")
    df = pd.read_csv(path, sep="\t")
    missing = {"sample", "circ_id", "family"} - set(df.columns)
    if missing:
        raise ValueError(
            f"{path}: missing required column(s) {sorted(missing)}; need "
            "'sample', 'circ_id', 'family'. sample values must be display "
            f"names (from --sample-map values, e.g. {sorted(set(sample_map.values()))})")
    return {(str(r["sample"]), str(r["circ_id"])): int(r["family"])
            for _, r in df.iterrows()}


def collect_cells(data: Dict[Tuple[str, int], pd.DataFrame]) -> dict:
    """Collapse replicates to per-(sample, circle, fragment) ratio values.

    Parameters
    ----------
    data : dict
        ``(display_name, rep) -> DataFrame`` as returned by :func:`load_inputs`.

    Returns
    -------
    dict
        ``(display_name, circ_id, fragment) -> {rep: ratio}`` where ratio is
        a float (possibly NaN for empty-background replicates).
    """
    cells: dict = {}
    for (sample, rep), df in data.items():
        for _, r in df.iterrows():
            key = (sample, str(r["circ_id"]), r["fragment"])
            cells.setdefault(key, {})[rep] = float(r["ratio_vs_flank"])
    return cells


def build_rows(cells: dict, family_map: Dict[Tuple[str, str], int],
               samples: List[str]) -> List[Tuple[Optional[int], tuple]]:
    """Order fragments into heatmap rows: family asc, then genomic coordinate.

    Keys without a family assignment are skipped (with a warning) so every
    rendered row carries family information.

    Parameters
    ----------
    cells : dict
        Output of :func:`collect_cells`.
    family_map : dict
        Output of :func:`load_family_map`.
    samples : list of str
        Display names in column order; rows whose sample is not listed are
        skipped.

    Returns
    -------
    list
        ``(family, (display_name, circ_id, fragment, {rep: ratio}))`` per
        heatmap row, top to bottom.
    """
    fam_units: Dict[Optional[int], list] = defaultdict(list)
    dropped = 0
    for (sample, circ, frag), reps in cells.items():
        if sample not in samples:
            continue
        fam = family_map.get((sample, circ))
        if fam is None:
            dropped += 1
            continue
        fam_units[fam].append((sample, circ, frag, reps))
    if dropped:
        logger.warning("%d fragment(s) without family assignment dropped", dropped)
    ordered: List[Tuple[Optional[int], tuple]] = []
    for fam in sorted(fam_units, key=lambda f: (f is None, f if f is not None else 0)):
        units = sorted(fam_units[fam],
                       key=lambda u: (u[2].split(":")[0],
                                      int(u[2].split(":")[1].split("-")[0])))
        for u in units:
            ordered.append((fam, u))
    return ordered


def ratio_matrix(ordered: List[Tuple[Optional[int], tuple]],
                 samples: List[str]) -> np.ndarray:
    """Build the log2 ratio matrix (rows x samples) from ordered rows.

    Each cell holds ``log2`` of the replicate-mean ratio for that
    (fragment-row, sample) pair; samples other than the circle's own are NaN
    (a circle is observed in one sample only). Values are clipped to
    ``[VMIN, VMAX]`` so extremes render as saturated colors instead of
    wrapping around the diverging colormap.

    Parameters
    ----------
    ordered : list
        Output of :func:`build_rows`.
    samples : list of str
        Display names in column order.

    Returns
    -------
    numpy.ndarray
        ``(n_rows, n_samples)`` float array, NaN = no estimate.
    """
    n_rows, n_cols = len(ordered), len(samples)
    mat = np.full((n_rows, n_cols), np.nan)
    for ri, (_, (sample, _, _, reps)) in enumerate(ordered):
        ci = samples.index(sample)
        vals = [v for v in reps.values() if np.isfinite(v)]
        if vals:
            mat[ri, ci] = np.log2(sum(vals) / len(vals))
    # NOTE: values are intentionally NOT clipped here -- annotations must
    # show the true ratio; callers clip only for color mapping.
    return mat


def family_runs(ordered: List[Tuple[Optional[int], tuple]]) -> List[Tuple[Optional[int], int, int]]:
    """Group consecutive rows of the same family into (family, first, last).

    Parameters
    ----------
    ordered : list
        Output of :func:`build_rows` (rows are already family-grouped).

    Returns
    -------
    list
        ``(family, row_start, row_end)`` with inclusive 0-based row indices.
    """
    row_fams = [fam for fam, _ in ordered]
    runs: List[Tuple[Optional[int], int, int]] = []
    start = 0
    for i in range(1, len(row_fams) + 1):
        if i == len(row_fams) or row_fams[i] != row_fams[start]:
            runs.append((row_fams[start], start, i - 1))
            start = i
    return runs


def draw_strength_bars(ax_bar: plt.Axes, mat: np.ndarray,
                       runs: List[Tuple[Optional[int], int, int]],
                       ordered: List[Tuple[Optional[int], tuple]]) -> None:
    """Draw one strength bar per family into the dedicated left axes.

    Bar color encodes the family-mean log2 ratio through the same diverging
    scale as the heatmap (NaN-only families render in the n.e. gray). The
    family tag (f1, f2, ...) is written inside the bar; text flips to white
    on dark fills. The axes is left with no ticks or spines.

    Parameters
    ----------
    ax_bar : matplotlib.axes.Axes
        Narrow axes left of the heatmap; its y limits are set to match the
        heatmap rows exactly.
    mat : numpy.ndarray
        Output of :func:`ratio_matrix`.
    runs : list
        Output of :func:`family_runs`.
    ordered : list
        Output of :func:`build_rows` (for row counts).
    """
    n_rows = mat.shape[0]
    cmap = plt.get_cmap("RdBu_r").copy()
    norm = plt.Normalize(vmin=VMIN, vmax=VMAX)
    for fam, r0, r1 in runs:
        vals = [v for v in mat[r0:r1 + 1, :].ravel() if np.isfinite(v)]
        if vals:
            mean_log2 = sum(vals) / len(vals)
            fill = cmap(norm(np.clip(mean_log2, VMIN, VMAX)))
        else:
            mean_log2, fill = None, NE_COLOR
        y0, y1 = r0 - 0.5, r1 + 0.5
        ax_bar.add_patch(Rectangle((0.0, y0), 1.0, y1 - y0,
                                   facecolor=fill, edgecolor="#333333",
                                   lw=0.4, zorder=3))
        ax_bar.text(0.5, (y0 + y1) / 2, f"f{fam}" if fam is not None else "-",
                    ha="center", va="center", fontsize=5.2,
                    color="white" if mean_log2 is not None and abs(mean_log2) > 1.4
                    else "#222222")
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(n_rows - 0.5, -0.5)
    ax_bar.set_xticks([])
    # NOTE: y ticks are NOT touched here -- plot_heatmap sets the fragment
    # labels on this axes and must win regardless of call order.
    for sp in ax_bar.spines.values():
        sp.set_visible(False)


def plot_heatmap(ax: plt.Axes, ax_bar: plt.Axes, mat: np.ndarray,
                 ordered: List[Tuple[Optional[int], tuple]],
                 runs: List[Tuple[Optional[int], int, int]],
                 samples: List[str]) -> plt.Axes:
    """Render the ratio heatmap plus annotations, labels and colorbar.

    Cell text shows the linear ratio (``2**log2``, one decimal, ``x``);
    NaN cells show ``n.e.`` when every replicate had an empty background
    window, otherwise ``-`` (replicate drop-out). y tick labels (fragment
    coordinates + circle id) live on ``ax_bar`` and are colored by family;
    the colorbar sits on the right at heatmap height.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Heatmap axes.
    ax_bar : matplotlib.axes.Axes
        Strength-bar axes (see :func:`draw_strength_bars`); receives the
        y tick labels.
    mat : numpy.ndarray
        Output of :func:`ratio_matrix`.
    ordered : list
        Output of :func:`build_rows`.
    runs : list
        Output of :func:`family_runs`.
    samples : list of str
        Display names in column order.

    Returns
    -------
    matplotlib.axes.Axes
        The colorbar axes (for layout adjustments by the caller).
    """
    n_rows, n_cols = mat.shape
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(NE_COLOR)
    norm = plt.Normalize(vmin=VMIN, vmax=VMAX)
    im = ax.imshow(np.clip(mat, VMIN, VMAX), aspect="auto", cmap=cmap,
                   norm=norm)

    # cell annotations: linear ratio text; white text on dark fills
    for ri, (_, (sample, circ, frag, reps)) in enumerate(ordered):
        for ci in range(n_cols):
            v = mat[ri, ci]
            if np.isnan(v):
                if samples[ci] == sample:
                    all_empty = all(not np.isfinite(x) for x in reps.values())
                    ax.text(ci, ri, "n.e." if all_empty else "-",
                            ha="center", va="center", fontsize=5.0,
                            color="#777777")
            else:
                ax.text(ci, ri, f"{2.0 ** v:.1f}x", ha="center", va="center",
                        fontsize=5.0,
                        color="white" if abs(v) > 1.4 else "#222222")

    # colorbar: right of the heatmap, same height
    cbar = plt.colorbar(im, ax=ax, location="right", fraction=0.03, pad=0.03,
                        ticks=[-3, -2, -1, 0, 1, 2, 3])
    cbar.ax.set_ylabel("log2 (fragment RNA / masked flank)", fontsize=5.5)
    cbar.ax.tick_params(labelsize=5.5)

    # y tick labels on the bar axes, colored by family; hide heatmap's own
    ax.set_yticks([])
    ax_bar.set_yticks(range(n_rows))
    ax_bar.set_yticklabels([f"{frag.split(':')[1]} ({circ})"
                            for _, (_, circ, frag, _) in ordered])
    for tick, (fam, _) in zip(ax_bar.get_yticklabels(), ordered):
        tick.set_color(family_color(fam))
    ax_bar.tick_params(axis="y", labelsize=4.8, length=0, pad=2.5)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(samples, fontsize=6.5, rotation=30, ha="right")
    ax.set_title("Circle fragment RNA / masked linear flank, by family",
                 fontsize=10, loc="left", pad=6)
    return cbar.ax


def run(datadir: str, out_prefix: str,
        sample_map: Dict[str, str], family_tsv: str,
        formats: Optional[List[PlotFormat]] = None) -> None:
    """Render the evidence-3 figure end to end.

    Parameters
    ----------
    datadir : str
        Directory with per-replicate ratio TSVs (``<prefix>_<rep>.tsv``).
    out_prefix : str
        Output path prefix; each requested format appends its extension.
    sample_map : dict of str -> str
        File-name prefix -> display sample name, in column order.
    family_tsv : str
        Path to the family assignment TSV (``sample``, ``circ_id``,
        ``family`` columns; ``sample`` uses display names).
    formats : list of str, optional
        Export formats (``png``/``pdf``/``svg``). Default ``["png"]``.

    Raises
    ------
    ValueError
        If no input TSVs match or the family table lacks columns.
    FileNotFoundError
        If ``datadir`` or ``family_tsv`` does not exist.
    """
    if formats is None:
        formats = ["png"]
    samples = list(sample_map.values())
    data = load_inputs(datadir, sample_map)
    family_map = load_family_map(family_tsv, sample_map)

    cells = collect_cells(data)
    ordered = build_rows(cells, family_map, samples)
    mat = ratio_matrix(ordered, samples)
    runs = family_runs(ordered)
    logger.info("plotted %d fragment rows, %d families, samples=%s",
                len(ordered), len(runs), samples)

    fig = plt.figure(figsize=(7.2, 5.8))
    # explicit placement: strength bars hug the heatmap left edge (no gap);
    # y labels sit left of the bars. No category legend: every cell carries
    # its exact ratio text, the colorbar defines the scale, and the old
    # >=2x/1-2x/<1x strip duplicated that information.
    ax_bar = fig.add_axes([0.235, 0.12, 0.014, 0.72])
    ax_b = fig.add_axes([0.249, 0.12, 0.50, 0.72])
    plot_heatmap(ax_b, ax_bar, mat, ordered, runs, samples)
    draw_strength_bars(ax_bar, mat, runs, ordered)
    for fmt in formats:
        path = f"{out_prefix}.{fmt}"
        fig.savefig(path, dpi=300)
        logger.info("saved %s", path)


def parse_sample_map(spec: str) -> Dict[str, str]:
    """Parse the ``--sample-map`` argument into an ordered dict.

    Parameters
    ----------
    spec : str
        Comma-separated ``PREFIX=Name`` items, e.g.
        ``"DMS_PO6=DMSO_P06,DMSO_P20=DMSO_P20"``. Order defines column order.

    Returns
    -------
    dict
        ``{prefix: display_name}`` in the given order.

    Raises
    ------
    ValueError
        If any item lacks ``=`` or a prefix/name is empty.
    """
    out: Dict[str, str] = {}
    for item in spec.split(","):
        prefix, sep, name = item.partition("=")
        if not sep or not prefix.strip() or not name.strip():
            raise ValueError(
                f"--sample-map item '{item}' invalid; expected PREFIX=Name")
        out[prefix.strip()] = name.strip()
    return out


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Build the CLI parser.

    Parameters
    ----------
    argv : list of str, optional
        Argument tokens; ``None`` uses ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        Parsed arguments: ``datadir``, ``out_prefix``, ``formats``,
        ``sample_map`` (dict), ``family_tsv``.
    """
    parser = argparse.ArgumentParser(
        description=(__doc__ or "").splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--datadir", required=True,
                        help="Directory with per-replicate ratio TSVs named "
                             "<prefix>_<rep>.tsv (written by rna_dna_ratio.py).")
    parser.add_argument("-o", "--out-prefix", default="rna_dna_ratio_evidence",
                        help="Output path prefix without extension "
                             "(default: %(default)s).")
    parser.add_argument("--sample-map", required=True,
                        help="Comma-separated PREFIX=Name items mapping TSV "
                             "name prefixes to display sample names, in column "
                             "order. Example: DMS_PO6=DMSO_P06,DMSO_P20=DMSO_P20")
    parser.add_argument("--family-tsv", required=True,
                        help="TSV with 'sample', 'circ_id', 'family' columns "
                             "(e.g. breakpoint_families output); 'sample' uses "
                             "display names.")
    parser.add_argument("-f", "--format", action="append", dest="formats",
                        metavar="FMT", choices=["png", "pdf", "svg"],
                        help="Image format, repeatable (default: png).")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point.

    Parameters
    ----------
    argv : list of str, optional
        Argument tokens; ``None`` uses ``sys.argv[1:]``.

    Returns
    -------
    int
        Process exit code (0 on success).
    """
    args = parse_args(argv)
    run(datadir=args.datadir, out_prefix=args.out_prefix,
        sample_map=parse_sample_map(args.sample_map),
        family_tsv=args.family_tsv, formats=args.formats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
