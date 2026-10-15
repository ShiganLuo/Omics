#!/usr/bin/env python3
"""Per-cycle ecDNA structure plots for Decoil reconstruction results.

Draws the circular structure of every ecDNA cycle detected by Decoil from
``reconstruct.ecDNA.filtered.bed`` (fragment-level rows: chr, start, end,
circ_id, fragment_id, strand, coverage, estimated_proportions) joined with
``summary.txt`` (topology, size, proportion per circ_id).

Detection guard: when the filtered ecDNA BED has no data lines (no ecDNA
passed ``--filter-score``), no structure plot is drawn for cycles; a
placeholder figure stating "No ecDNA detected" is emitted instead so
downstream file expectations stay stable.

CLI conventions follow ``decoil_plot.py``: standalone bin/ script, logging
to stderr, image formats repeatable via ``-f``. All figure text is English.

The script runs inside decoil.sif (matplotlib + pandas available).
"""
from __future__ import annotations

import argparse
import logging
import os
import math
import re
import sys
from typing import Dict, List, Optional, Tuple

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Wedge
    import pandas as pd
except ImportError as e:  # pragma: no cover
    sys.stderr.write(
        f"Missing dependency: {e}\nInstall matplotlib and pandas in the "
        "execution environment (provided by the decoil conda env).\n"
    )
    raise

DPI = 300

# One color per fragment chromosome (cycled when a cycle spans >N chromosomes)
CHR_COLORS = [
    "#7A93B5", "#B58B7A", "#8FB57A", "#B57A9E", "#7AB5A9",
    "#A9A17A", "#7A7FB5", "#B57A7F", "#7AB5B0", "#93B57A",
]

MAX_CYCLES_PER_FIGURE = 12  # grid cap; extra cycles go to additional pages

# Default oncogene panel for highlight (replaceable via --onco-genes file,
# one gene name per line). Small curated subset of common ecDNA drivers.
DEFAULT_ONCO_GENES = {
    "MYC", "MYCN", "EGFR", "ERBB2", "MDM2", "MDM4", "CCND1", "CCND2",
    "CCNE1", "CDK4", "CDK6", "ALK", "NTRK1", "NTRK2", "NTRK3", "ROS1",
    "MET", "FGFR1", "FGFR2", "FGFR3", "PDGFRA", "KIT", "KRAS", "NRAS",
    "BRAF", "PIK3CA", "AKT1", "AKT2", "REL", "JAK2", "GLI1", "GLI2",
    "TP53", "RB1", "PTEN", "BRCA1", "BRCA2", "ATRX", "TERT", "YAP1",
    "MCL1", "BCL2", "BCL2L1", "IKZF1", "IKZF3", "TERC", "TERT",
}


def setup_logging() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("decoil_structure")


logger = setup_logging()


def load_ecdna_bed(path: str) -> pd.DataFrame:
    """Load fragment-level rows of reconstruct.ecDNA.filtered.bed.

    Parameters
    ----------
    path : str
        Path to the filtered ecDNA BED (tab-separated, header line starting
        with ``#chr``; comment-tolerant).

    Returns
    -------
    pandas.DataFrame
        Columns: chrom, start, end, circ_id, fragment_id, strand, coverage,
        proportion. Empty frame when no ecDNA was detected.

    Raises
    ------
    ValueError
        If the file exists but has fewer than 8 columns.
    """
    try:
        df = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
    except pd.errors.EmptyDataError:
        # Header-only or blank file: no ecDNA passed the filter
        return pd.DataFrame(
            columns=["chrom", "start", "end", "circ_id", "fragment_id",
                     "strand", "coverage", "proportion"])
    if df.empty:
        return pd.DataFrame(
            columns=["chrom", "start", "end", "circ_id", "fragment_id",
                     "strand", "coverage", "proportion"])
    if df.shape[1] < 8:
        raise ValueError(
            f"{path}: expected >= 8 columns, got {df.shape[1]}")
    df = df.iloc[:, :8]
    df.columns = ["chrom", "start", "end", "circ_id", "fragment_id",
                  "strand", "coverage", "proportion"]
    for col in ("start", "end"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    return df.dropna(subset=["start", "end"])


def load_summary(path: str) -> pd.DataFrame:
    """Load summary.txt keyed by circ_id.

    Parameters
    ----------
    path : str
        Path to ``summary.txt`` (tab-separated, has header).

    Returns
    -------
    pandas.DataFrame
        Indexed by circ_id (str); columns size_mb, topology, proportion.
    """
    df = pd.read_csv(path, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]
    out = pd.DataFrame({
        "circ_id": df["circ_id"].astype(str),
        "size_mb": pd.to_numeric(df["size(MB)"], errors="coerce"),
        "topology": df["topology_name"].astype(str),
        "proportion": pd.to_numeric(df.get("estimated_proportions"),
                                    errors="coerce"),
    })
    return out.set_index("circ_id")


def load_gtf_genes(path: str) -> pd.DataFrame:
    """Load gene intervals (chrom, start, end, gene_name) from a GTF.

    Parameters
    ----------
    path : str
        Path to a GTF (GENCODE/Ensembl style) with ``gene`` feature rows.

    Returns
    -------
    pandas.DataFrame
        Columns: chrom, start (0-based), end, gene_name. One row per gene.
    """
    rows = []
    gene_re = re.compile(r'gene_name "([^"]+)"')
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            m = gene_re.search(fields[8])
            rows.append((fields[0], int(fields[3]) - 1, int(fields[4]),
                         m.group(1) if m else fields[8][:20], fields[6]))
    return pd.DataFrame(rows, columns=["chrom", "start", "end",
                                       "gene_name", "strand"])


def load_gtf_exons(path: str) -> pd.DataFrame:
    """Load exon intervals with gene/strand info from a GTF.

    Parameters
    ----------
    path : str
        Path to a GTF (GENCODE/Ensembl style) with ``exon`` feature rows.

    Returns
    -------
    pandas.DataFrame
        Columns: chrom, start (0-based), end, gene_name, strand.
        One row per exon.
    """
    rows = []
    gene_re = re.compile(r'gene_name "([^"]+)"')
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            m = gene_re.search(fields[8])
            rows.append((fields[0], int(fields[3]) - 1, int(fields[4]),
                         m.group(1) if m else fields[8][:20], fields[6]))
    return pd.DataFrame(rows, columns=["chrom", "start", "end",
                                       "gene_name", "strand"])


def genes_on_fragment(genes: pd.DataFrame, chrom: str, start: int,
                      end: int) -> pd.DataFrame:
    """Return genes overlapping one fragment interval, sorted by position.

    Parameters
    ----------
    genes : pandas.DataFrame
        Gene intervals (see :func:`load_gtf_genes`); empty frame = no
        annotation available.
    chrom : str
        Fragment chromosome.
    start, end : int
        Fragment 0-based interval.

    Returns
    -------
    pandas.DataFrame
        Overlapping gene rows sorted by start; empty when none overlap.
    """
    if genes is None or genes.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "gene_name"])
    sub = genes[(genes["chrom"] == chrom)
                & (genes["end"] > start) & (genes["start"] < end)]
    return sub.sort_values("start")


def load_bw_mean_coverage(bw_path: "Optional[str]") -> float:
    """Return the genome-wide mean coverage of a bigwig (0.0 when missing).

    Parameters
    ----------
    bw_path : str or None
        Path to the coverage bigwig.

    Returns
    -------
    float
        Mean coverage over all chromosomes with valid stats; 0.0 when the
        bigwig is unavailable.
    """
    if not bw_path or not os.path.isfile(bw_path):
        return 0.0
    try:
        import pyBigWig
    except ImportError:
        return 0.0
    try:
        bw = pyBigWig.open(bw_path)
        vals = []
        for c, clen in bw.chroms().items():
            # skip tiny scaffolds to keep the reference stable
            if clen < 1_000_000:
                continue
            try:
                v = bw.stats(c, 0, clen, type="mean")[0]
                if v is not None:
                    vals.append(float(v))
            except RuntimeError:
                continue
        bw.close()
        return float(np.mean(vals)) if vals else 0.0
    except Exception as e:  # noqa: BLE001
        logger.warning(f"Failed to read bigwig mean {bw_path}: {e}")
        return 0.0


def load_onco_genes(path: "Optional[str]") -> set:
    """Load the oncogene highlight panel.

    Parameters
    ----------
    path : str or None
        Text file with one gene name per line; None/missing -> built-in
        DEFAULT_ONCO_GENES panel.

    Returns
    -------
    set of str
        Uppercased gene symbols to highlight.
    """
    if path and os.path.isfile(path):
        with open(path) as f:
            genes = {line.strip().upper() for line in f if line.strip()}
        logger.info(f"Loaded {len(genes)} oncogene panel entries from {path}")
        return genes
    return set(DEFAULT_ONCO_GENES)


def load_tpm(path: "Optional[str]") -> Dict[str, float]:
    """Load per-gene expression (TPM) values.

    Parameters
    ----------
    path : str or None
        Two-column TSV (gene_name, tpm) with header; None/missing -> {}.

    Returns
    -------
    dict
        gene_name -> TPM (float).
    """
    if not path or not os.path.isfile(path):
        return {}
    try:
        df = pd.read_csv(path, sep="\t")
        cols = {c.lower(): c for c in df.columns}
        gcol = cols.get("gene_name") or cols.get("gene") or df.columns[0]
        vcol = cols.get("tpm") or cols.get("tpm_expr") or df.columns[-1]
        return dict(zip(df[gcol].astype(str),
                        pd.to_numeric(df[vcol], errors="coerce")))
    except Exception as e:  # noqa: BLE001 - expression is optional input
        logger.warning(f"Failed to load TPM file {path}: {e}")
        return {}


def fragment_mean_coverage(bw_path: "Optional[str]", chrom: str,
                           start: int, end: int,
                           bins: int = 50) -> "Optional[np.ndarray]":
    """Read mean coverage profile across a fragment from a bigwig.

    Parameters
    ----------
    bw_path : str or None
        Path to the coverage bigwig; None/unreadable -> None.
    chrom : str
        Chromosome name.
    start, end : int
        Fragment interval (0-based).
    bins : int
        Number of evenly spaced bins along the fragment.

    Returns
    -------
    numpy.ndarray or None
        Per-bin mean coverage; None when bigwig unavailable.
    """
    if not bw_path or not os.path.isfile(bw_path):
        return None
    try:
        import pyBigWig
    except ImportError:
        logger.warning("pyBigWig not available; skipping dosage ring")
        return None
    try:
        bw = pyBigWig.open(bw_path)
        chroms = bw.chroms()
        if chrom not in chroms:
            return None
        end_c = min(end, chroms[chrom])
        if end_c <= start:
            return None
        edges = np.linspace(start, end_c, bins + 1).astype(int)
        vals = []
        for b0, b1 in zip(edges[:-1], edges[1:]):
            if b1 <= b0:
                vals.append(0.0)
                continue
            try:
                stats = bw.stats(chrom, int(b0), int(b1), type="mean")[0]
                vals.append(float(stats) if stats is not None else 0.0)
            except RuntimeError:  # region outside zoom levels
                vals.append(0.0)
        bw.close()
        return np.array(vals)
    except Exception as e:  # noqa: BLE001 - coverage is optional input
        logger.warning(f"Failed to read bigwig {bw_path}: {e}")
        return None



def draw_cycle(ax_cov: "plt.Axes", ax_gene: "plt.Axes",
               fragments: pd.DataFrame, info: "Optional[pd.Series]",
               circ_id: str, genes: pd.DataFrame, onco_genes: set,
               bw_path: "Optional[str]", tpm: "Optional[Dict[str, float]]",
               bw_mean_cov: float, exons_df: "Optional[pd.DataFrame]"
               ) -> "Optional[pd.DataFrame]":
    """Draw one ecDNA cycle as two stacked linear tracks (decoil-viz style).

    Track 1 (ax_cov):  coverage profile of the aligned reads (area fill),
                       with fragment spans drawn as gray backbone bars on
                       top and breakpoint junctions as red dashed vertical
                       lines spanning both tracks.
    Track 2 (ax_gene): protein-coding gene models (exon blocks + intron
                       lines + strand arrows), one row per gene, stacked
                       to avoid overlap.  Truncated genes (spanning a
                       breakpoint) get an orange dashed outline.

    Both axes share the same reference-genome x axis; the fragments are
    placed side by side with 20 kb padding, exactly like decoil-viz's
    ``windows = reduce(nodes + 200000)``.

    Parameters
    ----------
    ax_cov, ax_gene : matplotlib.axes.Axes
        Coverage and gene-model axes (share x after this call).
    fragments : pandas.DataFrame
        Fragment rows of one circ_id, kept in BED file order.
    info : pandas.Series or None
        Matching summary.txt row (size_mb, topology, proportion).
    circ_id : str
        Cycle identifier used in the title.
    genes : pandas.DataFrame
        Gene intervals (see :func:`load_gtf_genes`); empty = no annotation.
    onco_genes : set
        Uppercased gene symbols to highlight.
    bw_path : str, optional
        Coverage bigwig; None = skip the coverage track content.
    tpm : dict, optional
        gene_name -> TPM; appended to oncogene labels.
    bw_mean_cov : float, optional
        Genome-wide mean coverage (reserved for future normalization).
    exons_df : pandas.DataFrame, optional
        Exon intervals (see :func:`load_gtf_exons`).

    Returns
    -------
    pandas.DataFrame or None
        Per-gene table (gene, status, dosage, tpm) for this cycle;
        None when no genes overlap.
    """
    exons_df = exons_df if exons_df is not None else pd.DataFrame(
        columns=["chrom", "start", "end", "gene_name", "strand"])

    # ---- layout: side-by-side windows, PAD = 20 kb between fragments ----
    PAD = 20000
    spans = []           # (chrom, f_start, f_end, x_start, x_end)
    cursor = 0
    for _, frag in fragments.iterrows():
        f_start, f_end = int(frag["start"]), int(frag["end"])
        length = f_end - f_start
        spans.append((frag["chrom"], f_start, f_end,
                      cursor, cursor + length))
        cursor += length + PAD
    x_min, x_max = 0, cursor - PAD

    def to_x(chrom: str, g_start: int, g_end: int):
        """Map a genomic interval to canvas x, summed over fragments.

        Parameters
        ----------
        chrom : str
            Chromosome of the interval; must match a fragment.
        g_start, g_end : int
            Genomic coordinates.

        Returns
        -------
        tuple of float or None
            (x_start, x_end) inside the fragment's window, or None when
            the interval lies on an unknown chromosome.
        """
        for fc, fs, fe, xs, xe in spans:
            if fc != chrom:
                continue
            ov_s, ov_e = max(g_start, fs), min(g_end, fe)
            if ov_e <= ov_s:
                continue
            frac_s = (ov_s - fs) / float(fe - fs)
            frac_e = (ov_e - fs) / float(fe - fs)
            return (xs + (xe - xs) * frac_s, xs + (xe - xs) * frac_e)
        return None

    # ---- topology colors (muted, decoil-viz style) ----
    prop = info.get("proportion") if info is not None else None
    frag_bars = ["#4A4A4A", "#6E6E6E", "#8A8A8A", "#A5A5A5"]

    # ================= Track 1: reconstruction + coverage ===============
    for i, (fc, fs, fe, xs, xe) in enumerate(spans):
        col = frag_bars[i % len(frag_bars)]
        # backbone bar: per-fragment mean coverage shades the gray
        fcov = fragment_mean_coverage(bw_path, fc, fs, fe, bins=30)
        mean_cov = float(np.mean(fcov)) if fcov is not None else 1.0
        ref = max(float(bw_mean_cov), 1.0)
        alpha = 0.55 + 0.45 * min(np.log10(max(mean_cov, 1.0)) /
                                  max(np.log10(ref) * 4.0, 1e-6), 1.0)
        ax_cov.add_patch(plt.Rectangle(
            (xs, 0.25), xe - xs, 0.5, facecolor=col, alpha=alpha,
            edgecolor="white", linewidth=0.8, zorder=2))
        ax_cov.text(xs + (xe - xs) / 2.0, 0.5,
                    f"{fc}\n{fs:,}-{fe:,}", ha="center", va="center",
                    fontsize=5.0, color="white", zorder=3)

    # coverage profile per fragment
    if bw_path and os.path.isfile(bw_path):
        for fc, fs, fe, xs, xe in spans:
            cov = fragment_mean_coverage(bw_path, fc, fs, fe, bins=150)
            if cov is None or len(cov) < 2:
                continue
            cmax = float(cov.max())
            if cmax <= 0:
                continue
            xs_grid = np.linspace(xs, xe, len(cov))
            ax_cov.fill_between(xs_grid, 0.95, 0.95 + 1.4 * cov / cmax,
                                step=None, facecolor="#E5C07B",
                                alpha=0.55, linewidth=0, zorder=1)
            ax_cov.plot(xs_grid, 0.95 + 1.4 * cov / cmax,
                        color="#C9A44A", linewidth=0.6, zorder=2)
    ax_cov.set_ylim(0, 3.6)
    ax_cov.set_xlim(x_min - PAD / 2, x_max + PAD / 2)
    ax_cov.set_yticks([])
    ax_cov.set_ylabel("reconstruction / coverage", fontsize=6)
    ax_cov.tick_params(axis="x", labelsize=5, length=2)

    # breakpoint junction lines: red dashed vertical between fragment ends
    for i in range(len(spans) - 1):
        xbp = (spans[i][4] + spans[i + 1][3]) / 2.0
        for a in (ax_cov, ax_gene):
            a.axvline(xbp, color="#C0392B", linestyle=(0, (4, 3)),
                      linewidth=0.9, alpha=0.8, zorder=1)

    # circularity arcs: connect fragment ends that are joined in the
    # ecDNA circle — adjacency arcs between consecutive fragments ABOVE
    # the track, plus the closing arc from the last fragment's end back
    # to the first fragment's start (this is what makes the molecule a
    # CIRCLE in the linear rendering; decoil-viz/gTrack draws the same
    # "join" lines above the track).
    arc_top = 2.55
    for i in range(len(spans) - 1):
        x_l = spans[i][4]     # right edge of fragment i
        x_r = spans[i + 1][3]  # left edge of fragment i+1
        xm = (x_l + x_r) / 2.0
        half = max((x_r - x_l) * 0.6, 4000)
        t = np.linspace(0, np.pi, 40)
        ax_cov.plot(xm + half * np.cos(t), arc_top + 0.30 * np.sin(t),
                    color="#C0392B", linewidth=1.4, alpha=0.9, zorder=3)
        ax_cov.plot([x_l, x_l], [0.75, arc_top + 0.02], color="#C0392B",
                    linewidth=0.7, linestyle=(0, (2, 2)), zorder=3)
        ax_cov.plot([x_r, x_r], [0.75, arc_top + 0.02], color="#C0392B",
                    linewidth=0.7, linestyle=(0, (2, 2)), zorder=3)
        ax_cov.text(xm, arc_top + 0.40, f"B{i + 1}", ha="center",
                    fontsize=4.5, color="#C0392B")
    # closing arc: last fragment's right end -> first fragment's left end
    if len(spans) >= 2:
        half_total = spans[-1][4] - spans[0][3]
        t = np.linspace(0, np.pi, 60)
        ax_cov.plot(spans[0][3] + half_total * (0.5 + 0.5 * np.cos(t)),
                    arc_top + 0.72 * np.sin(t) + 0.18,
                    color="#8A4FBF", linewidth=1.5, alpha=0.9, zorder=3)
        ax_cov.text(spans[0][3] + half_total / 2.0, arc_top + 1.0,
                    "circular join", ha="center", fontsize=5.0,
                    color="#8A4FBF", style="italic")

    # ================= Track 2: gene models =============================
    exons_all = []
    gene_rows = []
    if genes is not None and not genes.empty:
        frag_chroms = {fc for fc, *_ in spans}
        gpool = genes[genes["chrom"].isin(frag_chroms)].copy()
        # keep only genes overlapping at least one fragment window
        keep = []
        for _, g in gpool.iterrows():
            hit = to_x_hit = to_x(g["chrom"], int(g["start"]), int(g["end"]))
            if hit is not None:
                keep.append(g)
        gpool = pd.DataFrame(keep) if keep else gpool.iloc[0:0]

    else:
        gpool = pd.DataFrame(columns=["chrom", "start", "end", "gene_name"])

    # stack genes into rows: sort by x start, greedy row assignment
    records = []
    for _, g in gpool.iterrows():
        hit = to_x(g["chrom"], int(g["start"]), int(g["end"]))
        if hit is None:
            continue
        x0, x1 = hit
        contained = (int(g["start"]) >= min(fs for fc, fs, fe, *_ in spans
                                            if fc == g["chrom"])) and                     (int(g["end"]) <= max(fe for fc, fs, fe, *_ in spans
                                          if fc == g["chrom"]))
        records.append({"g": g, "x0": x0, "x1": x1,
                        "broken": not contained})
    records.sort(key=lambda r: r["x0"])

    row_end = []   # x end of the last gene placed in each row
    for rec in records:
        placed_row = None
        for ri, xe in enumerate(row_end):
            if rec["x0"] > xe + PAD * 0.2:
                placed_row = ri
                break
        if placed_row is None:
            placed_row = len(row_end)
            row_end.append(rec["x1"])
        else:
            row_end[placed_row] = rec["x1"]
        rec["row"] = placed_row

    n_rows = max((r["row"] for r in records), default=0) + 1
    ROW_H = 1.0
    for rec in records:
        g = rec["g"]
        gname = str(g["gene_name"])
        strand = str(g.get("strand", "+"))
        onco = gname.upper() in onco_genes
        broken = rec["broken"]
        row_y = rec["row"] * ROW_H
        col = "#C97B4A" if broken else "#4C6E91"
        lwd = 1.4 if onco else 1.0

        # long introns (>30 kb) are drawn with a double-slash break:
        # two short segments + skip, instead of a dot-dash line crossing
        # the whole window (decoil-viz convention for gapped features)
        INTRON_BREAK = 30000
        if rec["x1"] - rec["x0"] > INTRON_BREAK:
            seg = (rec["x1"] - rec["x0"]) * 0.12
            ax_gene.plot([rec["x0"], rec["x0"] + seg],
                         [row_y, row_y], color=col, linewidth=0.7,
                         alpha=0.8, zorder=2)
            ax_gene.plot([rec["x1"] - seg, rec["x1"]],
                         [row_y, row_y], color=col, linewidth=0.7,
                         alpha=0.8, zorder=2)
            mid_x = (rec["x0"] + rec["x1"]) / 2.0
            for dx in (-60, 60):
                ax_gene.plot([mid_x + dx - 14, mid_x + dx + 14],
                             [row_y - 0.18, row_y + 0.18],
                             color=col, linewidth=0.8, alpha=0.9, zorder=2)
        else:
            ax_gene.plot([rec["x0"], rec["x1"]], [row_y, row_y],
                         color=col, linewidth=0.7, alpha=0.8, zorder=2)
        # exon blocks clipped to the fragments
        exons = exons_df[(exons_df["chrom"] == g["chrom"])
                         & (exons_df["gene_name"] == gname)]
        if exons.empty:
            ax_gene.plot([rec["x0"], rec["x1"]],
                         [row_y, row_y], color=col,
                         linewidth=3.5, solid_capstyle="butt", zorder=3)
        for _, ex in exons.sort_values("start").iterrows():
            hit = to_x(g["chrom"], int(ex["start"]), int(ex["end"]))
            if hit is None:
                continue
            x0, x1 = hit
            ax_gene.add_patch(plt.Rectangle(
                (x0, row_y - 0.28), x1 - x0, 0.56,
                facecolor=col, edgecolor="white", linewidth=0.3,
                zorder=3))
        # strand arrow at 3' end
        if strand == "+":
            tip_x = rec["x1"]
            ax_gene.annotate("", xy=(tip_x + 8, row_y),
                             xytext=(tip_x - 8, row_y),
                             arrowprops=dict(arrowstyle="->,head_width=0.25,"
                                             "head_length=0.5",
                                             color=col, lw=lwd), zorder=4)
        else:
            tip_x = rec["x0"]
            ax_gene.annotate("", xy=(tip_x - 8, row_y),
                             xytext=(tip_x + 8, row_y),
                             arrowprops=dict(arrowstyle="->,head_width=0.25,"
                                             "head_length=0.5",
                                             color=col, lw=lwd), zorder=4)
        # label: alternate above/below the row center to absorb
        # neighbor collisions; horizontal text at the gene midpoint
        label = gname + ("*" if broken else "")
        if onco and gname in (tpm or {}):
            label += f" | {tpm[gname]:.1f}"
        side = 1 if rec["row"] % 2 == 0 else -1
        ly = row_y + side * 0.42
        ax_gene.text((rec["x0"] + rec["x1"]) / 2.0, ly, label,
                     ha="center", va="center", fontsize=5.0,
                     color="#6B4E9E" if onco else "#2C3E50",
                     fontweight="bold" if onco else "normal", zorder=4)

        cov_arr = fragment_mean_coverage(bw_path, g["chrom"],
                                         int(g["start"]), int(g["end"]),
                                         bins=20)
        dosage = (float(np.median(cov_arr))
                  if cov_arr is not None and len(cov_arr) else np.nan)
        gene_rows.append({
            "circ_id": circ_id,
            "gene": gname,
            "status": "truncated" if broken else "intact",
            "onco": onco,
            "strand": strand,
            "region": f"{g['chrom']}:{int(g['start'])}-{int(g['end'])}",
            "median_coverage": round(dosage, 1) if pd.notna(dosage) else "",
            "tpm": tpm.get(gname, ""),
        })

    ax_gene.set_xlim(x_min - PAD / 2, x_max + PAD / 2)
    ax_gene.set_ylim(-1.0, n_rows * ROW_H + 1.0)
    ax_gene.set_yticks([])
    ax_gene.set_ylabel("genes (protein-coding)", fontsize=6)
    ax_gene.tick_params(axis="x", labelsize=5, length=2)
    ax_gene.set_xlabel("reference genome coordinate (windowed per fragment)",
                       fontsize=6)
    # hide inner x ticks on coverage axis; share x
    ax_cov.tick_params(labelbottom=False)
    for spine in ("top", "right", "left"):
        ax_cov.spines[spine].set_visible(False)
        ax_gene.spines[spine].set_visible(False)

    title_topo = str(info["topology"]) if info is not None else "NA"
    size_kb = ((info["size_mb"] * 1000)
               if info is not None and pd.notna(info.get("size_mb")) else None)
    ax_cov.set_title(f"circ {circ_id}  ·  {title_topo}  ·  "
                     f"{size_kb:.1f} kb" if size_kb else
                     f"circ {circ_id}  ·  {title_topo}",
                     fontsize=8, loc="left")
    return pd.DataFrame(gene_rows) if gene_rows else None


def plot_placeholder(out_prefix: str, image_formats: List[str],
                     sample_id: str) -> None:
    """Emit a 'No ecDNA detected' placeholder figure.

    Parameters
    ----------
    out_prefix : str
        Output path prefix; ``.<fmt>`` appended per format.
    image_formats : list of str
        Image formats (png/pdf/...).
    sample_id : str
        Sample label used in the title.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.axis("off")
    ax.text(0.5, 0.55, "No ecDNA detected", ha="center", va="center",
            fontsize=16, color="#7F8C8D")
    ax.text(0.5, 0.40, f"sample: {sample_id}",
            ha="center", va="center", fontsize=10, color="#95A5A6")
    ax.text(0.5, 0.30,
            "reconstruct.ecDNA.filtered.bed contains no cycles",
            ha="center", va="center", fontsize=8, color="#95A5A6")
    for fmt in image_formats:
        path = f"{out_prefix}.{fmt}"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        logger.info(f"Saved {path} (placeholder: no ecDNA detected)")
    plt.close(fig)


def plot_structures(bed: pd.DataFrame, summary: pd.DataFrame,
                    genes: pd.DataFrame, onco_genes: set,
                    bw_path: "Optional[str]", tpm: Dict[str, float],
                    sample_id: str, out_prefix: str,
                    image_formats: List[str], max_per_fig: int,
                    one_file_per_cycle: bool = False,
                    exons_df: "Optional[pd.DataFrame]" = None) -> None:
    """Draw every detected cycle's structure.

    By default cycles are arranged in a grid (max_per_fig per page). With
    ``one_file_per_cycle`` each cycle is saved as its own image file
    (``<out_prefix>_circ<id>.<fmt>``), one ecDNA per figure, for
    slide/report use.

    Also writes a per-gene functional annotation table
    ``<out_prefix}_genes.tsv`` covering all cycles.

    Parameters
    ----------
    bed : pandas.DataFrame
        Fragment-level rows (see :func:`load_ecdna_bed`).
    summary : pandas.DataFrame
        Summary rows indexed by circ_id (see :func:`load_summary`).
    genes : pandas.DataFrame
        Gene intervals (see :func:`load_gtf_genes`); empty = no annotation.
    onco_genes : set
        Uppercased gene symbols to highlight.
    bw_path : str or None
        Coverage bigwig for the dosage ring; None = skip layer.
    tpm : dict
        gene_name -> TPM (empty = no expression data).
    sample_id : str
        Sample label used in figure suptitles.
    out_prefix : str
        Output path prefix; ``.<fmt>`` / ``.<page>.<fmt>`` appended (grid
        mode) or ``_circ<id>.<fmt>`` (per-cycle mode).
    image_formats : list of str
        Image formats (png/pdf/...).
    max_per_fig : int
        Maximum cycles drawn on one figure page (grid mode only).
    one_file_per_cycle : bool, optional
        Save one image per ecDNA cycle instead of a combined grid.
    """
    circ_ids = list(dict.fromkeys(bed["circ_id"].astype(str)))
    all_gene_rows = []
    bw_mean_cov = load_bw_mean_coverage(bw_path)
    if bw_path:
        logger.info(f"Genome-wide mean coverage: {bw_mean_cov:.2f}x "
                    "(dosage ring reference)")

    def _draw_and_collect(fig: "plt.Figure", cid: str) -> None:
        """Draw one cycle as two stacked track axes on the figure.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Target figure; two stacked axes are created internally.
        cid : str
            Cycle identifier.
        """
        frags = bed[bed["circ_id"].astype(str) == cid]
        info = summary.loc[cid] if cid in summary.index else None
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.4], hspace=0.12)
        ax_cov = fig.add_subplot(gs[0])
        ax_gene = fig.add_subplot(gs[1], sharex=ax_cov)
        rows = draw_cycle(ax_cov, ax_gene, frags, info, cid, genes,
                          onco_genes, bw_path, tpm, bw_mean_cov, exons_df)
        if rows is not None:
            all_gene_rows.append(rows)

    def _grid_page(page_ids: List[str], page: int, n_pages: int,
                   suffix: str) -> None:
        """Draw one combined grid page and save it.

        Parameters
        ----------
        page_ids : list of str
            Cycle ids on this page.
        page : int
            Zero-based page index.
        n_pages : int
            Total number of pages.
        suffix : str
            Filename suffix (``.p2`` when multi-page, else empty).
        """
        for k, cid in enumerate(page_ids):
            fig = plt.figure(figsize=(9.0, 5.2))
            _draw_and_collect(fig, cid)
            fig.suptitle(f"{sample_id}: ecDNA structures "
                         f"({len(circ_ids)} cycles)", fontsize=10, y=1.02)
            fig.tight_layout()
        for fmt in image_formats:
            path = f"{out_prefix}{suffix}.{fmt}"
            fig.savefig(path, dpi=DPI, bbox_inches="tight")
            logger.info(f"Saved {path}")
        plt.close(fig)

    if one_file_per_cycle:
        # One ecDNA, one figure: single large square canvas per cycle.
        # The combined grid page is still emitted afterwards so downstream
        # rules keep a single summary image.
        os.makedirs(os.path.dirname(os.path.abspath(out_prefix)),
                    exist_ok=True)
        for cid in circ_ids:
            fig = plt.figure(figsize=(9.0, 5.2))
            _draw_and_collect(fig, cid)
            for fmt in image_formats:
                path = f"{out_prefix}_circ{cid}.{fmt}"
                fig.savefig(path, dpi=DPI, bbox_inches="tight")
                logger.info(f"Saved {path}")
            plt.close(fig)
        n_pages = (len(circ_ids) + max_per_fig - 1) // max_per_fig
        if n_pages == 1:
            _grid_page(circ_ids, 0, 1, "")
        else:
            for page in range(n_pages):
                _grid_page(circ_ids[page * max_per_fig:
                                    (page + 1) * max_per_fig],
                           page, n_pages, f".p{page + 1}")
    else:
        n_pages = (len(circ_ids) + max_per_fig - 1) // max_per_fig
        if n_pages == 1:
            _grid_page(circ_ids, 0, 1, "")
        else:
            for page in range(n_pages):
                _grid_page(circ_ids[page * max_per_fig:
                                    (page + 1) * max_per_fig],
                           page, n_pages, f".p{page + 1}")

    # Per-gene functional table (all cycles)
    if all_gene_rows:
        table = pd.concat(all_gene_rows, ignore_index=True)
        table_path = f"{out_prefix}_genes.tsv"
        table.to_csv(table_path, sep="\t", index=False)
        logger.info(f"Saved {table_path} ({len(table)} gene annotations)")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Per-cycle ecDNA structure plots for Decoil results.")
    parser.add_argument("--sample", required=True, metavar="SAMPLE_ID",
                        help="sample id used in titles")
    parser.add_argument("--ecdna-bed", required=True,
                        metavar="BED",
                        help="reconstruct.ecDNA.filtered.bed")
    parser.add_argument("--summary", required=True, metavar="TSV",
                        help="decoil summary.txt")
    parser.add_argument("--gtf", default=None, metavar="GTF",
                        help="gene annotation GTF; when omitted or empty "
                             "cycles are drawn without gene annotation")
    parser.add_argument("--bw", default=None, metavar="BIGWIG",
                        help="coverage bigwig; enables the inner dosage "
                             "ring + per-gene median coverage")
    parser.add_argument("--onco-genes", default=None, metavar="FILE",
                        help="text file with one gene name per line to "
                             "highlight; default: built-in panel")
    parser.add_argument("--tpm", default=None, metavar="TSV",
                        help="two-column TSV (gene_name, tpm); adds |TPM to "
                             "oncogene labels and fills the tpm column "
                             "(reserved for future expression coupling)")
    parser.add_argument("-o", "--out-prefix", required=True,
                        help="output path prefix (extensions appended)")
    parser.add_argument("-f", "--format", action="append", dest="formats",
                        metavar="FMT", help="image format, repeatable; "
                                            "default: png")
    parser.add_argument("--max-per-fig", type=int,
                        default=MAX_CYCLES_PER_FIGURE,
                        help="max cycles per figure page")
    parser.add_argument("--one-file-per-cycle", action="store_true",
                        help="save one image per ecDNA cycle "
                             "(<out-prefix>_circ<id>.<fmt>) instead of a "
                             "combined grid page")
    args = parser.parse_args()

    image_formats = args.formats or ["png"]

    bed = load_ecdna_bed(args.ecdna_bed)
    summary = load_summary(args.summary)

    genes = pd.DataFrame(columns=["chrom", "start", "end", "gene_name"])
    exons = None
    if args.gtf:
        if os.path.isfile(args.gtf):
            genes = load_gtf_genes(args.gtf)
            logger.info(f"Loaded {len(genes)} genes from {args.gtf}")
            exons = load_gtf_exons(args.gtf)
            logger.info(f"Loaded {len(exons)} exons from {args.gtf}")
        else:
            logger.warning(f"GTF not found ({args.gtf}); "
                           "drawing without gene annotation")

    if bed.empty:
        # Detection guard: nothing passed the filter -> no structures to draw
        logger.info(f"{args.sample}: no ecDNA detected, drawing placeholder")
        plot_placeholder(args.out_prefix, image_formats, args.sample)
        return

    logger.info(f"{args.sample}: {bed['circ_id'].nunique()} ecDNA cycles, "
                f"{len(bed)} fragments")
    onco_genes = load_onco_genes(args.onco_genes)
    tpm = load_tpm(args.tpm)
    plot_structures(bed, summary, genes, onco_genes, args.bw, tpm,
                    args.sample, args.out_prefix,
                    image_formats, args.max_per_fig,
                    one_file_per_cycle=args.one_file_per_cycle,
                    exons_df=exons)


if __name__ == "__main__":
    main()
