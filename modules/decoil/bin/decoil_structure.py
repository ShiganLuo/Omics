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
import matplotlib.patches as patches
try:
    import gffutils
except ImportError:  # gffutils is only needed by create_db(); optional here
    gffutils = None
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import pandas as pd
import pyBigWig
# ---------------------------------------------------------------------------
# CGView-style circular renderer constants
# ---------------------------------------------------------------------------
R_RULER = 1.0           # outer ruler radius
R_COV_OUT = 0.90        # coverage ring outer
R_COV_IN = 0.72         # coverage ring inner (radial height 0.18)
R_GENE = 0.545          # outermost gene lane radius
LANE_W = 0.085          # gene lane width (exon blocks legible)
EXON_W = 0.068          # exon block radial width (0.8 * LANE_W)

# reference-style palette
CORAL = "#E8734F"       # fragment boundary markers
GOLD = "#E3B34D"        # coverage fill
GOLD_D = "#B98A2F"      # coverage outline
MAGENTA = "#A0457A"     # intact gene blocks
INDIGO = "#3A3A70"      # ticks / text
GRAY = "#8A94A0"        # secondary text
TRUNCATED = "#C97B4A"   # truncated gene colour
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

def create_db(gtf_file: str) -> gffutils.FeatureDB:
    """Create or load a gffutils database from a GTF/GFF3 annotation file.

    If the database file (``<gtf_file>.db``) already exists on disk it is
    reused; otherwise a new one is created.

    Parameters
    ----------
    gtf_file : str
        Path to a GTF or GFF3 gene annotation file.

    Returns
    -------
    gffutils.FeatureDB
        The parsed annotation database.
    """
    db_file = gtf_file + ".db"
    if not os.path.exists(db_file):
        gffutils.create_db(
            gtf_file,
            dbfn=db_file,
            force=True,
            keep_order=True,
            merge_strategy='merge',
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
    return gffutils.FeatureDB(db_file)

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


def load_gtf_features(path: str, feature_type) -> pd.DataFrame:
    """Load intervals of specific feature type(s) from a GTF.

    Parameters
    ----------
    path : str
        Path to a GTF (GENCODE/Ensembl style).
    feature_type : str or tuple of str
        Feature type(s) to load (e.g. "exon", "CDS", ("UTR",
        "five_prime_UTR", "three_prime_UTR")). A tuple loads all
        matching types in one pass.

    Returns
    -------
    pandas.DataFrame
        Columns: chrom, start (0-based), end, gene_name, strand.
        One row per feature.
    """
    types = {feature_type} if isinstance(feature_type, str) \
        else set(feature_type)
    rows = []
    gene_re = re.compile(r'gene_name "([^"]+)"')
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in types:
                continue
            m = gene_re.search(fields[8])
            rows.append((fields[0], int(fields[3]) - 1, int(fields[4]),
                         m.group(1) if m else fields[8][:20], fields[6]))
    return pd.DataFrame(rows, columns=["chrom", "start", "end",
                                       "gene_name", "strand"])

def load_bw_mean_coverage(bw_path: Optional[str]) -> float:
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


def load_onco_genes(path: Optional[str]) -> set:
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


def load_tpm(path: Optional[str]) -> Dict[str, float]:
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


def fragment_mean_coverage(bw_path: Optional[str], chrom: str,
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


def plot_gene_model_all_transcripts(
        ax:plt.Axes,
        gene: gffutils.Feature,
        transcripts: Dict[str, Dict[str, List[Tuple[int,int]]]],
    ) -> plt.Axes:
    """Plot a gene model showing all transcripts with optional mutation overlay.

    Draws exons (sky-blue rectangles), UTRs (green rectangles), intron
    connectors, a strand arrow, and coloured mutation/SV markers.  The
    figure is saved to *output* at 300 DPI in all requested formats.

    Parameters
    ----------
    gene : gffutils.Feature
        The gene feature, used for coordinate range and strand info.
    transcripts : dict
        Mapping of transcript ID to ``{'exons': [...], 'UTRs': [...]}``.
    mutations : list of dict, optional
        Each dict must contain ``pos`` (int) and ``type`` (str).
    figsize : tuple of float, optional
        Figure width and height in inches, by default ``(10, 2)``.
    output : str, optional
        Base path for the saved figure (without extension), by default
        ``"gene_model"``.
    """
    n_transcripts = len(transcripts)
    y_gap = 1.0 / max(n_transcripts, 1)
    
    type_colors = {"DEL":"#0AF0A7", "DUP":"#F00A0AFA", "INS":"#0A9CF0",
                   "INV":"#9FF00A", "BND":"#E4F00A", "TRA": "#F0660A"}

    for i, (tx_id, tx_data) in enumerate(transcripts.items()):
        y_base = i * y_gap
        exons = tx_data['exons']
        utrs = tx_data['UTRs']

        # Draw introns as lines
        for j in range(len(exons)-1):
            ax.plot([exons[j][1], exons[j+1][0]], [y_base, y_base], color='black', lw=1, zorder=0)

        # Exons
        for start, end in exons:
            ax.add_patch(patches.Rectangle((start, y_base-0.1), end-start, 0.2, color='skyblue', ec='red', zorder=2))
        
        # UTRs
        for start, end in utrs:
            ax.add_patch(patches.Rectangle((start, y_base-0.05), end-start, 0.1, color='lightgreen', ec='green', zorder=3))
        
        # Transcript label
        ax.text(gene.start-50, y_base, tx_id, fontsize=8, ha='right', va='center')

    # Strand arrow
    arrow_gap = (gene.end-gene.start)/20
    if gene.strand == '+':
        ax.arrow(gene.start, -0.1, gene.end-gene.start, 0, head_width=0.05, head_length=arrow_gap,
                 length_includes_head=True, color='black')
    elif gene.strand == '-':
        ax.arrow(gene.end, -0.1, gene.start-gene.end, 0, head_width=0.05, head_length=arrow_gap,
                 length_includes_head=True, color='black')

    ax.legend(handles=[patches.Patch(color=c, label=t) for t,c in type_colors.items()], bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlim(gene.start-100, gene.end+100)
    ax.set_ylim(-0.2, 1.2)
    ax.set_yticks([])
    ax.set_xlabel(f"{gene.chrom}:{gene.start}-{gene.end}")
    ax.set_title(f"Gene model: {gene.attributes.get('gene_name',[''])[0]}")
    return ax

def to_ang(chrom: str, g0: int, g1: int, spans: List[Tuple], circ_len: int) -> Optional[Tuple[float, float]]:
    """Map a genomic interval to circular angles (t2, t1) in degrees.

    pos 0 = 12 o'clock, clockwise; ``90 - 360*pos/CIRC_LEN``.
    Returns ``(t2, t1)`` where t2 <= t1 (ccw arc for matplotlib Wedge).
    """
    for fc, fs, fe, xs, xe, _st in spans:
        if fc != chrom:
            continue
        s, e = max(g0, fs), min(g1, fe)
        if e <= s:
            continue
        a0 = 90.0 - 360.0 * (s - fs + xs) / circ_len
        a1 = 90.0 - 360.0 * (e - fs + xs) / circ_len
        return a1, a0  # (t2, t1) ccw
    return None


def _arc_contains(start: float, span: float, ang: float) -> bool:
    """Test whether an angle lies inside a circular arc.

    Angles are in degrees with wraparound at 360.

    Parameters
    ----------
    start : float
        Arc start angle (degrees).
    span : float
        Arc span (degrees, non-negative).
    ang : float
        Angle to test (degrees).

    Returns
    -------
    bool
        True when ``ang`` falls within [start, start + span] (mod 360).
    """
    return (ang - start) % 360.0 <= span


def _ang_overlap(t2a, t1a, t2b, t1b):
    """True when two ccw arcs [t2, t1] overlap (with wraparound)."""
    span_a, span_b = t1a - t2a, t1b - t2b
    sa, sb = t2a % 360.0, t2b % 360.0
    return (_arc_contains(sa, span_a, sb)
            or _arc_contains(sa, span_a, (sb + span_b) % 360.0)
            or _arc_contains(sb, span_b, sa)
            or _arc_contains(sb, span_b, (sa + span_a) % 360.0))

def _mark_break(ax: plt.Axes, r_lane: float, e2: float, e1: float,
                ex_start: int, ex_end: int, col: str, exon_w: float,
                frag_edges: List[Tuple[int, int]], g0: int, g1: int,
                truncated: bool) -> None:
    """Add ``//`` break gaps when an exon touches a clipping boundary.

    A gene is truncated when it is not fully contained in one fragment;
    exons (or CDS blocks) ending exactly at a fragment clipping boundary
    get the double-slash break symbol to mark the cycle junction.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    r_lane : float
        Gene lane outer radius.
    e2, e1 : float
        Block arc angles (degrees, ccw from e2 to e1).
    ex_start, ex_end : int
        Genomic start/end of the exon or CDS block.
    col : str
        Gene colour (MAGENTA intact / TRUNCATED truncated).
    exon_w : float
        Radial width of the exon block.
    frag_edges : list of (int, int)
        (fragment_start, fragment_end) pairs of the gene's chromosome.
    g0, g1 : int
        Genomic start/end of the parent gene.
    truncated : bool
        Whether the gene is truncated (not contained in one fragment).
    """
    if not truncated:
        return
    for fs_e, fe_e in frag_edges:
        lo, hi = max(g0, fs_e), min(g1, fe_e)
        if hi <= lo:
            continue
        if abs(int(ex_start) - lo) < 2:
            _draw_break_gap_at(ax, r_lane, e2, col, exon_w)
        if abs(int(ex_end) - hi) < 2:
            _draw_break_gap_at(ax, r_lane, e1, col, exon_w)


def _draw_gene_structure(
    ax: plt.Axes,
    gene: pd.Series,
    transcripts: Dict[str, Dict[str, List[Tuple[int, int]]]],
    t2: float,
    t1: float,
    lane: int,
    contained: bool,
    spans: List[Tuple],
    circ_len: int,
    onco: bool,
) -> Tuple[float, Optional[float]]:
    """Draw one gene's structure on the circular map.

    Mirrors the signature of :func:`plot_gene_model_all_transcripts`:
    a single gene plus its ``transcripts`` structure dict in, drawing
    only that gene (one transcript per gene here; the dict is kept for
    structural parity).  Renders the intron connector arcs (dashed for
    truncated genes), the exon blocks (light full-span), the CDS
    blocks (solid full-height), the UTR blocks (half-height), the
    ``//`` break symbols where the gene is clipped at a cycle
    junction, and the 3' strand arrowhead.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    gene : pandas.Series
        The gene interval row (chrom, start, end, gene_name, strand).
    transcripts : dict
        Mapping of transcript ID to ``{'exons': [...], 'CDS': [...],
        'UTRs': [...]}`` in genomic coordinates; each entry a list of
        ``(start, end)`` tuples (0-based, half-open).  May contain
        empty lists when the GTF has no CDS/UTR annotation.
    t2, t1 : float
        Gene arc angles in degrees (ccw arc from t2 to t1).
    lane : int
        Radial lane index (0 = outermost gene lane).
    contained : bool
        Whether the gene lies fully inside one fragment (intact gene).
    spans : list of tuple
        Circular coordinate spans (chrom, f_start, f_end, x_start,
        x_end, strand) of the cycle fragments.
    circ_len : int
        Total circular coordinate length (sum of fragment lengths).
    onco : bool
        Whether this gene is an oncogene highlight (dark exon outline).

    Returns
    -------
    r_lane : float
        Outer radius of the lane the gene was drawn on.
    lab_ang : float or None
        Suggested label angle in degrees (the gene-start angle); None
        when the arc is too short to carry a label.
    """
    strand = str(gene.get("strand", "+"))
    r_lane = R_GENE - lane * LANE_W
    col = MAGENTA if contained else TRUNCATED
    g0, g1 = int(gene["start"]), int(gene["end"])

    # single-transcript rendering: merge the dict's feature lists
    tx = next(iter(transcripts.values())) if transcripts \
        else {"exons": [], "CDS": [], "UTRs": []}
    exons = tx.get("exons", [])
    cdss = tx.get("CDS", [])
    utrs = tx.get("UTRs", [])

    # intron connectors: thin arcs ONLY between consecutive exons
    # (gene_model.py convention: intron line from exon[j].end to
    # exon[j+1].start); dashed for truncated genes
    r_bb = r_lane - EXON_W / 2.0 + 0.007
    r_line = r_bb + 0.007
    if not exons:
        # no exon annotation: single backbone across the gene span
        if contained:
            ax.add_patch(Wedge((0, 0), r_line, t2, t1,
                               width=0.014,
                               facecolor=col, edgecolor="none",
                               alpha=0.9, zorder=4))
        else:
            ax.add_patch(Wedge((0, 0), r_line, t2, t1,
                               width=0.014,
                               facecolor="none", edgecolor=col,
                               linewidth=0.9,
                               linestyle=(0, (2.5, 1.8)),
                               alpha=0.9, zorder=4))
    else:
        spans_sorted = sorted(exons)
        for j in range(len(spans_sorted) - 1):
            gap0 = int(spans_sorted[j][1])
            gap1 = int(spans_sorted[j + 1][0])
            ghit = to_ang(gene["chrom"], gap0, gap1, spans, circ_len)
            if ghit is None:
                continue
            g2, g1a = ghit
            if contained:
                ax.add_patch(Wedge((0, 0), r_line, g2, g1a,
                                   width=0.014,
                                   facecolor=col, edgecolor="none",
                                   alpha=0.9, zorder=4))
            else:
                ax.add_patch(Wedge((0, 0), r_line, g2, g1a,
                                   width=0.014,
                                   facecolor="none", edgecolor=col,
                                   linewidth=0.9,
                                   linestyle=(0, (2.5, 1.8)),
                                   alpha=0.9, zorder=4))

    # gene model: exons as light full blocks, CDS full-height dark,
    # UTR half-height; colour = intact magenta / truncated orange
    frag_edges = [(fs, fe) for fc, fs, fe, _xs, _xe, _st in spans
                  if fc == gene["chrom"]]

    if not exons:
        # no exon annotation: draw the whole gene as one block
        ax.add_patch(Wedge((0, 0), r_lane, t2, t1, width=EXON_W,
                           facecolor=col, edgecolor="white",
                           linewidth=0.25, zorder=5))
        if not contained:
            _draw_break_gaps(ax, r_lane, t1, t2, col, EXON_W)
    else:
        for ex_start, ex_end in sorted(exons):
            ehit = to_ang(gene["chrom"], int(ex_start), int(ex_end), spans, circ_len)
            if ehit is None:
                continue
            e2, e1 = ehit
            ecol = "#3A3A70" if onco else "white"
            ax.add_patch(Wedge(
                (0, 0), r_lane, e2, e1, width=EXON_W,
                facecolor=col, edgecolor=ecol,
                linewidth=0.9 if onco else 0.25,
                alpha=0.35, zorder=5))
            _mark_break(ax, r_lane, e2, e1, ex_start, ex_end,
                        col, EXON_W, frag_edges, g0, g1, not contained)
        # UTR: half-height blocks centered in the exon band
        for ut_start, ut_end in sorted(utrs):
            uhit = to_ang(gene["chrom"], int(ut_start), int(ut_end), spans, circ_len)
            if uhit is None:
                continue
            u2, u1 = uhit
            ax.add_patch(Wedge(
                (0, 0), r_lane - EXON_W * 0.275, u2, u1,
                width=EXON_W * 0.45,
                facecolor=col, edgecolor="none", alpha=0.85,
                zorder=6))
        # CDS: full-height solid blocks
        for cd_start, cd_end in sorted(cdss):
            chit = to_ang(gene["chrom"], int(cd_start), int(cd_end), spans, circ_len)
            if chit is None:
                continue
            c2, c1 = chit
            ax.add_patch(Wedge(
                (0, 0), r_lane, c2, c1, width=EXON_W,
                facecolor=col, edgecolor="white",
                linewidth=0.25, zorder=6))
            _mark_break(ax, r_lane, c2, c1, cd_start, cd_end,
                        col, EXON_W, frag_edges, g0, g1, not contained)

    # strand arrowhead: small tangential arrow at the end of the
    # intron line (3' end), sized to stay visible at 300 dpi
    tip_ang = t2 if strand == "+" else t1
    rad = np.radians(tip_ang)
    r_tip = r_bb + 0.007  # center line of the intron band
    p0 = np.array([r_tip * np.cos(rad), r_tip * np.sin(rad)])
    tang = np.array([-np.sin(rad), np.cos(rad)]) * \
        (0.022 if strand == "+" else -0.022)
    norm = np.array([np.cos(rad), np.sin(rad)]) * 0.014
    ax.add_patch(plt.Polygon([p0 + tang, p0 + norm, p0 - norm],
                             closed=True, facecolor=col, zorder=6))

    # label angle: gene-start angle, None when the arc is too short
    start_ang = t1 if strand == "+" else t2
    arc_len = abs(t2 - t1)
    return r_lane, (start_ang if arc_len >= 1.2 else None)


def _draw_gene_ring(
    ax: plt.Axes,
    genes: pd.DataFrame,
    circ_id: str,
    spans: List[Tuple],
    circ_len: int,
    onco_genes: set,
    tpm: Optional[Dict[str, float]],
    bw_path: Optional[str],
    exons_df: pd.DataFrame,
    cds_df: pd.DataFrame,
    utr_df: pd.DataFrame,
) -> Tuple[List[dict], bool, bool]:
    """Draw the gene-structure ring: lanes, exon/CDS/UTR blocks, labels.

    Coordination layer: draws the fragment join ticks, maps each gene
    to its circular arc, assigns non-overlapping radial lanes greedily
    (longest arc first), assembles each gene's ``transcripts`` dict
    (same structure as :func:`plot_gene_model_all_transcripts` takes)
    from the exon/CDS/UTR interval tables, delegates the per-gene
    rendering to :func:`_draw_gene_structure`, and adds the radial
    gene labels (TPM appended for oncogenes, collision-avoided against
    previously placed labels).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    genes : pandas.DataFrame
        Gene intervals (see :func:`load_gtf_features`); empty = no genes.
    circ_id : str
        Cycle identifier used in the per-gene table rows.
    spans : list of tuple
        Circular coordinate spans (chrom, f_start, f_end, x_start,
        x_end, strand) of the cycle fragments.
    circ_len : int
        Total circular coordinate length (sum of fragment lengths).
    onco_genes : set
        Uppercased gene symbols to highlight.
    tpm : dict or None
        gene_name -> TPM; appended to oncogene labels.
    bw_path : str or None
        Coverage bigwig for per-gene median coverage; None = skip.
    exons_df, cds_df, utr_df : pandas.DataFrame
        Exon / CDS / UTR intervals (see :func:`load_gtf_features`).

    Returns
    -------
    gene_rows : list of dict
        Per-gene annotation rows (circ_id, gene, status, onco, strand,
        region, median_coverage, tpm).
    has_truncated : bool
        Whether any truncated gene was drawn.
    has_onco : bool
        Whether any oncogene-highlighted gene was drawn.
    """
    # Fragment join ticks: short coral marks just outside the outermost
    # gene lane at every fragment start/end (the circular join points).
    join_angles: List[float] = []
    for fc, fs, fe, xs, xe, st in spans:
        for pos in (xs, xe):
            join_angles.append(90.0 - 360.0 * pos / circ_len)
    for ja in join_angles:
        ra = np.radians(ja)
        rr0 = R_GENE + LANE_W * 0.6
        rr1 = rr0 - 0.045
        ax.plot([rr0 * np.cos(ra), rr1 * np.cos(ra)],
                [rr0 * np.sin(ra), rr1 * np.sin(ra)],
                color=CORAL, linewidth=1.6, zorder=4)

    gene_rows: List[dict] = []
    if genes is None or genes.empty:
        return gene_rows, False, False

    frag_chroms = {fc for fc, *_ in spans}
    gpool = genes[genes["chrom"].isin(frag_chroms)].copy()

    # Build gene records with circular angles
    recs: List[dict] = []
    for _, g in gpool.iterrows():
        hit = to_ang(g["chrom"], int(g["start"]), int(g["end"]), spans, circ_len)
        if hit is None:
            continue
        t2, t1 = hit
        # determine if gene is fully contained within a fragment
        contained = False
        for fc, fs, fe, xs, xe, _st in spans:
            if fc != g["chrom"]:
                continue
            if int(g["start"]) >= fs and int(g["end"]) <= fe:
                contained = True
                break
        recs.append({"g": g, "t2": t2, "t1": t1,
                     "contained": contained})

    # longest arc first: t2-t1 is negative (ccw), most negative = longest
    recs.sort(key=lambda r: r["t2"] - r["t1"])


    # lane assignment: greedy, longest-first
    placed: List[dict] = []
    for rec in recs:
        lane = 0
        while True:
            clash = any(o["lane"] == lane and _ang_overlap(
                rec["t2"], rec["t1"], o["t2"], o["t1"])
                for o in placed)
            if not clash:
                break
            lane += 1
        rec["lane"] = lane
        placed.append(rec)

    # draw each gene via the single-gene renderer, then add its label
    has_truncated = False
    has_onco = False
    tpm_map = tpm or {}
    for rec in placed:
        g = rec["g"]
        gname = str(g["gene_name"])
        strand = str(g.get("strand", "+"))
        onco = gname.upper() in onco_genes
        has_truncated = has_truncated or not rec["contained"]
        has_onco = has_onco or onco

        # this gene's transcripts dict (same structure as
        # plot_gene_model_all_transcripts takes); one transcript per
        # gene, keyed by the gene name
        exons = exons_df[(exons_df["chrom"] == g["chrom"])
                         & (exons_df["gene_name"] == gname)]
        cdss = cds_df[(cds_df["chrom"] == g["chrom"])
                      & (cds_df["gene_name"] == gname)]
        utrs = utr_df[(utr_df["chrom"] == g["chrom"])
                      & (utr_df["gene_name"] == gname)]
        transcripts: Dict[str, Dict[str, List[Tuple[int, int]]]] = {
            gname: {
                "exons": [(int(r["start"]), int(r["end"]))
                          for _, r in exons.iterrows()],
                "CDS": [(int(r["start"]), int(r["end"]))
                        for _, r in cdss.iterrows()],
                "UTRs": [(int(r["start"]), int(r["end"]))
                         for _, r in utrs.iterrows()],
            },
        }

        _, lab_ang = _draw_gene_structure(
            ax, g, transcripts, rec["t2"], rec["t1"], rec["lane"],
            rec["contained"], spans, circ_len, onco)

        # gene label at gene-start angle, colored same as arc, no leader line
        if lab_ang is not None:
            placed_angles = [o["lab_ang"] for o in placed
                             if "lab_ang" in o]
            la = lab_ang
            for pa in placed_angles:
                d = abs((la - pa + 180) % 360 - 180)
                if d < 20:
                    la = (pa - 22) % 360  # nudge apart
                    break
            rec["lab_ang"] = la
            rrad = np.radians(la)
            rr = R_RULER + 0.12
            col = MAGENTA if rec["contained"] else TRUNCATED
            label_text = gname
            if onco and gname in tpm_map:
                label_text += f" | {tpm_map[gname]:.1f}"
            ax.text(rr * np.cos(rrad), rr * np.sin(rrad), label_text,
                    fontsize=5.6, color=col,
                    ha="center", va="center", zorder=5,
                    fontweight="bold" if onco else "normal",
                    rotation=(np.degrees(rrad) + 90) % 180 - 90)

        # collect per-gene data row
        cov_arr = fragment_mean_coverage(bw_path, g["chrom"],
                                         int(g["start"]), int(g["end"]),
                                         bins=20)
        dosage = (float(np.median(cov_arr))
                  if cov_arr is not None and len(cov_arr) else np.nan)
        gene_rows.append({
            "circ_id": circ_id,
            "gene": gname,
            "status": "truncated" if not rec["contained"] else "intact",
            "onco": onco,
            "strand": strand,
            "region": f"{g['chrom']}:{int(g['start'])}-{int(g['end'])}",
            "median_coverage": round(dosage, 1) if pd.notna(dosage) else "",
            "tpm": tpm_map.get(gname, "") if tpm else "",
        })

    return gene_rows, has_truncated, has_onco


def draw_cycle_circular(
    ax: plt.Axes,
    fragments: pd.DataFrame,
    info: Optional[pd.Series],
    circ_id: str,
    genes: pd.DataFrame,
    onco_genes: set,
    bw_path: Optional[str],
    tpm: Optional[Dict[str, float]],
    bw_mean_cov: float,
    exons_df: Optional[pd.DataFrame],
    cds_df: Optional[pd.DataFrame] = None,
    utr_df: Optional[pd.DataFrame] = None,
) -> Optional[pd.DataFrame]:
    """Draw one ecDNA cycle as a CGView-style circular genome map.

    Renders the ecDNA as a circular figure with:

    * outer tick ruler (kb-scale on the circular coordinate)
    * ring 1  coverage profile (gold area, log-scale radial height)
    * ring 2  fragment boundary markers (coral dashed lines at breakpoints)
    * ring 3  gene blocks (magenta/violet), labels radiating outward
    * center  hub with circ_id / topology / size / coordinates

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Single axes on which to draw the circular map.
    fragments : pandas.DataFrame
        Fragment rows of one circ_id, kept in BED file order.
    info : pandas.Series or None
        Matching summary.txt row (size_mb, topology, proportion).
    circ_id : str
        Cycle identifier used in the hub text.
    genes : pandas.DataFrame
        Gene intervals (see :func:`load_gtf_features`); empty = no annotation.
    onco_genes : set
        Uppercased gene symbols to highlight.
    bw_path : str, optional
        Coverage bigwig; None = skip the coverage ring.
    tpm : dict, optional
        gene_name -> TPM; appended to oncogene labels.
    bw_mean_cov : float
        Genome-wide mean coverage (used as log-scale reference floor).
    exons_df : pandas.DataFrame, optional
        Exon intervals (see :func:`load_gtf_features`).
    cds_df : pandas.DataFrame, optional
        CDS intervals; drawn as full-height coding blocks inside exons.
    utr_df : pandas.DataFrame, optional
        UTR intervals; drawn as half-height blocks inside exons.

    Returns
    -------
    pandas.DataFrame or None
        Per-gene table (gene, status, dosage, tpm) for this cycle;
        None when no genes overlap.
    """
    if exons_df is None:
        exons_df = pd.DataFrame(
            columns=["chrom", "start", "end", "gene_name", "strand"])
    if cds_df is None:
        cds_df = pd.DataFrame(
            columns=["chrom", "start", "end", "gene_name", "strand"])
    if utr_df is None:
        utr_df = pd.DataFrame(
            columns=["chrom", "start", "end", "gene_name", "strand"])

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.75, 1.55)

    # ---- circular coordinate: concatenate fragments end to end ----
    spans: List[Tuple] = []  # (chrom, f_start, f_end, x_start, x_end, strand)
    cursor = 0
    for _, frag in fragments.iterrows():
        f_start, f_end = int(frag["start"]), int(frag["end"])
        length = f_end - f_start
        spans.append((frag["chrom"], f_start, f_end,
                      cursor, cursor + length, str(frag.get("strand", "+"))))
        cursor += length
    circ_len = cursor
    if circ_len <= 0:
        return None

    for kb in range(0, int(circ_len / 1000) + 1, 10):
        angle_deg = 90.0 - 360.0 * kb * 1000 / circ_len
        a = np.radians(angle_deg)

        long_tick = kb % 50 == 0
        r0 = R_RULER
        tick_len = 0.045 if long_tick else 0.025
        r1 = r0 - tick_len

        ax.plot(
            [r0 * np.cos(a), r1 * np.cos(a)],
            [r0 * np.sin(a), r1 * np.sin(a)],
            color=INDIGO,
            linewidth=2.2 if long_tick else 1.4,
            solid_capstyle="butt",
        )

        # 主刻度标签距离更远，次刻度标签距离更近
        label_r = r0 + (0.055 if long_tick else 0.045)

        # 文字沿圆周切线方向旋转，并保持正向可读
        rotation = angle_deg - 90.0
        rotation = (rotation + 180.0) % 360.0 - 180.0
        if rotation > 90.0:
            rotation -= 180.0
        elif rotation < -90.0:
            rotation += 180.0

        ax.text(
            label_r * np.cos(a),
            label_r * np.sin(a),
            f"{kb} kb",
            fontsize=6 if long_tick else 5,
            color=INDIGO,
            ha="center",
            va="center",
            rotation=rotation,
            rotation_mode="anchor",
        )


    # ================= ring 1: coverage (log-scale radial) =============
    if bw_path and os.path.isfile(bw_path):
        try:
            bwf = pyBigWig.open(bw_path)
            n = 600
            cov = np.zeros(n)
            for i in range(n):
                pos = circ_len * i / n
                for fc, fs, fe, xs, xe, _st in spans:
                    if xs <= pos < xe:
                        g = fs + (pos - xs)
                        bin_size = max(int(circ_len // n), 1)
                        v = bwf.stats(fc, int(g), int(g) + bin_size,
                                      type="mean")
                        cov[i] = v[0] if v and v[0] else 0.0
                        break
            bwf.close()
            cmax = float(cov.max()) if cov.max() else 1.0
            # log scale with floor: 0.5x -> ~12% height, genome mean -> ~45%
            ref = max(float(bw_mean_cov), 0.3)
            h = np.clip(np.log10(np.maximum(cov, 0.3)) /
                        max(np.log10(max(cmax, ref)), 1e-6), 0.06, 1.0)
            # angular frame matching to_ang: th = pi/2 - 2*pi*i/n
            th = np.pi / 2 - 2 * np.pi * np.arange(n) / n
            r_out = R_COV_IN + (R_COV_OUT - R_COV_IN) * h
            ox = r_out * np.cos(th)
            oy = r_out * np.sin(th)
            ix = R_COV_IN * np.cos(th[::-1])
            iy = R_COV_IN * np.sin(th[::-1])
            ax.add_patch(plt.Polygon(
                np.column_stack([np.append(ox, ix), np.append(oy, iy)]),
                facecolor=GOLD, alpha=0.75, edgecolor="none", zorder=2))
            ax.plot(np.append(ox, ox[0]), np.append(oy, oy[0]),
                    color=GOLD_D, linewidth=0.7, zorder=3)
        except Exception as exc:  # noqa: BLE001
            ax.text(0, 0.85, f"coverage unavailable ({exc})",
                    fontsize=6, color=GRAY, ha="center")

    # ================= ring 2: genes (radial lanes + exons) ===========
    gene_rows, has_truncated, has_onco = _draw_gene_ring(
        ax, genes, circ_id, spans, circ_len, onco_genes, tpm,
        bw_path, exons_df, cds_df, utr_df)

    # ================= hub + legend =================
    _draw_hub(ax, info, circ_id, spans, circ_len)
    _draw_legend(ax, has_truncated=has_truncated, has_onco=has_onco)

    return pd.DataFrame(gene_rows) if gene_rows else None


def _draw_hub(ax: plt.Axes, info: Optional[pd.Series], circ_id: str,
              spans: List[Tuple], circ_len: int) -> None:
    """Draw the central hub circle with cycle metadata."""
    ax.add_patch(plt.Circle((0, 0), 0.40, facecolor="white",
                            edgecolor="#D5DBE1", linewidth=1.0, zorder=6))
    ax.text(0, 0.10, f"circ {circ_id}", ha="center", fontsize=11,
            color=INDIGO, fontweight="bold", zorder=7)
    topology = str(info["topology"]) if info is not None else "NA"
    ax.text(0, -0.02, topology, ha="center", fontsize=7,
            color=GRAY, zorder=7)
    ax.text(0, -0.14, f"{circ_len / 1000:.1f} kb", ha="center", fontsize=8,
            color=GRAY, zorder=7)
    # coordinate summary from first/last fragment
    if spans:
        first = spans[0]
        last = spans[-1]
        if first[0] == last[0]:
            coord_text = f"{first[0]}:{first[1]:,}-{last[2]:,}"
        else:
            coord_text = f"{first[0]}:{first[1]:,}…{last[0]}:{last[2]:,}"
        ax.text(0, -0.26, coord_text, ha="center",
                fontsize=5.5, color=GRAY, zorder=7)



def _draw_break_gaps(ax: plt.Axes, r_lane: float,
                     t1: float, t2: float, col: str, exon_w: float,
                     n_gaps: int = 2) -> None:
    """Draw double-slash break gaps at both ends of a truncated gene block.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    r_lane : float
        Gene lane outer radius.
    t1, t2 : float
        Block arc angles (degrees, ccw from t2 to t1).
    col : str
        Colour of the slash strokes.
    exon_w : float
        Radial width of the exon block.
    n_gaps : int
        Number of slash strokes per gap.
    """
    for ta in (t1, t2):
        _draw_break_gap_at(ax, r_lane, ta, col, exon_w, n_gaps)


def _draw_break_gap_at(ax: plt.Axes, r_lane: float, ta: float,
                       col: str, exon_w: float,
                       n_gaps: int = 2) -> None:
    """Draw a double-slash break gap at one arc angle end.

    Two short parallel strokes crossing the exon block at angle ``ta``
    (each rotated 45 deg to the local tangent), the standard "//" break
    symbol used for truncated features.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    r_lane : float
        Gene lane outer radius.
    ta : float
        Arc end angle (degrees).
    col : str
        Stroke colour.
    exon_w : float
        Radial width of the exon block.
    n_gaps : int
        Number of slash strokes.
    """
    rad = np.radians(ta)
    r_mid = r_lane - exon_w / 2.0
    cx, cy = r_mid * np.cos(rad), r_mid * np.sin(rad)
    # local tangent direction (clockwise) and radial direction
    tang = np.array([-np.sin(rad), np.cos(rad)])
    norm = np.array([np.cos(rad), np.sin(rad)])
    # slash strokes: slightly longer than the block half-width, offset
    # along the tangent so they read as "//" rather than "||"
    slash_len = exon_w * 0.75
    for i in range(n_gaps):
        off = (i - (n_gaps - 1) / 2.0) * 0.008
        # each slash: 45 deg between tangent and radial
        d = (tang + norm) / np.sqrt(2.0)
        p0 = np.array([cx, cy]) + tang * off - d * slash_len / 2.0
        p1 = np.array([cx, cy]) + tang * off + d * slash_len / 2.0
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
                color=col, linewidth=1.2, solid_capstyle="round",
                zorder=7)


def _legend_block(ax: plt.Axes, x: float, y: float, w: float,
                  h_frac: float, col: str) -> None:
    """Draw one exon/CDS/UTR block of the miniature legend gene model.

    Block height scales with ``h_frac``: 1.0 = full CDS height,
    0.5 = half-height UTR.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x : float
        Left edge of the block.
    y : float
        Centre line of the model.
    w : float
        Block width.
    h_frac : float
        Height as a fraction of the full block height (1.0 = CDS).
    col : str
        Fill colour.
    """
    ax.add_patch(plt.Rectangle(
        (x, y - 0.026 * h_frac), w, 0.052 * h_frac,
        facecolor=col, edgecolor="white", linewidth=0.25, zorder=8))


def _draw_legend_gene(ax: plt.Axes, col: str, x0: float, y: float,
                      truncated: bool) -> None:
    """Draw a miniature gene model as one legend entry.

    Mirrors the on-circle rendering exactly: light exon blocks, solid
    CDS full-height blocks, half-height UTR blocks, thin intron line,
    3' arrowhead; truncated adds dashed intron and // break symbols.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    col : str
        Gene colour (MAGENTA intact / TRUNCATED truncated).
    x0, y : float
        Left edge and centre line of the sample.
    truncated : bool
        Draw the truncated variant (dashed intron + // breaks).
    """
    line_style: Dict = dict(color=col, linewidth=0.9, zorder=8)
    if truncated:
        line_style["linestyle"] = (0, (2.5, 1.8))
    # intron line spanning the model
    ax.plot([x0, x0 + 0.10], [y, y], **line_style)

    # layout: [UTR][CDS][intron][CDS][UTR] over the 0.10 span
    xs = [0.000, 0.018, 0.040, 0.062, 0.082, 0.100]
    _legend_block(ax, x0 + xs[0], y, xs[1] - xs[0], 0.5, col)  # UTR half-height
    _legend_block(ax, x0 + xs[1], y, xs[2] - xs[1], 1.0, col)  # CDS full
    _legend_block(ax, x0 + xs[3], y, xs[4] - xs[3], 1.0, col)  # CDS full
    _legend_block(ax, x0 + xs[4], y, xs[5] - xs[4], 0.5, col)  # UTR half-height
    if truncated:
        # // break symbols at both outer ends (gene cut by the junction)
        for xx, direction in ((x0, -1.0), (x0 + 0.100, 1.0)):
            d = np.array([0.7071, 0.7071]) * direction
            for k in (-1, 1):
                p0 = np.array([xx + 0.005 * k, y])
                p1 = p0 + d * 0.024
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
                        color=col, linewidth=1.1,
                        solid_capstyle="round", zorder=8)
    # 3' arrowhead at the right end (tail of the model gene)
    tip = np.array([x0 + 0.106, y])
    tang = np.array([0.030, 0.0])
    norm = np.array([0.0, 0.022])
    ax.add_patch(plt.Polygon(
        [tip + tang, tip - norm, tip + norm],
        closed=True, facecolor=col, zorder=8))


def _draw_legend(ax: plt.Axes, has_truncated: bool = True,
                 has_onco: bool = False) -> None:
    """Draw a frameless legend outside the circular map.

    Entries are full miniature gene models (not abstract line swatches)
    so every legend item maps 1:1 onto what is drawn on the circle.
    Generated dynamically: the truncated-gene entry only appears when
    the cycle contains truncated genes, the oncogene entry only when
    highlighted genes are present.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    has_truncated : bool
        Whether truncated genes occur in this cycle.
    has_onco : bool
        Whether oncogene-highlighted blocks occur in this cycle.
    """
    lx0, ly0 = 0.80, -1.02
    entries = [("patch", GOLD, "read coverage (log scale)")]
    if has_onco:
        entries.append(("gene", INDIGO, "oncogene (dark outline)"))
    entries.append(("gene", MAGENTA, "intact gene"))
    if has_truncated:
        entries.append(("gene_trunc", TRUNCATED,
                        "truncated gene (// = cycle junction)"))
    entries.append(("patch", CORAL, "cycle junction (fragment join)"))

    for i, (kind, col, txt) in enumerate(entries):
        yy = ly0 - i * 0.13
        if kind == "patch":
            ax.add_patch(plt.Rectangle(
                (lx0, yy - 0.018), 0.07, 0.036,
                facecolor=col, edgecolor="none", zorder=8))
        elif kind == "gene":
            _draw_legend_gene(ax, col, lx0, yy, truncated=False)
        elif kind == "gene_trunc":
            _draw_legend_gene(ax, col, lx0, yy, truncated=True)
        ax.text(
            lx0 - 0.10, yy, txt,
            fontsize=5.8, color="#3B4A5A",
            va="center", ha="right", zorder=8)


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
    fig, ax = plt.subplots(figsize=(6, 6))
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

def _draw_and_collect(
        fig: plt.Figure, cid: str, bed: pd.DataFrame,
        summary: pd.DataFrame, genes: pd.DataFrame, onco_genes: set,
        bw_path: Optional[str], tpm: Dict[str, float],
        bw_mean_cov: float, exons_df: Optional[pd.DataFrame],
        cds_df: Optional[pd.DataFrame],
        utr_df: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
    """Draw one cycle on a full-figure axes and return its gene table.

    Allocates an axes covering the whole figure, delegates rendering to
    :func:`draw_cycle_circular`, and hands back the per-gene annotation
    rows for aggregation by the caller.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to add the axes to.
    cid : str
        Cycle identifier to draw.
    bed : pandas.DataFrame
        Fragment-level rows (see :func:`load_ecdna_bed`).
    summary : pandas.DataFrame
        Summary rows indexed by circ_id (see :func:`load_summary`).
    genes : pandas.DataFrame
        Gene intervals (see :func:`load_gtf_features`).
    onco_genes : set
        Uppercased gene symbols to highlight.
    bw_path : str or None
        Coverage bigwig; None = skip the dosage ring.
    tpm : dict
        gene_name -> TPM.
    bw_mean_cov : float
        Genome-wide mean coverage (log-scale reference floor).
    exons_df, cds_df, utr_df : pandas.DataFrame or None
        Exon / CDS / UTR intervals (see :func:`load_gtf_features`).

    Returns
    -------
    pandas.DataFrame or None
        Per-gene table for this cycle; None when no genes overlap.
    """
    frags = bed[bed["circ_id"].astype(str) == cid]
    info = summary.loc[cid] if cid in summary.index else None
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    return draw_cycle_circular(
        ax, frags, info, cid, genes, onco_genes,
        bw_path, tpm, bw_mean_cov, exons_df, cds_df, utr_df)


def _grid_page(
        page_ids: List[str], page: int, n_pages: int, n_cycles: int,
        suffix: str, bed: pd.DataFrame, summary: pd.DataFrame,
        genes: pd.DataFrame, onco_genes: set, bw_path: Optional[str],
        tpm: Dict[str, float], bw_mean_cov: float,
        sample_id: str, out_prefix: str, image_formats: List[str],
        exons_df: Optional[pd.DataFrame],
        cds_df: Optional[pd.DataFrame],
        utr_df: Optional[pd.DataFrame]) -> List[pd.DataFrame]:
    """Render one grid page of cycle structures to image files.

    Lays out up to 3 cycles per row on a 9x9 figure, delegates each
    subplot to :func:`draw_cycle_circular`, adds the sample suptitle and
    saves in every requested format.

    Parameters
    ----------
    page_ids : list of str
        Cycle identifiers on this page.
    page, n_pages : int
        Zero-based page index and total page count (for the suptitle).
    n_cycles : int
        Total cycle count across all pages (for the suptitle).
    suffix : str
        Filename suffix between ``out_prefix`` and the extension
        ("" for a single page, ".p1", ".p2", ... otherwise).
    bed : pandas.DataFrame
        Fragment-level rows (see :func:`load_ecdna_bed`).
    summary : pandas.DataFrame
        Summary rows indexed by circ_id (see :func:`load_summary`).
    genes : pandas.DataFrame
        Gene intervals (see :func:`load_gtf_features`).
    onco_genes : set
        Uppercased gene symbols to highlight.
    bw_path : str or None
        Coverage bigwig; None = skip the dosage ring.
    tpm : dict
        gene_name -> TPM.
    bw_mean_cov : float
        Genome-wide mean coverage (log-scale reference floor).
    sample_id : str
        Sample label used in the suptitle.
    out_prefix : str
        Output path prefix; ``<suffix>.<fmt>`` appended per format.
    image_formats : list of str
        Image formats (png/pdf/...).
    exons_df, cds_df, utr_df : pandas.DataFrame or None
        Exon / CDS / UTR intervals (see :func:`load_gtf_features`).

    Returns
    -------
    list of pandas.DataFrame
        Per-gene tables of the cycles on this page (may be empty).
    """
    fig = plt.figure(figsize=(9, 9))
    gene_rows: List[pd.DataFrame] = []
    for k, cid in enumerate(page_ids):
        # Allocate an independent axes for each cycle.
        ncols = min(3, len(page_ids))
        nrows = math.ceil(len(page_ids) / ncols)
        ax = fig.add_subplot(nrows, ncols, k + 1)
        frags = bed[bed["circ_id"].astype(str) == cid]
        info = summary.loc[cid] if cid in summary.index else None

        rows = draw_cycle_circular(
            ax, frags, info, cid, genes, onco_genes,
            bw_path, tpm, bw_mean_cov, exons_df, cds_df, utr_df)

        if rows is not None:
            gene_rows.append(rows)

    fig.suptitle(
        f"{sample_id}: ecDNA structures "
        f"({n_cycles} cycles; page {page + 1}/{n_pages})",
        fontsize=10)

    for fmt in image_formats:
        path = f"{out_prefix}{suffix}.{fmt}"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        logger.info("Saved %s", path)

    plt.close(fig)
    return gene_rows


def plot_structures(bed: pd.DataFrame, summary: pd.DataFrame,
                    genes: pd.DataFrame, onco_genes: set,
                    bw_path: Optional[str], tpm: Dict[str, float],
                    sample_id: str, out_prefix: str,
                    image_formats: List[str], max_per_fig: int,
                    one_file_per_cycle: bool = False,
                    exons_df: Optional[pd.DataFrame] = None,
                    cds_df: Optional[pd.DataFrame] = None,
                    utr_df: Optional[pd.DataFrame] = None) -> None:
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
        Gene intervals (see :func:`load_gtf_features`); empty = no annotation.
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
    all_gene_rows: List[pd.DataFrame] = []
    bw_mean_cov = load_bw_mean_coverage(bw_path)

    if bw_path:
        logger.info(f"Genome-wide mean coverage: {bw_mean_cov:.2f}x "
                    "(dosage ring reference)")

    if one_file_per_cycle:
        # One ecDNA, one figure: single large square canvas per cycle.
        # The combined grid page is still emitted afterwards so downstream
        # rules keep a single summary image.
        os.makedirs(os.path.dirname(os.path.abspath(out_prefix)),
                    exist_ok=True)
        for cid in circ_ids:
            fig = plt.figure(figsize=(9, 9))
            rows = _draw_and_collect(
                fig, cid, bed, summary, genes, onco_genes,
                bw_path, tpm, bw_mean_cov, exons_df, cds_df, utr_df)
            if rows is not None:
                all_gene_rows.append(rows)
            for fmt in image_formats:
                path = f"{out_prefix}_circ{cid}.{fmt}"
                fig.savefig(path, dpi=DPI, bbox_inches="tight")
                logger.info(f"Saved {path}")
            plt.close(fig)
        n_pages = (len(circ_ids) + max_per_fig - 1) // max_per_fig
        if n_pages == 1:
            all_gene_rows.extend(_grid_page(
                circ_ids, 0, 1, len(circ_ids), "", bed, summary,
                genes, onco_genes, bw_path, tpm, bw_mean_cov,
                sample_id, out_prefix, image_formats,
                exons_df, cds_df, utr_df))
        else:
            for page in range(n_pages):
                all_gene_rows.extend(_grid_page(
                    circ_ids[page * max_per_fig:
                             (page + 1) * max_per_fig],
                    page, n_pages, len(circ_ids), f".p{page + 1}",
                    bed, summary, genes, onco_genes, bw_path, tpm,
                    bw_mean_cov, sample_id, out_prefix, image_formats,
                    exons_df, cds_df, utr_df))
    else:
        n_pages = (len(circ_ids) + max_per_fig - 1) // max_per_fig
        if n_pages == 1:
            all_gene_rows.extend(_grid_page(
                circ_ids, 0, 1, len(circ_ids), "", bed, summary,
                genes, onco_genes, bw_path, tpm, bw_mean_cov,
                sample_id, out_prefix, image_formats,
                exons_df, cds_df, utr_df))
        else:
            for page in range(n_pages):
                all_gene_rows.extend(_grid_page(
                    circ_ids[page * max_per_fig:
                             (page + 1) * max_per_fig],
                    page, n_pages, len(circ_ids), f".p{page + 1}",
                    bed, summary, genes, onco_genes, bw_path, tpm,
                    bw_mean_cov, sample_id, out_prefix, image_formats,
                    exons_df, cds_df, utr_df))

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
    cds = None
    utrs = None
    if args.gtf:
        if os.path.isfile(args.gtf):
            genes = load_gtf_features(args.gtf, "gene")
            logger.info(f"Loaded {len(genes)} genes from {args.gtf}")
            exons = load_gtf_features(args.gtf, "exon")
            logger.info(f"Loaded {len(exons)} exons from {args.gtf}")
            cds = load_gtf_features(args.gtf, "CDS")
            logger.info(f"Loaded {len(cds)} CDS from {args.gtf}")
            utrs = load_gtf_features(
                args.gtf, ("UTR", "five_prime_UTR", "three_prime_UTR"))
            logger.info(f"Loaded {len(utrs)} UTR from {args.gtf}")
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
                    exons_df=exons, cds_df=cds, utr_df=utrs)


if __name__ == "__main__":
    main()
