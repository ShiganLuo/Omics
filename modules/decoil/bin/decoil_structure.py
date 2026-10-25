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
import colorsys
from dataclasses import dataclass
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
# ---------------------------------------------------------------------------
# Radial band allocation.  RingBands is the single source of truth for the
# figure layout: every ring owns one radial band, bands are strictly
# ordered (hub -> topology -> gene -> coverage -> ruler -> gene labels),
# and every ring-drawing function receives its band edges through this
# object so rings can never overlap.  Instance fields carry defaults;
# construct RingBands(...) with keyword overrides to retune the layout.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RingBands:
    """Radial band allocation for one circular cycle figure.

    Each ring occupies the annulus between its inner and outer radius;
    a validation check in ``__post_init__`` guarantees the bands are
    disjoint and ordered from the hub outward.  Derived properties
    expose the per-ring geometry (lane pitch, exon height, lane count,
    micro-gene dot radius) so no draw function needs module globals.

    Parameters
    ----------
    hub : float
        Hub disc radius (cycle metadata text).
    topo_in, topo_out : float
        eccDNA topology motif band.
    gene_in, gene_out : float
        Gene structure band; gene lanes stack inward from ``gene_out``.
    lane_pitch : float
        Radial pitch between gene lanes; also the exon block height.
    cov_in, cov_out : float
        Coverage profile band (log-scale radial height).
    ruler : float
        Ruler radius (ticks point inward toward ``cov_out``).
    label_in : float
        Gene-label band inner edge; labels step outward on collision.
    label_pitch : float
        Radial pitch between label collision rings.

    Raises
    ------
    ValueError
        When the bands are not strictly ordered hub-outward.
    """
    hub: float = 0.34
    topo_in: float = 0.345
    topo_out: float = 0.40
    gene_in: float = 0.45
    gene_out: float = 0.68
    lane_pitch: float = 0.075
    cov_in: float = 0.70
    cov_out: float = 0.85
    ruler: float = 0.95
    label_in: float = 1.05
    label_pitch: float = 0.09

    def __post_init__(self) -> None:
        edges = (self.hub, self.topo_in, self.topo_out, self.gene_in,
                 self.gene_out, self.cov_in, self.cov_out, self.ruler,
                 self.label_in)
        if any(a >= b for a, b in zip(edges, edges[1:])):
            raise ValueError(
                f"ring bands not strictly ordered: {edges}")

    @property
    def exon_h(self) -> float:
        """Radial height of an exon block (one full lane pitch)."""
        return self.lane_pitch

    @property
    def max_lanes(self) -> int:
        """Number of gene lanes that fit between gene_in and gene_out."""
        return max(1, int(round((self.gene_out - self.gene_in)
                                / self.lane_pitch)))

    @property
    def dot_r(self) -> float:
        """Micro-gene dot radius: the free corridor inside the gene band."""
        return (self.topo_out + self.gene_in) / 2.0

    def lane_outer(self, lane: int) -> float:
        """Outer radius of a gene lane (lane 0 = outermost)."""
        return self.gene_out - lane * self.lane_pitch


INTRON_W = 0.006        # intron line radial width (~8% of exon height)
MIN_GENE_DEG = 0.15     # genes with total arc below this render as dots
TOPO_GAP_DEG = 0.8      # white excision-gap width (simple_circle) — as
                        # narrow as a breakpoint should read
TOPO_TICK_DEG = 3.0     # orange SV tick arc width (simple_sv/mixed_sv)
TOPO_EPS_DEG = 0.25     # flush-join overlap to hide antialias hairlines

# palette mirrors topology_legend.png
TOPO_BASE = "#5B7DA3"   # backbone ring
TOPO_SV = "#F28E2B"     # SV / deletion junction ticks
TOPO_INV = "#59A14F"    # inversion arc (single arrow tip)
TOPO_DUP_HI = "#C2474F" # first occurrence of a duplicated region
TOPO_DUP_LO = "#EFB9BE" # repeated occurrence(s) (A / A' pair)
TOPO_FB_A = "#F5C08A"   # foldback cone, entry side
TOPO_FB_B = "#E8862B"   # foldback cone, return side

# reference-style palette
CORAL = "#E8734F"       # fragment boundary markers
GOLD = "#E3B34D"        # coverage fill
GOLD_D = "#B98A2F"      # coverage outline
MAGENTA = "#8E3A6B"     # intact gene blocks (deep, high contrast)
INDIGO = "#3A3A70"      # ticks / text
GRAY = "#8A94A0"        # secondary text
TRUNCATED = "#C97B4A"   # legacy truncated tint (kept for fallback)

# hues assigned during the current render — gene_color probes against
# this registry so every gene on a ring gets a distinct colour; reset
# per cycle by calling reset_gene_hues()
GENE_HUES: List[float] = []
# name -> colour cache: a gene's model and its label both call
# gene_color(name); the cache guarantees one single colour per gene
# and keeps the hue registry free of duplicates
GENE_COLOR_CACHE: Dict[str, str] = {}


def reset_gene_hues() -> None:
    """Clear the assigned-hue registry before rendering a new cycle."""
    GENE_HUES.clear()
    GENE_COLOR_CACHE.clear()


def gene_color(name: str) -> str:
    """Return a stable, unique-per-name colour.

    Colour = hue from the name hash, adjusted so that hues already
    taken by other genes in this render are avoided (linear probe for
    the largest hue distance); saturation/value are fixed at levels
    that survive white outlines and dark labels.  Uniqueness holds
    within a render because the caller registers each assigned hue in
    ``GENE_HUES``; across renders the same name always yields the same
    colour (hash + same probe sequence over the same registry state).

    Parameters
    ----------
    name : str
        Gene symbol (any stable identifier works).

    Returns
    -------
    str
        Hex colour string.
    """
    h = 0
    for ch in name:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFF
    hue = (h % 360) / 360.0
    if name in GENE_COLOR_CACHE:
        return GENE_COLOR_CACHE[name]
    if GENE_HUES:
        best_hue, best_d = hue, -1.0
        n = len(GENE_HUES)
        for i in range(4 * n):
            cand = (hue + i * 0.5 / n) % 1.0
            d = min(abs((cand - t + 0.5) % 1.0 - 0.5) for t in GENE_HUES)
            if d > best_d:
                best_d, best_hue = d, cand
        hue = best_hue
    GENE_HUES.append(hue)
    r, g, b = colorsys.hsv_to_rgb(hue, 0.55, 0.62)
    colour = "#%02X%02X%02X" % (int(r * 255), int(g * 255), int(b * 255))
    GENE_COLOR_CACHE[name] = colour
    return colour


DPI = 300

# legend anchor (data coords): the bottom-right corner of the canvas
# (xlim -1.55..1.80, ylim -1.30..1.55); text runs right-aligned up to
# the swatch, swatch right edge at 1.75 — flush with the canvas corner
LEGEND_XY = (1.68, -1.20)

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


def load_gtf_all_features(path: str) -> Dict[str, pd.DataFrame]:
    """Load gene/exon/UTR intervals from a GTF in a single pass.

    Reads the GTF file once and splits the feature types needed by the
    structure plots.  CDS is deliberately not loaded: the plots show
    exon/intron/UTR gene structure, not coding-sequence detail.

    Parameters
    ----------
    path : str
        Path to a GTF (GENCODE/Ensembl style).

    Returns
    -------
    dict
        Keys "genes", "exons", "utrs" -> DataFrames with columns
        chrom, start (0-based), end, gene_name, strand.  UTRs collect
        ``UTR``, ``five_prime_UTR`` and ``three_prime_UTR`` rows.
    """
    type_map = {"gene": "genes", "exon": "exons",
                "UTR": "utrs", "five_prime_UTR": "utrs",
                "three_prime_UTR": "utrs"}
    buckets: Dict[str, List[Tuple]] = {"genes": [], "exons": [],
                                       "utrs": []}
    gene_re = re.compile(r'gene_name "([^"]+)"')
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            key = type_map.get(fields[2])
            if key is None:
                continue
            m = gene_re.search(fields[8])
            buckets[key].append(
                (fields[0], int(fields[3]) - 1, int(fields[4]),
                 m.group(1) if m else fields[8][:20], fields[6]))
    cols = ["chrom", "start", "end", "gene_name", "strand"]
    return {k: pd.DataFrame(v, columns=cols) for k, v in buckets.items()}

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


def _resolve_label_angle(la: float, occupied: List[Tuple[float, int]],
                         min_sep: float = 20.0,
                         max_iter: int = 40) -> Tuple[float, int]:
    """Keep the label at its gene-start angle; resolve collisions radially.

    The label must stay aligned with the gene's 5' start, so the angle
    is never moved.  Each occupant is an ``(angle, ring)`` pair; when
    the desired angle falls inside ``min_sep`` of an occupant on the
    same radial ring, the label moves one ring outward.  Occupants on
    different rings never clash — one radial step (0.10 R ≈ 80 px)
    is far more than a text row height.

    Parameters
    ----------
    la : float
        Gene-start label angle (degrees).
    occupied : list of (float, int)
        ``(angle, ring)`` of residents — other gene labels and ruler
        ticks (ticks at ring 0).
    min_sep : float
        Angular window within which an occupant on the same ring
        forces a radial step.
    max_iter : int
        Maximum radial steps.

    Returns
    -------
    tuple (float, int)
        ``(angle, ring_index)`` — angle unchanged; ring 0 is the
        innermost label radius.
    """
    la = la % 360.0
    ring = 0
    for _ in range(max_iter):
        # clash rules:
        #   - same ring + angle window < min_sep  → one ring outward
        #   - adjacent ring (|Δring| == 1) + angle window < 5°   → also
        #     outward: a rotated label's text width projects far enough
        #     radially to touch the neighbouring ring at near-equal angles
        clash = False
        for pa, pr in occupied:
            d = abs((la - pa + 180.0) % 360.0 - 180.0)
            if pr == ring and d < min_sep:
                clash = True
                break
            if abs(pr - ring) == 1 and d < 5.0:
                clash = True
                break
        if not clash:
            break
        ring += 1
    return la, ring


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
    bands: RingBands,
    clip_edges: Optional[List[Tuple[int, int]]] = None,
) -> Tuple[float, Optional[float]]:
    """Draw one gene's structure on the circular map.

    Mirrors the signature of :func:`plot_gene_model_all_transcripts`:
    a single gene plus its ``transcripts`` structure dict in, drawing
    only that gene (one transcript per gene here; the dict is kept for
    structural parity).  Renders the intron connector arcs (solid,
    centered in the exon band), the solid exon blocks, the paler
    half-height UTR blocks, the ``//`` break symbols where the gene is
    clipped at a cycle junction, and the 3' strand arrowhead.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    gene : pandas.Series
        The gene interval row (chrom, start, end, gene_name, strand).
    transcripts : dict
        Mapping of transcript ID to ``{'exons': [...],
        'UTRs': [...]}`` in genomic coordinates; each entry a list of
        ``(start, end)`` tuples (0-based, half-open).  May contain
        empty lists when the GTF has no UTR annotation.
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
    bands : RingBands
        Radial band allocation; supplies the lane radii and exon height.
    clip_edges : list of (int, int), optional
        ``(fragment_start, fragment_end)`` genomic boundaries of the
        cycle fragments on this gene's chromosome.  A gene-span end
        gets a ``//`` break symbol only when it lies within 2 bp of one
        of these boundaries — i.e. only where the circle actually cut
        the gene.  An end far from any boundary (gene reaching across
        the ring start back into intact sequence) is not a breakpoint
        and must not be marked.  None/empty = skip break symbols.

    Returns
    -------
    r_lane : float
        Outer radius of the lane the gene was drawn on.
    lab_ang : float or None
        Suggested label angle in degrees (the gene-start angle); None
        when the arc is too short to carry a label.
    """
    strand = str(gene.get("strand", "+"))
    r_lane = bands.lane_outer(lane)
    exon_h = bands.exon_h
    # per-gene colour: distinct hue per gene, keyed by gene name so the
    # same gene keeps its colour across cycles of the same figure.
    # Truncation is signalled by the // break marks, not by hue.
    col = gene_color(str(gene.get("gene_name", gene.name)))

    # single-transcript rendering: merge the dict's feature lists
    tx = next(iter(transcripts.values())) if transcripts \
        else {"exons": [], "UTRs": []}
    exons = tx.get("exons", [])
    utrs = tx.get("UTRs", [])

    # intron connectors: dashed thin arcs between consecutive exons,
    # centered exactly in the exon band so the line touches the exon
    # blocks at both ends
    r_line = r_lane - exon_h / 2.0
    if not exons:
        # no exon annotation: single backbone across the gene span
        ax.add_patch(Wedge((0, 0), r_line, t2, t1,
                           width=INTRON_W,
                           facecolor="none", edgecolor=col,
                           linewidth=0.9,
                           linestyle=(0, (2.5, 1.8)),
                           zorder=4))
    else:
        spans_sorted = sorted(exons)
        for j in range(len(spans_sorted) - 1):
            gap0 = int(spans_sorted[j][1])
            gap1 = int(spans_sorted[j + 1][0])
            ghit = to_ang(gene["chrom"], gap0, gap1, spans, circ_len)
            if ghit is None:
                continue
            g2, g1a = ghit
            ax.add_patch(Wedge((0, 0), r_line, g2, g1a,
                               width=0,
                               facecolor="none", edgecolor=col,
                               linewidth=0.9,
                               linestyle=(0, (6, 4)),
                               zorder=4))

    # gene model: exons as solid full-height blocks; UTR portions inside
    # each exon are overlaid paler (same band, same width), the standard
    # browser convention: UTR reads as the pale ends of the exon block
    if not exons:
        # no exon annotation: draw the whole gene as one block
        ax.add_patch(Wedge((0, 0), r_lane, t2, t1, width=exon_h,
                           facecolor=col, edgecolor="white",
                           linewidth=0.25, zorder=5))
    else:
        # exons render as full-height rectangles; a gene whose whole
        # visible arc is too small to show blocks (below MIN_GENE_DEG)
        # renders as a single dot at its midpoint so tiny genes still
        # occupy visible space instead of vanishing next to their name
        gene_sweep = abs(t1 - t2)
        if gene_sweep < MIN_GENE_DEG:
            mid = np.radians((t1 + t2) / 2.0)
            r_mid_dot = r_lane - exon_h / 2.0
            ax.plot([r_mid_dot * np.cos(mid)],
                    [r_mid_dot * np.sin(mid)],
                    marker="o", markersize=4.5, color=col,
                    markeredgecolor="white", markeredgewidth=0.8,
                    zorder=6)
            return r_lane, ((t1 if strand == "+" else t2)
                            if gene_sweep >= 0.2 else None)
        for ex_start, ex_end in sorted(exons):
            ehit = to_ang(gene["chrom"], int(ex_start), int(ex_end), spans, circ_len)
            if ehit is None:
                continue
            e2, e1 = ehit
            ecol = "#3A3A70" if onco else "white"
            ax.add_patch(Wedge(
                (0, 0), r_lane, e2, e1, width=exon_h,
                facecolor=col, edgecolor=ecol,
                linewidth=0.9 if onco else 0.25,
                zorder=5))
        # UTR: paler overlay on the exon band, same radius and width —
        # reads as the non-coding ends of the exon blocks
        for ut_start, ut_end in sorted(utrs):
            uhit = to_ang(gene["chrom"], int(ut_start), int(ut_end), spans, circ_len)
            if uhit is None:
                continue
            u2, u1 = uhit
            ax.add_patch(Wedge(
                (0, 0), r_lane, u2, u1, width=exon_h,
                facecolor=col, edgecolor="white", alpha=0.35,
                linewidth=0.3, zorder=6))

    # ``//`` break symbols: only on a truncated gene, and only at the
    # gene-span end(s) that lie on a real fragment clipping boundary.
    # A truncated gene whose far end wraps across the ring start sits
    # on intact reference sequence there — that end gets no break mark.
    # A gene clipped *inside* its span (fragment ends mid-gene) has its
    # visible span end exactly at the fragment boundary: the effective
    # clipped end is min/max(gene end, fragment end), so compare the
    # span-end angle against that effective coordinate.
    if not contained and clip_edges:
        # to_ang semantics: (t2, t1) = (high-coord end, low-coord end),
        # independent of strand — t2 always faces the gene's genomic end
        g_start, g_end = int(gene["start"]), int(gene["end"])
        ends = [(t1, g_start),   # low-coordinate end
                (t2, g_end)]     # high-coordinate end
        # a gene whose span reaches past a fragment boundary is clipped
        # there: treat the fragment boundary itself as the visible end
        eff_ends = []
        for ang, gpos in ends:
            hit = any(abs(gpos - fs) < 2 or abs(gpos - fe) < 2
                      for fs, fe in clip_edges)
            if hit:
                eff_ends.append(ang)
            else:
                # is this end beyond a fragment boundary (clipped mid-gene)?
                for fs, fe in clip_edges:
                    if fs < gpos <= g_end and gpos > fe:
                        eff_ends.append(ang)   # 3' side runs past fragment end
                        break
                    if g_end >= fs and gpos < fs:
                        eff_ends.append(ang)   # start side before fragment start
                        break
        for ang in eff_ends:
            _draw_break_gap_at(ax, r_lane, ang, col, exon_h)

    # strand direction is expressed by the tapered first exon above —
    # no separate arrow element (a fixed-size arrow overwhelmed short
    # genes: Gm13222's arrow was 20× the gene's own arc length)

    # label angle: for short genes the 5' start angle (name stays next
    # to the gene origin); for long genes (arc > 90°) the arc MIDPOINT,
    # because a 100+ kb gene's model stretches across a large part of
    # the ring and a name pinned to one end sits far from the model's
    # visual centre, leaving the reader unable to match name to model
    start_ang = t1 if strand == "+" else t2
    arc_len = abs(t2 - t1)
    if arc_len > 90.0:
        mid_ang = (t1 + t2) / 2.0 if t1 > t2 else ((t1 + t2) / 2.0 + 180.0)
        start_ang = mid_ang
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
    utr_df: pd.DataFrame,
    bands: RingBands,
    tick_angles: Optional[List[float]] = None,
) -> Tuple[List[dict], bool, bool]:
    """Draw the gene-structure ring: lanes, exon/UTR blocks, labels.

    Coordination layer: maps each gene to its circular arc, assigns
    non-overlapping radial lanes greedily (longest arc first),
    assembles each gene's ``transcripts`` dict (same structure as
    :func:`plot_gene_model_all_transcripts` takes) from the exon/UTR
    interval tables, delegates the per-gene rendering to
    :func:`_draw_gene_structure`, and adds the radial gene labels (TPM
    appended for oncogenes, collision-avoided against previously
    placed labels).  The fragment joins are visible as the white gaps
    in the coverage ring and the ``//`` break symbols; no separate
    join marks are drawn.

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
    exons_df, utr_df : pandas.DataFrame
        Exon / UTR intervals (see :func:`load_gtf_features`).

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
    gene_rows: List[dict] = []
    if genes is None or genes.empty:
        return gene_rows, False, False

    frag_chroms = {fc for fc, *_ in spans}
    gpool = genes[genes["chrom"].isin(frag_chroms)].copy()

    # per-chromosome fragment clipping boundaries (genomic coords):
    # a gene end within 2 bp of one of these is really cut by the circle
    frag_edges_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for fc, fs, fe, _xs, _xe, _st in spans:
        frag_edges_by_chrom.setdefault(fc, []).append((int(fs), int(fe)))

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

    # micro genes (whole-gene arc below MIN_GENE_DEG) leave the gene
    # ring entirely: their dot + name move to the outer mark band just
    # inside the ruler.  On the gene lanes they would either vanish
    # (sub-pixel) or sit on top of the host gene covering them.
    micro_recs: List[dict] = []
    lane_recs: List[dict] = []
    for rec in recs:
        if abs(rec["t1"] - rec["t2"]) < MIN_GENE_DEG:
            micro_recs.append(rec)
        else:
            lane_recs.append(rec)
    tpm_map = tpm or {}

    # ---- pass 1: lane genes (models + labels) -----------------------
    # lane-gene label angles register into lane_occupied so that the
    # micro-gene names placed afterwards avoid them too
    max_lane = bands.max_lanes
    placed: List[dict] = []
    micro_rows: List[dict] = []
    for rec in lane_recs:
        lane = 0
        while lane < max_lane:
            clash = any(o["lane"] == lane and _ang_overlap(
                rec["t2"], rec["t1"], o["t2"], o["t1"])
                for o in placed)
            if not clash:
                break
            lane += 1
        rec["lane"] = lane
        placed.append(rec)

    # ---- pass 2: lane gene models + labels (angles register first) ----
    has_truncated = False
    has_onco = False
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
        utrs = utr_df[(utr_df["chrom"] == g["chrom"])
                      & (utr_df["gene_name"] == gname)]
        transcripts: Dict[str, Dict[str, List[Tuple[int, int]]]] = {
            gname: {
                "exons": [(int(r["start"]), int(r["end"]))
                          for _, r in exons.iterrows()],
                "UTRs": [(int(r["start"]), int(r["end"]))
                         for _, r in utrs.iterrows()],
            },
        }

        _, lab_ang = _draw_gene_structure(
            ax, g, transcripts, rec["t2"], rec["t1"], rec["lane"],
            rec["contained"], spans, circ_len, onco, bands,
            clip_edges=frag_edges_by_chrom.get(g["chrom"]))

        # gene label pinned to the gene-start (5') angle; collisions are
        # resolved radially (outward rings) so the name stays aligned
        # with its gene origin
        if lab_ang is not None:
            # residents: gene labels (their own ring) + ruler ticks (ring 0)
            occupied = ([(o["lab_ang"], o.get("label_ring", 0))
                         for o in placed if "lab_ang" in o]
                        + [(ta, 0) for ta in (tick_angles or [])])
            la, ring = _resolve_label_angle(lab_ang, occupied)
            rec["lab_ang"] = la
            rec["label_ring"] = ring
            rrad = np.radians(la)
            rr = bands.label_in + ring * bands.label_pitch
            rec["label_r"] = rr
            rec["label_x"] = rr * np.cos(rrad)
            rec["label_y"] = rr * np.sin(rrad)
            col = gene_color(gname)
            label_text = gname
            if onco and gname in tpm_map:
                label_text += f" | {tpm_map[gname]:.1f}"
            ax.text(rr * np.cos(rrad), rr * np.sin(rrad), label_text,
                    fontsize=5.6, color=col,
                    ha="center", va="center", zorder=5,
                    fontweight="bold" if onco else "normal",
                    # labels on radial rings (ring>0) stay horizontal:
                    # a rotated long name projects its full text length
                    # onto the radial direction and collides with the
                    # neighbouring ring; horizontal text only occupies
                    # its glyph height radially
                    rotation=(0 if ring > 0
                              else (np.degrees(rrad) + 90) % 180 - 90))

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
            "label_x": rec.get("label_x"),
            "label_y": rec.get("label_y"),
        })

    # ---- pass 3: micro genes (dot + name in a dedicated radial band) --
    # micro genes own a fixed radial band strictly between the gene
    # lanes and the ruler; lane genes never enter it, so the two can
    # never crowd each other regardless of angles:
    #   dot band : free corridor between topology band and gene band
    #   name band: [label_in, label_in + n*label_pitch]  (outside ruler)
    for rec in micro_recs:
        g = rec["g"]
        gname = str(g["gene_name"])
        onco = gname.upper() in onco_genes
        col = gene_color(gname)
        mid_ang = ((rec["t1"] + rec["t2"]) / 2.0) % 360.0
        occupied = ([(o["lab_ang"], o.get("label_ring", 0))
                     for o in placed if "lab_ang" in o]
                    + [(ta, 0) for ta in (tick_angles or [])])
        la, micro_ring = _resolve_label_angle(mid_ang, occupied)
        rec["lab_ang"] = la
        rec["label_ring"] = micro_ring
        rrad = np.radians(la)
        # dot radius: the free corridor between the topology band and
        # the gene band — a fixed zone no other element uses
        r_dot = bands.dot_r
        r_lab = bands.label_in + micro_ring * bands.label_pitch
        ax.plot([r_dot * np.cos(rrad)], [r_dot * np.sin(rrad)],
                marker="o", markersize=4.0, color=col,
                markeredgecolor="white", markeredgewidth=0.8, zorder=7)
        ax.text(r_lab * np.cos(rrad), r_lab * np.sin(rrad), gname,
                fontsize=5.6, color=col,
                ha="center", va="center", zorder=5,
                fontweight="bold" if onco else "normal",
                rotation=(np.degrees(rrad) + 90) % 180 - 90)
        cov_arr = fragment_mean_coverage(bw_path, g["chrom"],
                                         int(g["start"]), int(g["end"]),
                                         bins=20)
        dosage = (float(np.median(cov_arr))
                  if cov_arr is not None and len(cov_arr) else np.nan)
        micro_rows.append({
            "circ_id": circ_id,
            "gene": gname,
            "status": "small_gene",
            "onco": onco,
            "strand": str(g.get("strand", "+")),
            "region": f"{g['chrom']}:{int(g['start'])}-{int(g['end'])}",
            "median_coverage": round(dosage, 1) if pd.notna(dosage) else "",
            "tpm": tpm_map.get(gname, "") if tpm else "",
        })

    gene_rows.extend(micro_rows)
    return gene_rows, has_truncated, has_onco


def _draw_ruler_ring(ax: plt.Axes, spans: List[Tuple], circ_len: int,
                     bands: RingBands) -> List[float]:
    """Draw the coordinate ruler ring and return major-tick angles.

    Adaptive scale: the major step is the largest of the preferred
    steps giving <= 12 labels around the ring; minor ticks subdivide
    each major interval 5:1.  Ticks point inward from ``bands.ruler``
    toward the coverage band; major ticks carry tangential kb labels
    just outside the ruler.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    spans : list of tuple
        Circular coordinate spans (unused for geometry, kept for
        signature parity with the other ring functions).
    circ_len : int
        Total circular coordinate length.
    bands : RingBands
        Radial band allocation (ruler radius).

    Returns
    -------
    list of float
        Major-tick angles in degrees — gene labels must avoid these.
    """
    total_kb = circ_len / 1000.0
    for step in (5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000):
        if total_kb / step <= 12:
            major_step = step
            break
    else:
        major_step = 5000
    minor_step = major_step / 5.0

    # major-tick label angles — gene labels must not land on these
    tick_angles = [(90.0 - 360.0 * k / total_kb) % 360.0
                   for k in np.arange(0, total_kb + 1e-9, major_step)]

    r0 = bands.ruler
    for k in np.arange(0, total_kb + minor_step * 0.1, minor_step):
        # skip a trailing minor tick that would crowd the 0 kb major
        if k > total_kb - minor_step / 2.0 and k > 0:
            continue
        angle_deg = 90.0 - 360.0 * k / total_kb
        a = np.radians(angle_deg)

        # major when k is (numerically) a multiple of major_step
        long_tick = abs(k / major_step - round(k / major_step)) < 1e-6
        tick_len = 0.045 if long_tick else 0.025
        r1 = r0 - tick_len

        ax.plot(
            [r0 * np.cos(a), r1 * np.cos(a)],
            [r0 * np.sin(a), r1 * np.sin(a)],
            color=INDIGO,
            linewidth=2.2 if long_tick else 1.4,
            solid_capstyle="butt",
        )

        # only major ticks carry labels
        if not long_tick:
            continue
        label_r = r0 + 0.055

        # 文字沿圆周切线方向旋转，并保持正向可读
        rotation = angle_deg - 90.0
        rotation = (rotation + 180.0) % 360.0 - 180.0
        if rotation > 90.0:
            rotation -= 180.0
        elif rotation < -90.0:
            rotation += 180.0

        label = f"{int(round(k))} kb" if major_step >= 1 else f"{k:.1f} kb"
        ax.text(
            label_r * np.cos(a),
            label_r * np.sin(a),
            label,
            fontsize=6,
            color=INDIGO,
            ha="center",
            va="center",
            rotation=rotation,
            rotation_mode="anchor",
        )
    return tick_angles


def _draw_coverage_ring(ax: plt.Axes, spans: List[Tuple], circ_len: int,
                        bw_path: Optional[str], bw_mean_cov: float,
                        bands: RingBands,
                        topology: str = "") -> None:
    """Draw the coverage profile ring (log-scale radial height).

    Samples the bigwig evenly along the circular coordinate and fills
    the band between ``bands.cov_in`` and a log-scaled outer radius
    capped at ``bands.cov_out``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    spans : list of tuple
        Circular coordinate spans of the cycle fragments.
    circ_len : int
        Total circular coordinate length.
    bw_path : str or None
        Coverage bigwig; None/missing = draw nothing.
    bw_mean_cov : float
        Genome-wide mean coverage (log-scale reference floor).
    bands : RingBands
        Radial band allocation (coverage band edges).
    topology : str
        topology_name of the cycle.  Backbone topologies carry the
        shared white excision wedge at the coordinate origin; flush-
        joined topologies (multi_region/duplications) keep their
        seamless colour join there instead.
    """
    if not bw_path or not os.path.isfile(bw_path):
        return
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
        r_out = bands.cov_in + (bands.cov_out - bands.cov_in) * h
        ox = r_out * np.cos(th)
        oy = r_out * np.sin(th)
        ix = bands.cov_in * np.cos(th[::-1])
        iy = bands.cov_in * np.sin(th[::-1])
        ax.add_patch(plt.Polygon(
            np.column_stack([np.append(ox, ix), np.append(oy, iy)]),
            facecolor=GOLD, alpha=0.75, edgecolor="none", zorder=2))
        ax.plot(np.append(ox, ox[0]), np.append(oy, oy[0]),
                color=GOLD_D, linewidth=0.7, zorder=3)
        # breakpoint slits: a cycle of n fragments has exactly n
        # breakpoints — the wrap-around join at the coordinate origin
        # (12 o'clock) plus the n-1 inter-fragment junctions.  Cut a
        # white slit of the shared TOPO_GAP_DEG width at each, so every
        # junction reads identically in every ring (the topology band
        # marks the same angles with its gap/orange ticks/colour joins).
        for slit_ang in ([90.0]
                         + [_t1 for _t2, _t1 in _frag_angles(spans,
                                                             circ_len)[1:]]):
            ax.add_patch(Wedge((0, 0), bands.cov_out + 0.02,
                               slit_ang - TOPO_GAP_DEG / 2.0,
                               slit_ang + TOPO_GAP_DEG / 2.0,
                               width=bands.cov_out - bands.cov_in + 0.02,
                               facecolor="white", edgecolor="none",
                               zorder=4))
    except Exception as exc:  # noqa: BLE001
        ax.text(0, bands.cov_out + 0.05, f"coverage unavailable ({exc})",
                fontsize=6, color=GRAY, ha="center")


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
    utr_df: Optional[pd.DataFrame] = None,
    bands: Optional[RingBands] = None,
) -> Optional[pd.DataFrame]:
    """Draw one ecDNA cycle as a CGView-style circular genome map.

    Renders the ecDNA as a circular figure with:

    * outer tick ruler (kb-scale on the circular coordinate)
    * ring 1  coverage profile (gold area, log-scale radial height)
    * ring 2  gene blocks (magenta/violet) with intron connectors;
      fragment joins visible as coverage-ring gaps and // breaks
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
    if utr_df is None:
        utr_df = pd.DataFrame(
            columns=["chrom", "start", "end", "gene_name", "strand"])
    if bands is None:
        bands = RingBands()

    ax.set_aspect("equal")
    ax.axis("off")
    # canvas spans the ruler ring, gene labels (labels step outward
    # from bands.label_in) and the right-side legend (text to x≈1.68)
    ax.set_xlim(-1.55, 1.80)
    ax.set_ylim(-1.30, 1.55)

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

    # ---- ring dispatch (hub -> topology -> gene -> coverage -> ruler) --
    # every ring draws strictly inside its own RingBands band, so the
    # call order below is purely cosmetic; the band edges guarantee the
    # rings never overlap
    tick_angles = _draw_ruler_ring(ax, spans, circ_len, bands)
    _draw_coverage_ring(ax, spans, circ_len, bw_path, bw_mean_cov, bands,
                        topology=str(info["topology"])
                        if info is not None else "")
    gene_rows, has_truncated, has_onco = _draw_gene_ring(
        ax, genes, circ_id, spans, circ_len, onco_genes, tpm,
        bw_path, exons_df, utr_df, bands, tick_angles=tick_angles)
    _draw_topology_ring(
        ax, str(info["topology"]) if info is not None else "",
        spans, circ_len, bands)
    _draw_hub(ax, info, circ_id, spans, circ_len, bands)
    # legend pinned flush to the bottom-right canvas corner; gene labels
    # near the corner can only push it DOWN (grow the canvas), never
    # left — the corner position itself is fixed
    if gene_rows:
        bottom_labels = [r for r in gene_rows
                         if r.get("label_y") is not None
                         and r["label_y"] < -0.4
                         and r["label_x"] > 0.3]
        lowest_y = (min(r["label_y"] for r in bottom_labels)
                    if bottom_labels else None)
    else:
        lowest_y = None
    lx0, ly0 = LEGEND_XY
    if lowest_y is not None and lowest_y - 0.25 < ly0:
        ly0 = lowest_y - 0.25
    _draw_legend(ax, has_truncated=has_truncated, has_onco=has_onco,
                 lx0=lx0, ly0=ly0)
    # re-fit canvas: legend extends to ly0 - 5*0.13 - 0.1; matplotlib
    # clips patches outside the axes, so the limits must cover it
    n_entries = 5 + (1 if has_onco else 0) + (1 if has_truncated else 0)
    needed_ymin = ly0 - (n_entries - 1) * 0.13 - 0.10
    if needed_ymin < -1.30:
        ax.set_ylim(-1.30, 1.55)
        ax.set_ylim(needed_ymin, 1.55)
        # keep aspect equal after the change
        ax.set_aspect("equal")

    return pd.DataFrame(gene_rows) if gene_rows else None


def _frag_angles(spans: List[Tuple], circ_len: int) -> List[Tuple[float, float]]:
    """Convert circular-coordinate spans to (t2, t1) ccw arc degrees.

    Mirrors the mapping in :func:`to_ang`: pos 0 at 12 o'clock,
    clockwise, ``90 - 360*pos/circ_len``.  Each fragment becomes one
    (t2, t1) arc with t2 <= t1 (matplotlib Wedge convention).

    Parameters
    ----------
    spans : list of tuple
        Circular coordinate spans (chrom, f_start, f_end, x_start,
        x_end, strand) of the cycle fragments.
    circ_len : int
        Total circular coordinate length.

    Returns
    -------
    list of (float, float)
        One (t2, t1) degree pair per fragment, BED order.
    """
    arcs = []
    for fc, fs, fe, xs, xe, _st in spans:
        a0 = 90.0 - 360.0 * xs / circ_len
        a1 = 90.0 - 360.0 * xe / circ_len
        arcs.append((a1, a0))
    return arcs


def _topo_wedge(ax: plt.Axes, t2: float, t1: float, color: str,
                bands: RingBands, zorder: int = 8) -> None:
    """Draw one topology-band arc segment (full band height)."""
    ax.add_patch(Wedge((0, 0), bands.topo_out, t2, t1,
                       width=bands.topo_out - bands.topo_in,
                       facecolor=color, edgecolor="none", zorder=zorder))


def _draw_topology_ring(ax: plt.Axes, topology: str,
                        spans: List[Tuple], circ_len: int,
                        bands: RingBands) -> None:
    """Draw the topology mini-ring at the cycle's real fragment angles.

    A thin band hugging the hub renders the topology motif using the
    legend's visual language (topology_legend.png), but positioned at
    the REAL fragment coordinates so junction gaps, colour joins and
    motif marks align radially with the coverage/gene rings outside:

    * simple_circle            — backbone with one white excision gap
    * simple_sv                — backbone + orange tick at the deletion
    * mixed_sv                 — same as simple_sv + green inversion arc
    * multi_region_intra_chr   — per-fragment shades of one hue
    * multi_region_inter_chr   — per-fragment chromosome colours
    * simple_duplications      — deep red first copy, pale red repeats
    * foldback                  — backbone + converging cones at the
      foldback junction

    Single-fragment cycles (the common simple_circle case) have no
    interior junction inside the ring — their excision point sits at
    the ring's coordinate origin (pos 0 / 12 o'clock), where the BED
    span wraps; the white gap is drawn there.  Unrecognised topology
    names leave the band empty (hub text still names the topology).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    topology : str
        topology_name from summary.txt.
    spans : list of tuple
        Circular coordinate spans of the cycle fragments.
    circ_len : int
        Total circular coordinate length.
    """
    arcs = _frag_angles(spans, circ_len)
    n = len(arcs)
    topo = (topology or "").strip().lower()
    r_hi, r_lo = bands.topo_out, bands.topo_in
    band_h = r_hi - r_lo

    def backbone() -> None:
        """Backbone ring with the white excision gap centred at 12 o'clock.

        The gap is centred exactly on the coordinate origin (math angle
        90 deg) with the shared ``TOPO_GAP_DEG`` width — the coverage
        ring cuts the identical wedge there, so the breakpoint reads
        with the same width in every ring.
        """
        _topo_wedge(ax, 90.0 + TOPO_GAP_DEG / 2.0,
                    90.0 + 360.0 - TOPO_GAP_DEG / 2.0,
                    TOPO_BASE, bands)

    def sv_ticks() -> None:
        """Orange full-band tick at each inter-fragment junction.

        The junction between consecutive fragments i, i+1 sits at the
        circular coordinate ``xs_{i+1}``, i.e. the **t1** end of arc
        i+1 (its low-coordinate end) — never at the ring origin.
        """
        for _t2, t1 in arcs[1:]:
            _topo_wedge(ax, t1 - TOPO_TICK_DEG / 2.0,
                        t1 + TOPO_TICK_DEG / 2.0, TOPO_SV, bands, zorder=9)

    # topology flag decode: duplications and foldbacks may combine with
    # any base layout; decoil names them directly, but the band renders
    # the dominant motif per topology_name
    if topo == "simple_circle":
        # one backbone ring with a single full-band white gap at the
        # excision breakpoint: the coordinate origin (pos 0 = 12
        # o'clock = math angle 90 deg), where the fragment span wraps
        backbone()
        return

    if n == 0:
        return

    if topo == "simple_sv":
        # backbone over all fragments; the SV (the small deleted piece)
        # is the gap BETWEEN two consecutive fragments: mark the flanking
        # junction with an orange full-band tick.  The backbone gap sits
        # at the coordinate origin (12 o'clock) like simple_circle.
        backbone()
        sv_ticks()
        return

    if topo == "mixed_sv":
        # like simple_sv, plus the inversion arm: fragments alternate
        # strand in a mixed cycle — the inverted run renders green with
        # a single arrow tip at its high-coordinate end
        backbone()
        sv_ticks()
        strands = [st for *_x, st in spans]
        run_start = None
        for i, st in enumerate(strands):
            if st == "-":
                if run_start is None:
                    run_start = i
            elif run_start is not None:
                t2 = arcs[run_start][0]
                t1 = arcs[i - 1][1]
                # green body across the inverted run
                _topo_wedge(ax, t2 + TOPO_EPS_DEG, t1, TOPO_INV, bands,
                            zorder=9)
                # single arrow tip at the run's t2 end (ccw arc tip)
                tip_len = min(band_h, 0.06)
                body_hi = r_lo + tip_len
                tip_ang = np.radians(t2)
                tri = [
                    (r_hi * np.cos(tip_ang), r_hi * np.sin(tip_ang)),
                    (r_hi * np.cos(tip_ang - 0.09),
                     r_hi * np.sin(tip_ang - 0.09)),
                    (r_lo * np.cos(tip_ang - 0.02),
                     r_lo * np.sin(tip_ang - 0.02)),
                    (body_hi * np.cos(tip_ang), body_hi * np.sin(tip_ang)),
                    (r_lo * np.cos(tip_ang + 0.02),
                     r_lo * np.sin(tip_ang + 0.02)),
                    (r_hi * np.cos(tip_ang + 0.09),
                     r_hi * np.sin(tip_ang + 0.09)),
                ]
                ax.add_patch(plt.Polygon(tri, closed=True,
                                         facecolor=TOPO_INV,
                                         edgecolor="none", zorder=10))
                run_start = None
        if run_start is not None:
            # inverted run reaches the last fragment; render body only
            t2, t1 = arcs[run_start]
            _topo_wedge(ax, t2 + TOPO_EPS_DEG, t1, TOPO_INV, bands, zorder=9)
        return

    if topo in ("multi_region_intra_chr", "multi_region_inter_chr"):
        # per-fragment coloured segments separated by white slits of the
        # shared TOPO_GAP_DEG width — every fragment junction (including
        # the wrap-around at the coordinate origin) is a real breakpoint
        # and reads identically in the coverage ring; intra = shades of
        # one hue, inter = one colour per source chromosome
        if topo == "multi_region_intra_chr":
            seq = [TOPO_BASE, "#8CA6C0", TOPO_BASE, "#8CA6C0"]
        else:
            # one colour per distinct chromosome (legend palette)
            chrom_order = list(dict.fromkeys(fc for fc, *_ in spans))
            cmap = {c: CHR_COLORS[(i * 2) % len(CHR_COLORS)]
                    for i, c in enumerate(chrom_order)}
            seq = [cmap[fc] for fc, *_ in spans]
        for i, (t2, t1) in enumerate(arcs):
            col = seq[i % len(seq)]
            # shrink both arc ends half a slit: t2 faces the next
            # junction, t1 the previous one — leaves one slit width of
            # white at every junction, angles wrap harmlessly
            half = min(TOPO_GAP_DEG / 2.0, (t1 - t2) / 4.0)
            _topo_wedge(ax, t2 + half, t1 - half, col, bands)
        return

    if topo == "simple_duplications":
        # same fragment twice on the ring: deep red = first occurrence
        # (A), pale red = repeat (A'); junction slits as everywhere else
        seen: Dict[Tuple, int] = {}
        counts = [(fc, int(fs), int(fe)) for fc, fs, fe, *_ in spans]
        for key in counts:
            seen[key] = seen.get(key, 0) + 1
        for (t2, t1), key in zip(arcs, counts):
            col = TOPO_DUP_HI if seen[key] == 1 else TOPO_DUP_LO
            half = min(TOPO_GAP_DEG / 2.0, (t1 - t2) / 4.0)
            _topo_wedge(ax, t2 + half, t1 - half, col, bands)
        return

    if topo == "foldback":
        # foldback join: two palindromic arms converge at the junction;
        # render the backbone plus converging cones at the junction
        # angle where the arms meet (adjacent-fragment boundary with a
        # foldback link — for single-arm foldbacks the junction sits at
        # the coordinate origin)
        backbone()
        juncs = [arcs[0][0]] if n == 1 else [arcs[1][0]]
        for ja in juncs:
            for sign in (-1, 1):
                a0 = np.radians(ja + sign * TOPO_TICK_DEG)
                a1 = np.radians(ja)
                r_tip = r_lo + 0.01
                tri = [
                    (r_hi * np.cos(a0), r_hi * np.sin(a0)),
                    (r_hi * np.cos(a1), r_hi * np.sin(a1)),
                    (r_tip * np.cos(a1), r_tip * np.sin(a1)),
                ]
                col = TOPO_FB_A if sign < 0 else TOPO_FB_B
                ax.add_patch(plt.Polygon(tri, closed=True,
                                         facecolor=col,
                                         edgecolor="none", zorder=9))
        return

    # unrecognised topology: leave the band empty; hub text names it


def _draw_hub(ax: plt.Axes, info: Optional[pd.Series], circ_id: str,
              spans: List[Tuple], circ_len: int,
              bands: RingBands) -> None:
    """Draw the central hub circle with cycle metadata.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    info : pandas.Series or None
        Summary row (topology) for this cycle; None = 'NA'.
    circ_id : str
        Cycle identifier for the hub title.
    spans : list of tuple
        Circular coordinate spans (first/last fragment coordinates).
    circ_len : int
        Total circular coordinate length.
    bands : RingBands
        Radial band allocation (hub disc radius).
    """
    ax.add_patch(plt.Circle((0, 0), bands.hub, facecolor="white",
                            edgecolor="#D5DBE1", linewidth=1.0, zorder=6))
    # text rows scale with the hub radius so the disc never clips them
    ax.text(0, 0.085, f"circ {circ_id}", ha="center", fontsize=10.5,
            color=INDIGO, fontweight="bold", zorder=7)
    topology = str(info["topology"]) if info is not None else "NA"
    ax.text(0, -0.02, topology, ha="center", fontsize=6.6,
            color=GRAY, zorder=7)
    ax.text(0, -0.125, f"{circ_len / 1000:.1f} kb", ha="center", fontsize=7.5,
            color=GRAY, zorder=7)
    # coordinate summary: one line per fragment ("+"-joined), each with
    # its chromosome — a multi-region cycle must show every piece, not
    # a first..last range that hides the interior breakpoints
    if spans:
        parts = [f"{fc}:{int(fs):,}-{int(fe):,}"
                 for fc, fs, fe, _xs, _xe, _st in spans]
        coord_text = " + ".join(parts)
        # cap: >3 fragments collapse to count (text would overflow the hub)
        if len(parts) > 3:
            coord_text = (f"{parts[0]} + {parts[1]} + … "
                          f"({len(parts)} segments)")
        # long text wraps to two rows so it never spills past the disc
        if len(coord_text) > 46:
            mid = coord_text.rfind(" + ", 0, len(coord_text) // 2 + 8)
            if mid > 0:
                coord_text = (coord_text[:mid] + " +\n"
                              + coord_text[mid + 3:])
        # deterministic font fit: figure sizes are fixed by the callers
        # (9.4x8 per-cycle, 9x9 grid), so data-units-per-inch is known;
        # digits at fontsize F occupy ~0.60*F points of width -> font
        # sized so the widest row spans <= 80% of the hub diameter.
        # Row count scales the row pitch; the block is vertically
        # centred in the lower half of the disc so 1-3 rows all fit.
        fig_w_in = ax.figure.get_size_inches()[0]
        ax_w = ax.get_xlim()[1] - ax.get_xlim()[0]
        data_per_pt = ax_w / (fig_w_in * 72.0)
        widest = max(len(row) for row in coord_text.split("\n"))
        avail_w = 1.6 * bands.hub
        fit_size = avail_w / (widest * 0.60 * data_per_pt)
        n_rows = 1 + coord_text.count("\n")
        row_pitch = max(min(fit_size, 5.0), 2.2) * 1.6 * data_per_pt
        y_top = -0.16 - (n_rows - 1) * row_pitch / 2.0
        ax.text(0, y_top, coord_text, ha="center", va="center",
                fontsize=max(min(fit_size, 5.0), 2.2),
                linespacing=1.5, color=GRAY, zorder=7)



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
    # slash strokes: compact marks inside the block height — oversized
    # slashes visually drowned thin intron-backbone genes
    slash_len = exon_w * 0.45
    for i in range(n_gaps):
        off = (i - (n_gaps - 1) / 2.0) * 0.005
        # each slash: 45 deg between tangent and radial
        d = (tang + norm) / np.sqrt(2.0)
        p0 = np.array([cx, cy]) + tang * off - d * slash_len / 2.0
        p1 = np.array([cx, cy]) + tang * off + d * slash_len / 2.0
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
                color=col, linewidth=0.9, solid_capstyle="round",
                zorder=7)


def _draw_legend(ax: plt.Axes, has_truncated: bool = True,
                 has_onco: bool = False,
                 lx0: Optional[float] = None,
                 ly0: Optional[float] = None) -> None:
    """Draw a frameless legend of functional elements outside the map.

    Horizontal straight swatches, one per element, grouped by role:

    * data layer   — coverage fill (gold)
    * gene elements— exon (solid bar), intron (thin bar), UTR (pale
      bar); all in the gene colour, shape/opacity carries the meaning
    * state marks  — oncogene outline, ``//`` junction break (only
      shown when present in this cycle)

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    has_truncated : bool
        Whether truncated genes occur in this cycle.
    has_onco : bool
        Whether oncogene-highlighted blocks occur in this cycle.
    """
    # legend anchored in the empty bottom-right quadrant (LEGEND_XY),
    # overridden by the caller when gene labels extend further
    if lx0 is None or ly0 is None:
        lx0, ly0 = LEGEND_XY
    bar_w = 0.07            # legend swatch width
    bar_h = 0.03            # legend swatch height (visual size, cosmetic)
    gene_col = "#1A1A1A"    # gene-structure swatches in black; per-gene
                            # colours on the ring carry identity, the
                            # legend only explains shape
    entries: List[Tuple[str, str, str]] = [
        ("patch", GOLD, "coverage (log scale)"),
        ("exon", gene_col, "exon"),
        ("intron", gene_col, "intron"),
        ("utr", gene_col, "UTR"),
        ("micro_exon", gene_col, "small gene (dot)"),
    ]
    if has_onco:
        entries.append(("onco_outline", gene_col, "oncogene"))
    if has_truncated:
        entries.append(("break", gene_col, "junction break (//)"))

    for i, (kind, col, txt) in enumerate(entries):
        yy = ly0 - i * 0.13
        if kind == "patch":
            ax.add_patch(plt.Rectangle(
                (lx0, yy - bar_h / 2.0), bar_w, bar_h,
                facecolor=col, edgecolor="none", zorder=8))
        elif kind == "exon":
            ax.add_patch(plt.Rectangle(
                (lx0, yy - bar_h / 2.0), bar_w, bar_h,
                facecolor=col, edgecolor="white", linewidth=0.25,
                zorder=8))
        elif kind == "intron":
            # dashed line, same style as the on-circle intron arcs
            ax.plot([lx0, lx0 + bar_w], [yy, yy],
                    color=col, linewidth=1.0,
                    linestyle=(0, (6, 4)),
                    solid_capstyle="butt", zorder=8)
        elif kind == "micro_exon":
            # dot glyph: tiny exons render as dots on the ring
            ax.plot([lx0 + bar_w / 2.0], [yy],
                    marker="o", markersize=3.5, color=col,
                    markeredgecolor="none", zorder=8)
        elif kind == "utr":
            # pale fill as on the circle: gene colour under white
            ax.add_patch(plt.Rectangle(
                (lx0, yy - bar_h / 2.0), bar_w, bar_h,
                facecolor=col, edgecolor="none", zorder=8))
            ax.add_patch(plt.Rectangle(
                (lx0, yy - bar_h / 2.0), bar_w, bar_h,
                facecolor="white", edgecolor="none", alpha=0.65,
                zorder=9))
        elif kind == "onco_outline":
            ax.add_patch(plt.Rectangle(
                (lx0, yy - bar_h / 2.0), bar_w, bar_h,
                facecolor=col, edgecolor="#3A3A70", linewidth=0.9,
                zorder=8))
        elif kind == "strand":
            # taper glyph matching the on-circle convention: exon block
            # narrowing to a point along the reading direction
            ax.add_patch(plt.Polygon(
                [(lx0 + bar_w, yy),
                 (lx0 + bar_w * 0.45, yy + bar_h / 2.0),
                 (lx0 + bar_w * 0.45, yy - bar_h / 2.0),
                 (lx0, yy - bar_h / 2.0),
                 (lx0, yy + bar_h / 2.0)],
                closed=True, facecolor=col, edgecolor="none", zorder=8))
        elif kind == "break":
            # just the two diagonal slashes, no background swatch
            for k in (-1, 1):
                x_mid = lx0 + bar_w / 2.0 + k * 0.005
                dx, dy = 0.008, bar_h * 0.55
                ax.plot([x_mid - dx, x_mid + dx],
                        [yy - dy, yy + dy],
                        color=col, linewidth=0.9,
                        solid_capstyle="round", zorder=9)
        ax.text(
            lx0 - 0.10, yy, txt,
            fontsize=5.8, color="#3B4A5A",
            va="center", ha="right", zorder=8)


# Canonical columns of the per-gene functional table; used for EVERY write
# (including empty ones) so the declared <out-prefix>_genes.tsv output always
# exists with a stable header — snakemake deletes declared outputs of failed
# jobs, and a missing table previously failed the whole rule on samples with
# no detected ecDNA.
GENE_TABLE_COLUMNS = ["circ_id", "gene", "status", "onco", "strand", "region",
                      "median_coverage", "tpm", "label_x", "label_y"]


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
    exons_df, utr_df : pandas.DataFrame or None
        Exon / UTR intervals (see :func:`load_gtf_features`).

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
        bw_path, tpm, bw_mean_cov, exons_df, utr_df)


def _grid_page(
        page_ids: List[str], page: int, n_pages: int, n_cycles: int,
        suffix: str, bed: pd.DataFrame, summary: pd.DataFrame,
        genes: pd.DataFrame, onco_genes: set, bw_path: Optional[str],
        tpm: Dict[str, float], bw_mean_cov: float,
        sample_id: str, out_prefix: str, image_formats: List[str],
        exons_df: Optional[pd.DataFrame],
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
    exons_df, utr_df : pandas.DataFrame or None
        Exon / UTR intervals (see :func:`load_gtf_features`).

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
            bw_path, tpm, bw_mean_cov, exons_df, utr_df)

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
            reset_gene_hues()   # per-cycle unique colour assignment
            # canvas ratio must match the data ratio (xlim span 3.35 vs
            # ylim span 2.85) so aspect='equal' doesn't letterbox the
            # ring into a corner with vast empty margins
            fig = plt.figure(figsize=(9.4, 8.0))
            rows = _draw_and_collect(
                fig, cid, bed, summary, genes, onco_genes,
                bw_path, tpm, bw_mean_cov, exons_df, utr_df)
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
                exons_df, utr_df))
        else:
            for page in range(n_pages):
                all_gene_rows.extend(_grid_page(
                    circ_ids[page * max_per_fig:
                             (page + 1) * max_per_fig],
                    page, n_pages, len(circ_ids), f".p{page + 1}",
                    bed, summary, genes, onco_genes, bw_path, tpm,
                    bw_mean_cov, sample_id, out_prefix, image_formats,
                    exons_df, utr_df))
    else:
        n_pages = (len(circ_ids) + max_per_fig - 1) // max_per_fig
        if n_pages == 1:
            all_gene_rows.extend(_grid_page(
                circ_ids, 0, 1, len(circ_ids), "", bed, summary,
                genes, onco_genes, bw_path, tpm, bw_mean_cov,
                sample_id, out_prefix, image_formats,
                exons_df, utr_df))
        else:
            for page in range(n_pages):
                all_gene_rows.extend(_grid_page(
                    circ_ids[page * max_per_fig:
                             (page + 1) * max_per_fig],
                    page, n_pages, len(circ_ids), f".p{page + 1}",
                    bed, summary, genes, onco_genes, bw_path, tpm,
                    bw_mean_cov, sample_id, out_prefix, image_formats,
                    exons_df, utr_df))

    # Per-gene functional table (all cycles) — always written (empty table
    # with stable header when no gene rows) so the declared output exists.
    table = (pd.concat(all_gene_rows, ignore_index=True) if all_gene_rows
             else pd.DataFrame(columns=GENE_TABLE_COLUMNS))
    table = table.reindex(columns=GENE_TABLE_COLUMNS)
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
    utrs = None
    if args.gtf:
        if os.path.isfile(args.gtf):
            feats = load_gtf_all_features(args.gtf)
            genes = feats["genes"]
            exons = feats["exons"]
            utrs = feats["utrs"]
            logger.info(f"Loaded {len(genes)} genes, {len(exons)} exons, "
                        f"{len(utrs)} UTR from {args.gtf}")
        else:
            logger.warning(f"GTF not found ({args.gtf}); "
                           "drawing without gene annotation")

    if bed.empty:
        # Detection guard: nothing passed the filter -> no structures to draw
        logger.info(f"{args.sample}: no ecDNA detected, drawing placeholder")
        plot_placeholder(args.out_prefix, image_formats, args.sample)
        # Keep the output contract: write an empty genes table too, or the
        # missing <out-prefix>_genes.tsv fails the snakemake rule (which then
        # also deletes the freshly drawn placeholder figure).
        pd.DataFrame(columns=GENE_TABLE_COLUMNS).to_csv(
            f"{args.out_prefix}_genes.tsv", sep="\t", index=False)
        return

    logger.info(f"{args.sample}: {bed['circ_id'].nunique()} ecDNA cycles, "
                f"{len(bed)} fragments")
    onco_genes = load_onco_genes(args.onco_genes)
    tpm = load_tpm(args.tpm)
    plot_structures(bed, summary, genes, onco_genes, args.bw, tpm,
                    args.sample, args.out_prefix,
                    image_formats, args.max_per_fig,
                    one_file_per_cycle=args.one_file_per_cycle,
                    exons_df=exons, utr_df=utrs)


if __name__ == "__main__":
    main()
