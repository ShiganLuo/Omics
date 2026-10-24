#!/usr/bin/env python3
"""Breakpoint-sharing family figure for ecDNA circles.

One row per family: sample-colored circle tracks, shared/private breakpoint
lines, and (optionally) GENCODE gene models below.

布局(与 2026-10-23 15:35 稳定版一致):所有家族行等高,标题在每行上方
内部预留带中,基因模型/标签在轨道下方的固定带内,底部为全局元素图例。
"""
from __future__ import annotations

import argparse
import colorsys
import csv
import itertools
import logging

logger = logging.getLogger(__name__)
import os
import re
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ---------------------------------------------------------------------------
# 数据加载
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse and validate command-line arguments for the plotting workflow.

    Args:
        argv: Optional command-line tokens. When ``None``, argparse reads
            ``sys.argv``; supplying a list is useful for tests and wrappers.

    Returns:
        An ``argparse.Namespace`` containing the input directory, sample IDs,
        breakpoint-matching settings, optional GTF path, output prefix, and
        requested image formats.

    Raises:
        SystemExit: If required options are missing or an option is invalid.
    """
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("-d","--decoil-dir", required=True)
    parser.add_argument("-s", "--samples", action="append", required=True)
    parser.add_argument("-t", "--tolerance", type=int, default=200)
    parser.add_argument("-m", "--min-shared", type=int, default=2)
    parser.add_argument("-g", "--gtf", default=None,
                        help="GENCODE GTF annotation file (enables gene models)")
    parser.add_argument("-o", "--out-prefix", default="breakpoint_families")
    parser.add_argument("-f", "--format", action="append", dest="formats",
                        help="图片格式,可重复。默认 png。")
    return parser.parse_args(argv)


def load_circles(decoil_dir: str, samples: List[str]) -> Dict[Tuple[str, str], List[Tuple[str, int, int]]]:
    """Read Decoil's filtered fragment BED files for the requested samples.

    For each sample, this function reads
    ``<decoil_dir>/<sample>/decoil/reconstruct.ecDNA.filtered.bed``. Blank and
    comment lines are ignored; rows with fewer than four columns are skipped.
    The first four columns are interpreted as chromosome, BED start, BED end,
    and circle ID. Coordinates remain 0-based, half-open as in BED.

    Args:
        decoil_dir: Root directory containing one subdirectory per sample.
        samples: Sample IDs whose Decoil BED files should be loaded.

    Returns:
        A dictionary keyed by ``(sample_id, circle_id)``. Each value is a list
        of ``(chromosome, start, end)`` fragment tuples. A circle can have
        multiple fragments, retained in input-file order.

    Raises:
        FileNotFoundError: If any requested sample BED file is missing.
        ValueError: If a usable row contains non-integer start/end coordinates.
    """
    circles = {}
    for s in samples:
        path = os.path.join(decoil_dir, s, "decoil",
                            "reconstruct.ecDNA.filtered.bed")
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                circles.setdefault((s, f[3]), []).append(
                    (f[0], int(f[1]), int(f[2])))
    return circles


def load_genes(gtf_path: Optional[str]) -> List[tuple]:
    """Parse gene structures from a GTF annotation into plotting records.

    Only ``gene``, ``exon``, and UTR feature rows are used. GTF's 1-based,
    closed coordinates are converted to 0-based, half-open intervals. Features
    are grouped by ``gene_id``; exon and UTR intervals are deduplicated and
    sorted. If multiple gene IDs share a gene symbol, the ID is appended to
    the symbol so that each plotted record has a unique key. Gzipped GTF files
    are supported when the path ends in ``.gz``.

    Args:
        gtf_path: GTF filename, optionally gzip-compressed. ``None`` or an
            empty string disables gene annotation and returns an empty list.

    Returns:
        A list of tuples with fields ``(name, start, end, tag, exons, utrs,
        strand, chromosome)``. ``exons`` and ``utrs`` are lists of coordinate
        pairs; ``tag`` initially equals ``"full"`` and can later be changed to
        ``"partial"`` when the gene extends outside the visible fragments.

    Raises:
        FileNotFoundError: If a non-empty GTF path does not exist.
    """
    if not gtf_path:
        return []
    if not os.path.isfile(gtf_path):
        raise FileNotFoundError(f"GTF file does not exist: {gtf_path}")

    attribute_pattern = re.compile(r'([^;\s]+)\s+"([^"]*)"')
    records: Dict[str, Dict[str, Any]] = {}
    opener = __import__("gzip").open if gtf_path.endswith(".gz") else open

    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                logger.warning("Skipping malformed GTF line: %s", line.rstrip())
                continue

            chrom, _, feature, start, end, _, strand, _, attributes = fields
            if feature not in {"gene", "exon", "UTR", "five_prime_UTR", "three_prime_UTR"}:
                continue
            attrs = dict(attribute_pattern.findall(attributes))
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue

            # GTF is 1-based inclusive; convert to BED-compatible coordinates.
            start0, end0 = int(start) - 1, int(end)
            item = records.setdefault(gene_id, {
                "name": attrs.get("gene_name", gene_id),
                "chrom": chrom,
                "start": start0,
                "end": end0,
                "strand": strand,
                "exons": [],
                "utrs": [],
            })
            if item["chrom"] != chrom:
                logger.warning("Ignoring inconsistent chromosome for gene %s", gene_id)
                continue
            item["start"] = min(item["start"], start0)
            item["end"] = max(item["end"], end0)

            if feature == "exon":
                item["exons"].append((start0, end0))
            elif feature in {"UTR", "five_prime_UTR", "three_prime_UTR"}:
                item["utrs"].append((start0, end0))

    # Keep names unique because plotting lane assignments use the name as a key.
    name_counts: Dict[str, int] = defaultdict(int)
    for item in records.values():
        name_counts[item["name"]] += 1

    genes = []
    for gene_id, item in records.items():
        name = item["name"]
        if name_counts[name] > 1:
            name = f"{name}|{gene_id}"
        exons = sorted(set(item["exons"]))
        utrs = sorted(set(item["utrs"]))
        genes.append((name, item["start"], item["end"], "full",
                      exons, utrs, item["strand"], item["chrom"]))

    logger.info("Loaded %d genes from GTF: %s", len(genes), gtf_path)
    return genes


# ---------------------------------------------------------------------------
# 家族聚类
# ---------------------------------------------------------------------------

def edge_set(frags: List[Tuple[str, int, int]]) -> Set[Tuple[str, int]]:
    """Extract the unique endpoints of all fragment intervals.

    Args:
        frags: Fragment tuples ``(chromosome, start, end)`` using BED
            coordinates. Each fragment contributes both its start and end.

    Returns:
        A set of ``(chromosome, coordinate)`` pairs. Repeated endpoints are
        collapsed, so the result represents distinct breakpoint coordinates,
        not the number of fragments supporting them.
    """
    out = set()
    for ch, s, e in frags:
        out.add((ch, s))
        out.add((ch, e))
    return out


def shared_breakpoints(e1: Set[Tuple[str, int]], e2: Set[Tuple[str, int]], tol: int) -> List[Tuple[str, int, int]]:
    """Find endpoint pairs that match between two circles within a tolerance.

    An endpoint pair is considered a match only when both endpoints are on
    the same chromosome and their absolute coordinate difference is no larger
    than ``tol``. This function returns pairwise matches; it does not cluster
    nearby coordinates or deduplicate multiple matches to one endpoint.

    Args:
        e1: Unique endpoints for the first circle.
        e2: Unique endpoints for the second circle.
        tol: Maximum allowed absolute coordinate difference, in base pairs.

    Returns:
        A list of ``(chromosome, coordinate_from_e1, coordinate_from_e2)``
        tuples, one for each qualifying endpoint pair.
    """
    return [(c1, x, y) for c1, x in e1 for c2, y in e2
            if c1 == c2 and abs(x - y) <= tol]


def group_families(circles: Dict[Tuple[str, str], List[Tuple[str, int, int]]], tol: int, min_shared: int) -> List[List[Tuple[str, str]]]:
    """Cluster circles into breakpoint-sharing families across samples.

    Every circle starts in its own union-find set. Circle pairs from the same
    sample are deliberately not compared. For each cross-sample pair, the
    circles are linked when the number of pairwise endpoint matches returned
    by :func:`shared_breakpoints` reaches ``min_shared``. Union-find makes
    this relation transitive: a family can contain circles whose endpoints do
    not directly match each other, provided a chain of links connects them.

    Args:
        circles: Mapping from ``(sample_id, circle_id)`` to fragment tuples.
        tol: Maximum endpoint-coordinate difference for a match, in base pairs.
        min_shared: Minimum number of pairwise endpoint matches required to
            link two circles. This is a count of matches, not necessarily a
            count of distinct breakpoint clusters.

    Returns:
        Families represented as lists of ``(sample_id, circle_id)`` keys,
        sorted by decreasing family size and then by member key for stability.
    """
    uf = {k: k for k in circles}

    def find(x: Tuple[str, str]) -> Tuple[str, str]:
        """Return the union-find representative while applying path compression."""
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x

    for (s1, c1), (s2, c2) in itertools.combinations(circles, 2):
        if s1 == s2:
            continue
        hits = shared_breakpoints(edge_set(circles[(s1, c1)]),
                                  edge_set(circles[(s2, c2)]), tol)
        if len(hits) >= min_shared:
            ra, rb = find((s1, c1)), find((s2, c2))
            if ra != rb:
                uf[ra] = rb
    fams = defaultdict(list)
    for k in circles:
        fams[find(k)].append(k)
    return sorted(fams.values(), key=lambda m: (-len(m), min(m)))


# ---------------------------------------------------------------------------
# 布局计算(纯函数)
# ---------------------------------------------------------------------------

def adjust_color_saturation(color: Any, saturation: float) -> str:
    """Scale a color's HSV saturation while preserving its hue and value.

    Args:
        color: Any Matplotlib-compatible color specification.
        saturation: Saturation multiplier in the inclusive range [0, 1]. A
            value of 1 preserves the original color; 0 produces a gray tone
            with the same value (brightness). Intermediate values soften
            highly saturated colors without blending them toward white.

    Returns:
        The adjusted color as a hexadecimal RGB string.

    Raises:
        ValueError: If ``saturation`` is outside [0, 1].
    """
    if not 0.0 <= saturation <= 1.0:
        raise ValueError("color_saturation must be between 0 and 1")
    red, green, blue = matplotlib.colors.to_rgb(color)
    hue, original_saturation, value = colorsys.rgb_to_hsv(red, green, blue)
    adjusted = colorsys.hsv_to_rgb(
        hue, original_saturation * saturation, value
    )
    return matplotlib.colors.to_hex(adjusted)


def gene_hues(fam_genes: list, color_saturation: float = 0.55) -> Dict[str, str]:
    """Assign a deterministic colormap color to each distinct gene name.

    The discrete ``turbo`` colormap is sampled at the number of unique names,
    so there is no hard-coded maximum number of genes. The saturation of each
    sampled color is scaled by ``color_saturation`` to avoid overly vivid
    annotations. Very large gene sets can still contain perceptually similar
    colors.

    Args:
        fam_genes: Gene annotation tuples associated with one family; the first
            tuple field is used as the unique display name.
        color_saturation: Saturation multiplier in [0, 1]. Lower values make
            gene colors more muted; 1 retains the original colormap saturation.

    Returns:
        A dictionary mapping each gene name to a hexadecimal color string.
        An empty input produces an empty dictionary.
    """
    names = list(dict.fromkeys(g[0] for g in fam_genes))
    if not names:
        return {}
    cmap = plt.get_cmap("tab20", len(names))
    return {
        name: adjust_color_saturation(cmap(i), color_saturation)
        for i, name in enumerate(names)
    }


def plan_lanes(fam_genes: list) -> Tuple[Dict[str, int], int]:
    """Assign gene models to vertical lanes so overlapping spans do not collide.

    Genes are processed from longest to shortest. Each gene is placed in the
    first lane whose previously assigned intervals do not overlap it; a new
    lane is created if no existing lane is available. Lane indices start at 0.

    Args:
        fam_genes: Gene tuples whose first three fields are name, start, and
            end. The intervals are expected to be clipped to the plotted span.

    Returns:
        A pair ``(lane_by_gene, lane_count)``. ``lane_by_gene`` maps each gene
        name to its integer lane index; ``lane_count`` is the number of lanes.
    """
    order = sorted(fam_genes, key=lambda g: -(g[2] - g[1]))
    lane_iv: List[list] = []
    lane_of: dict = {}
    for g in order:
        placed = False
        for li, ivs in enumerate(lane_iv):
            if all(g[2] < s or g[1] > e for s, e in ivs):
                ivs.append((g[1], g[2]))
                lane_of[g[0]] = li
                placed = True
                break
        if not placed:
            lane_of[g[0]] = len(lane_iv)
            lane_iv.append([(g[1], g[2])])
    return lane_of, len(lane_iv)


def plan_labels(fam_genes: list, lo: int, hi: int, pad: float, panel_width_inches: float, label_font_size: float, char_width_factor: float, color_saturation: float) -> Tuple[list, int]:
    """Place gene labels and wrap horizontally overlapping labels into rows.

    Approximate label widths are estimated from the displayed genomic span,
    panel width, font size, and a configurable character-width factor. Genes
    are considered longest-first. Each label is centered on the visible gene
    interval, clamped to the horizontal plotting limits, and assigned to the
    first row where its estimated span does not overlap another label.
    Partial genes receive a trailing ``*`` in the displayed label.

    Args:
        fam_genes: Gene tuples with name, clipped start/end, and tag fields.
        lo: Left edge of the family's genomic fragment span.
        hi: Right edge of the family's genomic fragment span.
        pad: Extra genomic padding used on both sides of the x-axis.
        panel_width_inches: Approximate physical width available for this panel.
        label_font_size: Label font size in points.
        char_width_factor: Multiplier approximating average character width
            relative to the font size.

    Returns:
        ``(labels, row_count)``. Each label is ``(row_index, x_left, text,
        color)``; ``row_count`` is the number of rows required. Empty input
        returns an empty label list and zero rows.
    """
    x0_lim, x1_lim = lo - pad * 2, hi + pad * 2
    data_per_pt = (x1_lim - x0_lim) / (panel_width_inches * 72.0)
    char_w = label_font_size * char_width_factor * data_per_pt
    hues = gene_hues(fam_genes, color_saturation)
    rows: List[List[Tuple[float, float]]] = []
    out = []
    for name, gs, ge, tag, *_ in sorted(fam_genes,
                                        key=lambda g: -(g[2] - g[1])):
        txt = f"{name}{'*' if tag == 'partial' else ''}"
        txt_w = len(txt) * char_w
        vis_s, vis_e = max(gs, lo), min(ge, hi)
        center = (vis_s + vis_e) / 2
        x_l = center - txt_w / 2
        x_l = max(x0_lim, min(x_l, x1_lim - txt_w))
        x_r = x_l + txt_w
        placed = False
        for ri, spans in enumerate(rows):
            if all(x_r < s or x_l > e for s, e in spans):
                spans.append((x_l, x_r))
                out.append((ri, x_l, txt, hues[name]))
                placed = True
                break
        if not placed:
            rows.append([(x_l, x_r)])
            out.append((len(rows) - 1, x_l, txt, hues[name]))
    return out, len(rows)


def plan_breakpoints(members: List[Tuple[str, str]], circles: Dict[Tuple[str, str], List[Tuple[str, int, int]]]) -> Tuple[list, list]:
    """Separate approximate shared and private endpoints within a family.

    Endpoint coordinates are traversed in sorted order. For each unprocessed
    endpoint, members with an endpoint on the same chromosome within 200 bp
    are counted as carriers; nearby endpoints are marked as used to avoid
    repeatedly reporting the same local cluster. The tolerance is currently
    fixed at 200 bp in this helper and is independent of the family-grouping
    tolerance passed to :func:`group_families`.

    Args:
        members: Circle keys belonging to the family.
        circles: Mapping from circle keys to their fragment intervals.

    Returns:
        ``(shared, private)``. Shared records contain chromosome, a
        representative coordinate, carrier count, and carrier circle keys.
        Private records contain chromosome, coordinate, and carrier keys.
    """
    edge_sets = {m: edge_set(circles[m]) for m in members}
    all_edges = set().union(*edge_sets.values())
    shared, private = [], []
    used = set()
    for ch, pos in sorted(all_edges):
        if (ch, pos) in used:
            continue
        carriers = {}
        for m, es in edge_sets.items():
            for (c2, y) in es:
                if c2 == ch and abs(pos - y) <= 200:
                    carriers[m] = y
        used.update((c2, y) for (c2, y) in all_edges
                    if c2 == ch and abs(pos - y) <= 200)
        if len(carriers) >= 2:
            shared.append((ch, min(carriers.values()), len(carriers),
                           sorted(carriers)))
        else:
            private.append((ch, pos, sorted(carriers)))
    return shared, private


def plan_family(members: List[Tuple[str, str]], circles: Dict[Tuple[str, str], List[Tuple[str, int, int]]], fam_genes: list, head: float, lane_pitch: float, label_row_height: float, track_to_genes: float, axis_padding: float, band_base: float, panel_width_inches: float, label_font_size: float, char_width_factor: float, color_saturation: float) -> Dict[str, Any]:
    """Calculate genomic bounds, gene-label layout, and vertical space for a family.

    The family's x-span is the minimum fragment start through the maximum
    fragment end across all member circles. Genes are assigned lanes and
    labels are wrapped into rows when annotations are present. The calculated
    gene band reserves enough vertical space for gene tracks, labels, and the
    gap below the circle tracks; the returned y-limits include both the gene
    band and the header area.

    Args:
        members: Circle keys belonging to this family.
        circles: Circle-to-fragment mapping used to derive family bounds.
        fam_genes: Gene annotation tuples already clipped to visible fragment
            intervals, including any metadata needed for drawing.
        head: Vertical space reserved above the tracks for the family header.
        lane_pitch: Vertical distance between adjacent gene lanes.
        label_row_height: Vertical distance between adjacent label rows.
        track_to_genes: Gap between circle tracks and the first gene lane.
        axis_padding: Reserved vertical space around the x-axis/labels.
        band_base: Additional fixed height added to the gene band.
        panel_width_inches: Physical width used to estimate label collisions.
        label_font_size: Gene-label font size in points.
        char_width_factor: Character-width estimate used by label placement.
        color_saturation: Saturation multiplier used for gene-label colors.

    Returns:
        A dictionary containing members, chromosome and coordinate bounds,
        padding, track count, gene records, lane assignments, labels, gene-band
        height, and y-axis limits. Drawing functions consume this dictionary.
    """
    all_frags = [f for m in members for f in circles[m]]
    lo = min(f[1] for f in all_frags)
    hi = max(f[2] for f in all_frags)
    pad = (hi - lo) * 0.03
    if fam_genes:
        lane_of, n_lanes = plan_lanes(fam_genes)
        labels, n_label_rows = plan_labels(fam_genes, lo, hi, pad, panel_width_inches, label_font_size, char_width_factor, color_saturation)
    else:
        lane_of, n_lanes, labels, n_label_rows = {}, 0, [], 0
    gene_band = (track_to_genes + n_lanes * lane_pitch
                 + n_label_rows * label_row_height
                 + axis_padding + band_base)
    return dict(members=members, chrom=all_frags[0][0], lo=lo, hi=hi,
                pad=pad, n_tracks=len(members), fam_genes=fam_genes,
                lane_of=lane_of, n_lanes=n_lanes, labels=labels,
                n_label_rows=n_label_rows, gene_band=gene_band,
                ylim=(-(gene_band + 0.2), len(members) + head))


# ---------------------------------------------------------------------------
# 绘制(每个视觉层一个函数)
# ---------------------------------------------------------------------------

def draw_tracks(ax: Any, layout: Dict[str, Any], circles: Dict[Tuple[str, str], List[Tuple[str, int, int]]], sample_colors: Dict[str, str]) -> None:
    """Draw one horizontal fragment track for each circle in a family.

    Each fragment is rendered as a rectangle at the y-position assigned to
    its circle. All fragments from a sample share that sample's color, while
    a thin dark border keeps adjacent fragment intervals distinguishable. A
    right-aligned sample/circle label is placed to the left of the genomic
    span.

    Args:
        ax: Matplotlib axes receiving the rectangles and labels.
        layout: Family layout dictionary returned by :func:`plan_family`.
        circles: Mapping from circle keys to chromosome/start/end fragments.
        sample_colors: Mapping from sample IDs to hexadecimal color strings.

    Returns:
        None. The supplied axes are modified in place.
    """
    lo, pad = layout["lo"], layout["pad"]
    n_tr = layout["n_tracks"]
    for ti, m in enumerate(layout["members"]):
        y0 = n_tr - 1 - ti
        s, c = m
        color = sample_colors[s]
        for ch, fs, fe in circles[m]:
            ax.add_patch(Rectangle((fs, y0 + 0.18), fe - fs, 0.64,
                                   facecolor=color, edgecolor="#333333",
                                   lw=0.5, zorder=3))
        ax.text(lo - pad * 1.5, y0 + 0.5, f"{s}:circ{c}", fontsize=5.0,
                ha="right", va="center", color=color)



def draw_breakpoints(
    ax: Any,
    layout: Dict[str, Any],
    shared: list,
    private: list,
) -> None:
    """Draw shared and private breakpoint markers without text labels.

    Args:
        ax: Matplotlib axes receiving breakpoint markers.
        layout: Family layout containing ordered circle members.
        shared: Shared breakpoint records.
        private: Private breakpoint records.

    Returns:
        None. The supplied axes are modified in place.
    """
    members = layout["members"]
    n_tr = layout["n_tracks"]

    for ch, pos, n, carriers in shared:
        idx = [n_tr - 1 - members.index(m) for m in carriers]
        y_bot = min(idx) + 0.12
        y_top = max(idx) + 0.88

        ax.plot(
            [pos, pos],
            [y_bot, y_top],
            color="#222222",
            lw=1.4,
            ls="--",
            zorder=4,
        )

    for ch, pos, carriers in private:
        m = carriers[0]
        ti = n_tr - 1 - members.index(m)

        ax.plot(
            [pos, pos],
            [ti + 0.02, ti + 0.98],
            color="#aaaaaa",
            lw=0.7,
            ls="--",
            zorder=4,
        )

def _gene_entry(layout: Dict[str, Any], name: str) -> Optional[tuple]:
    """Look up a family's annotation tuple by its unique display name.

    Args:
        layout: Family layout dictionary containing ``fam_genes``.
        name: Unique gene display name to locate.

    Returns:
        The matching gene tuple, or ``None`` when the name is not present.
    """
    for g in layout["fam_genes"]:
        if g[0] == name:
            return g
    return None


def draw_gene_models(ax: Any, layout: Dict[str, Any], track_to_genes: float, bar_height: float, lane_pitch: float, color_saturation: float) -> None:
    """Render clipped gene structures below the circle-fragment tracks.

    Exons are drawn as solid gene-colored blocks, intronic intervals as thin
    connecting lines, and UTR intervals as narrower gray overlays. Exon and
    UTR intervals are intersected with each visible fragment block before
    drawing, so no annotation is drawn through uncovered genomic gaps. Intron
    connectors are only drawn between exons belonging to the same continuous
    visible block. A strand-direction mark is shown only when the annotated
    3-prime end remains inside a visible block.

    Args:
        ax: Matplotlib axes receiving the gene models.
        layout: Family layout containing clipped gene records, visible blocks,
            and the lane assignment for each gene.
        track_to_genes: Vertical gap between the circle tracks and gene lanes.
        bar_height: Height of each exon block in y-axis data units.
        lane_pitch: Vertical distance between adjacent gene lanes.
        color_saturation: Saturation multiplier used for gene-model colors.

    Returns:
        None. The supplied axes are modified in place. If the family has no
        gene annotations, the function exits without drawing anything.
    """
    fam_genes = layout["fam_genes"]
    if not fam_genes:
        return
    lo, hi = layout["lo"], layout["hi"]
    lane_of = layout["lane_of"]
    hues = gene_hues(fam_genes, color_saturation)
    min_w = (hi - lo) * 0.0025
    for name, gs, ge, tag, *_ in fam_genes:
        lane = lane_of[name]
        mid = -(track_to_genes + 0.5 * bar_height) - lane * lane_pitch
        col = hues[name]
        entry = _gene_entry(layout, name)
        exons = entry[4] if entry and len(entry) > 4 else []
        utrs = entry[5] if entry and len(entry) > 5 else []
        strand = entry[6] if entry and len(entry) > 6 else "+"
        visible_blocks = entry[8] if entry and len(entry) > 8 else [(lo, hi)]
        ex_c = []
        for s, e in exons:
            for block_start, block_end in visible_blocks:
                s2, e2 = max(s, block_start), min(e, block_end)
                if e2 <= s2:
                    continue
                if e2 - s2 < min_w:
                    mp = (s2 + e2) / 2
                    s2, e2 = mp - min_w / 2, mp + min_w / 2
                ex_c.append((s2, e2, block_start, block_end))
        ex_c.sort()
        # Connect exons only when they lie in the same uninterrupted fragment block.
        for left, right in zip(ex_c, ex_c[1:]):
            if left[2:] == right[2:]:
                ax.plot([left[1], right[0]], [mid, mid],
                        color=col, lw=0.6, zorder=2)
        # Do not fabricate an exon when the GTF has no exon annotation.
        for s2, e2, _, _ in ex_c:
            ax.add_patch(Rectangle((s2, mid - bar_height / 2), e2 - s2, bar_height,
                                   facecolor=col, edgecolor="white",
                                   lw=0.3, zorder=3))
        for s, e in utrs:
            for block_start, block_end in visible_blocks:
                s2, e2 = max(s, block_start), min(e, block_end)
                if e2 <= s2:
                    continue
                ax.add_patch(Rectangle((s2, mid - bar_height / 4), e2 - s2, bar_height / 2,
                                       facecolor="#cfd8dc", edgecolor="none",
                                       zorder=4))
        # Draw the strand arrow only if the true 3' end remains inside a fragment.
        true_three_prime = (entry[10] if strand == "+" else entry[9]) if entry and len(entry) > 10 else (ge if strand == "+" else gs)
        if any(block_start <= true_three_prime <= block_end
               for block_start, block_end in visible_blocks):
            ax_size = (hi - lo) * 0.004
            if strand == "+":
                x3, x_tip = true_three_prime - ax_size, true_three_prime
            else:
                x3, x_tip = true_three_prime + ax_size, true_three_prime
            for xs_, ys_, ye_ in ([x3, x_tip], mid, mid), \
                    ([x_tip, x3], mid, mid + bar_height * 0.16), \
                    ([x_tip, x3], mid, mid - bar_height * 0.16):
                ax.plot(xs_, [ys_, ye_], color=col, lw=0.7, zorder=5,
                        path_effects=[pe.Stroke(linewidth=1.2,
                                                foreground="white"),
                                      pe.Normal()])


def draw_labels(ax: Any, layout: Dict[str, Any], track_to_genes: float, lane_pitch: float, label_row_height: float, label_font_size: float) -> None:
    """Draw gene labels at the positions calculated by :func:`plan_labels`.

    Labels are placed below the gene-model lanes. Their text, x-coordinate,
    row, and color come from the precomputed layout; clipping is enabled to
    prevent text from extending outside the axes region.

    Args:
        ax: Matplotlib axes receiving the text labels.
        layout: Family layout dictionary containing label tuples and lane count.
        track_to_genes: Vertical gap between circle tracks and gene lanes.
        lane_pitch: Vertical distance between adjacent gene lanes.
        label_row_height: Vertical distance between adjacent label rows.
        label_font_size: Gene-label font size in points.

    Returns:
        None. The supplied axes are modified in place. If there are no labels,
        the function exits without drawing anything.
    """
    labels = layout["labels"]
    if not labels:
        return
    n_lanes = layout["n_lanes"]
    base = -(track_to_genes + n_lanes * lane_pitch) - 0.12
    for ri, x_l, txt, colr in labels:
        ly = base - ri * label_row_height
        ax.text(x_l, ly, txt, fontsize=label_font_size, ha="left", va="top",
                color=colr, clip_on=True)


def draw_header(ax: Any, layout: Dict[str, Any], fi: int) -> None:
    """Add the family title and its genomic/sample summary above the tracks.

    The header reports the family index, chromosome interval, total number of
    circles, and number of distinct samples represented by those circles.

    Args:
        ax: Matplotlib axes receiving the header text.
        layout: Family layout dictionary containing member IDs and genomic span.
        fi: One-based family index used in the displayed title.

    Returns:
        None. The supplied axes are modified in place.
    """
    members = layout["members"]
    lo, pad = layout["lo"], layout["pad"]
    n_tr = layout["n_tracks"]
    ax.text(lo - pad * 8, n_tr + 0.18,
            f"family {fi}   {layout['chrom']}:{layout['lo']:,}"
            f"-{layout['hi']:,}   "
            f"({len(members)} circles, "
            f"{len({s for s, _ in members})} samples)",
            fontsize=7.5, va="bottom", ha="left")


def style_axes(ax: Any, x_min: float, x_max: float) -> None:
    """Apply panel axis styling and constrain the bottom spine to the fragment span.

    Y-axis ticks and the top/left/right spines are hidden. The bottom spine is
    retained but its visible bounds are set to the family's minimum and maximum
    fragment coordinates; x tick labels remain visible without tick marks.

    Args:
        ax: Matplotlib axes to format.
        x_min: Minimum fragment coordinate for this family.
        x_max: Maximum fragment coordinate for this family.

    Returns:
        None. The supplied axes are modified in place.
    """
    ax.set_yticks([])
    for spine in ("left", "top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_bounds(x_min, x_max)
    ax.tick_params(axis="x", labelsize=5.2, length=0, pad=3,
                   labelbottom=True, bottom=False)
    ax.ticklabel_format(axis="x", style="plain")


# ---------------------------------------------------------------------------
# 总控:计算布局 → 分配空间 → 调用各绘制函数
# ---------------------------------------------------------------------------

def plot(
    families: List[List[Tuple[str, str]]],
    circles: Dict[Tuple[str, str], List[Tuple[str, int, int]]],
    samples: List[str],
    out_prefix: str,
    formats: Optional[List[str]],
    genes: Optional[List[tuple]] = None,
    head: float = 0.62,
    lane_pitch: float = 0.55,
    bar_height: float = 0.34,
    label_row_height: float = 0.30,
    track_to_genes: float = 0.16,
    axis_padding: float = 0.30,
    band_base: float = 0.20,
    label_font_size: float = 4.2,
    char_width_factor: float = 0.72,
    panel_width_inches: float = 9.0 * 0.68,
    figure_width_inches: float = 9.0,
    scale_inches_per_unit: float = 0.25,
    row_gap_inches: float = 0.24,
    color_saturation: float = 0.9,
) -> None:
    """Compute layouts, draw every family panel, and save the combined figure.

    The plotting pipeline first intersects each gene with the union of the
    family's fragment intervals, retaining separate visible blocks so gene
    models cannot bridge uncovered gaps. It then calculates panel layouts,
    assigns dynamic colors to samples, draws tracks/breakpoints/annotations,
    adds a shared legend, and writes one image per requested format.

    Args:
        families: Breakpoint-sharing families; each family is a list of
            ``(sample_id, circle_id)`` keys, usually from :func:`group_families`.
        circles: Mapping from each circle key to its fragment intervals.
        samples: Sample IDs used to assign consistent colors across panels.
        out_prefix: Output path prefix; the format extension is appended.
        formats: Image formats such as ``["png", "pdf"]``. ``None`` or an
            empty list defaults to PNG.
        genes: Optional records from :func:`load_genes`. If absent or empty,
            gene lanes and labels are omitted.
        head: Vertical headroom above each family's circle tracks.
        lane_pitch: Vertical spacing between adjacent gene-model lanes.
        bar_height: Height of exon rectangles in y-axis data units.
        label_row_height: Vertical spacing between gene-label rows.
        track_to_genes: Gap from circle tracks to the first gene lane.
        axis_padding: Extra reserved vertical space around the axis and labels.
        band_base: Additional fixed height reserved for the gene band.
        label_font_size: Gene-label font size in points.
        char_width_factor: Approximate character width as a fraction of font
            size; used to estimate label collisions before rendering.
        panel_width_inches: Effective width used to estimate label widths.
        figure_width_inches: Width of the complete output figure in inches.
        scale_inches_per_unit: Converts each panel's calculated vertical span
            into physical figure height.
        row_gap_inches: Additional physical separation between family panels.
        color_saturation: Saturation multiplier in [0, 1] applied consistently
            to both sample-track and gene colors. Use a lower value (for
            example, 0.35) for muted colors; 1 preserves the colormap's
            original saturation and 0 converts the colors to grayscale.

    Returns:
        None. Files named ``<out_prefix>.<format>`` are written to disk.

    Raises:
        OSError: If an input file cannot be read or an output file cannot be
            written. Matplotlib may also raise for unsupported output formats.
    """
    if not 0.0 <= color_saturation <= 1.0:
        raise ValueError("color_saturation must be between 0 and 1")

    # pass 1:所有家族的布局一次性算完(纯计算)
    layouts = []
    for fi, members in enumerate(families, 1):
        sorted_members = sorted(members)
        family_fragments = [fragment for member in sorted_members
                            for fragment in circles[member]]
        family_chrom = family_fragments[0][0]
        family_lo = min(fragment[1] for fragment in family_fragments)
        family_hi = max(fragment[2] for fragment in family_fragments)
        # Merge the family's actual fragment intervals. Gene models may only be
        # drawn inside these blocks, not across internal breakpoint gaps.
        fragment_intervals = sorted(
            (start, end) for chrom, start, end in family_fragments
            if chrom == family_chrom
        )
        coverage = []
        for start, end in fragment_intervals:
            if coverage and start <= coverage[-1][1]:
                coverage[-1] = (coverage[-1][0], max(coverage[-1][1], end))
            else:
                coverage.append((start, end))

        fam_genes = []
        for gene in (genes or []):
            name, gene_start, gene_end, tag, exons, utrs, strand, chrom = gene
            if chrom != family_chrom or gene_end <= family_lo or gene_start >= family_hi:
                continue

            visible_blocks = [
                (max(gene_start, start), min(gene_end, end))
                for start, end in coverage
                if min(gene_end, end) > max(gene_start, start)
            ]
            if not visible_blocks:
                continue
            clipped_start = min(start for start, _ in visible_blocks)
            clipped_end = max(end for _, end in visible_blocks)
            covered_bp = sum(end - start for start, end in visible_blocks)
            is_partial = covered_bp < gene_end - gene_start

            clipped_exons = [
                (max(exon_start, block_start), min(exon_end, block_end))
                for exon_start, exon_end in exons
                for block_start, block_end in visible_blocks
                if min(exon_end, block_end) > max(exon_start, block_start)
            ]
            clipped_utrs = [
                (max(utr_start, block_start), min(utr_end, block_end))
                for utr_start, utr_end in utrs
                for block_start, block_end in visible_blocks
                if min(utr_end, block_end) > max(utr_start, block_start)
            ]
            fam_genes.append((name, clipped_start, clipped_end,
                              "partial" if is_partial else tag,
                              clipped_exons, clipped_utrs, strand, chrom,
                              visible_blocks, gene_start, gene_end))

        layouts.append(plan_family(
            sorted_members, circles, fam_genes, head, lane_pitch,
            label_row_height, track_to_genes, axis_padding, band_base,
            panel_width_inches, label_font_size, char_width_factor,
            color_saturation
        ))


    n_fam = len(families)
    head_span = head + 0.30
    for m in layouts:
        span = m["n_tracks"] + 0.2 + m["gene_band"] + head_span
        m["row_weight"] = span
    height_ratios = [m["row_weight"] for m in layouts]
    # Reserve physical space between panels for the lower x ticks and the next header.
    fig_h = (sum(height_ratios) * scale_inches_per_unit
             + n_fam * row_gap_inches)
    fig, axes = plt.subplots(n_fam, 1, figsize=(figure_width_inches, fig_h),
                             squeeze=False,
                             gridspec_kw={"height_ratios": height_ratios,
                                          "hspace": 0.30,
                                          "bottom": 0.06, "top": 0.99})

    # Assign a distinct colormap entry to every sample, without a fixed sample list.
    sample_names = list(dict.fromkeys(samples))
    sample_cmap = plt.get_cmap("Set2", max(1, len(sample_names)))
    sample_colors = {
        name: adjust_color_saturation(sample_cmap(i), color_saturation)
        for i, name in enumerate(sample_names)
    }

    # pass 2:绘制
    for fi, (ax, layout) in enumerate(zip(axes[:, 0], layouts), 1):
        ax.set_ylim(*layout["ylim"])
        ax.set_xlim(layout["lo"] - layout["pad"] * 2,
                    layout["hi"] + layout["pad"] * 2)
        shared, private = plan_breakpoints(layout["members"], circles)
        draw_header(ax, layout, fi)
        draw_tracks(ax, layout, circles, sample_colors)
        draw_breakpoints(ax, layout, shared, private)
        draw_gene_models(ax, layout, track_to_genes, bar_height, lane_pitch, color_saturation)
        draw_labels(ax, layout, track_to_genes, lane_pitch, label_row_height, label_font_size)
        style_axes(ax, layout["lo"], layout["hi"])

    elem_handles = [
        plt.Line2D([0], [0], color="#222222", lw=1.4, ls="--",
                   label="shared breakpoint (n circles)"),
        plt.Line2D([0], [0], color="#aaaaaa", lw=0.7, ls="--",
                   label="private breakpoint (single circle)"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#bbbbbb", edgecolor="#333333",
                      lw=0.5, label="ecDNA fragment"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#4a708b", edgecolor="white",
                      lw=0.4, label="exon (one color per gene; * = partial)"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#cfd8dc", edgecolor="none",
                      label="UTR"),
        plt.Line2D([0], [0], color="#4a708b", lw=0.6,
                   label="intron (chevron = strand 3' end)"),
    ]
    fig.legend(elem_handles, [h.get_label() for h in elem_handles],
               fontsize=5.8, frameon=False, loc="lower left",
               bbox_to_anchor=(0.08, 0.01), ncol=3)
    for fmt in (formats or ["png"]):
        path = f"{out_prefix}.{fmt}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        logger.info("saved %s", path)
    plt.close(fig)


def main(argv: Optional[List[str]] = None) -> int:
    """Execute the command-line pipeline from input loading to figure output.

    The workflow loads Decoil fragments, groups circles into families using
    breakpoint similarity, optionally parses a GTF, and passes the results to
    :func:`plot`. Logging configuration is established here rather than at
    module import time, so importing this module does not reconfigure the
    caller's root logger.

    Args:
        argv: Optional command-line tokens; ``None`` uses ``sys.argv``.

    Returns:
        Integer process status: 0 when the workflow completes successfully.

    Raises:
        SystemExit: Propagated from argument parsing for invalid CLI input.
        OSError: If required input files cannot be read or outputs cannot be
            written.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )
    args = parse_args(argv)
    circles = load_circles(args.decoil_dir, args.samples)
    logger.info("Loaded %d circles", len(circles))
    families = group_families(circles, args.tolerance, args.min_shared)
    logger.info("%d families", len(families))
    genes = load_genes(args.gtf)
    if genes:
        logger.info("Loaded %d gene models from GTF", len(genes))
    plot(families, circles, args.samples, args.out_prefix, args.formats,
         genes=genes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
