#!/usr/bin/env python3
"""Paper-style fiber-seq locus browser — Peter et al., Cell Rep Methods 2024,
100911 Fig 2H/1D scheme, per the paper's figure legends, methods, and pixels.

What the paper draws (verified against figure legends + methods + the figure
itself at 500 dpi):
  - Fiber rows are UCSC BED12-style: a BLACK connecting line carries the
    molecule; its visible segments between blocks are nucleosomes ('n').
    Purple rectangles ('L') mark linker DNA.  MSPs scoring as regulatory
    elements ('P' promoter / 'RE') are colored by FIRE FDR tier:
    orange 0.10>=q>0.05, bright red 0.05>=q>0.01, dark red-brown q<=0.01.
  - (brown) m5CpG track per cell type, stacked at the top of the Fiber-seq
    section: dark red-brown filled histogram.
  - (black) FIRE score track per sample: S_g = -(50/C_g)*sum(log10(q_i)) over
    FIRE elements overlapping g (FIRE pipeline README); unscored (<3 FIREs)
    as in the pipeline.  Drawn as black bars ("black blocks").
  - Fiber coverage track: multi-color stacked element coverages — nucleosome
    coverage (grey, top), linker coverage (purple, middle), FIRE coverage
    (red, bottom), i.e. the FIRE pipeline's element coverages.
  - Each fiber stack is bracketed at the left with a green label box
    ("sample Single Fibers").
  - No coordinate axes on subplots (UCSC screenshot style); a kb scale bar
    sits under the last fiber stack.  Real data only: no decorative banner
    and no hand-drawn schematic key.

Per-call FIRE q: the Ma/AQ tag spec (molecular-annotation crate) writes
`fire.Q:` sections whose per-call linear qualities live in the `Aq` B-array,
Aq = precision*255 (u8); q = 1 - Aq/255.  The FIRE pipeline caps precision at
0.99 which rounds to 252 in u8, so Aq>=252 is the q<=0.01 tier.
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# palette sampled from the paper figure (500 dpi pixels)
NUC_LINE = "#1A1A1A"           # black connecting line = nucleosome 'n'
LINKER = "#8A78B6"             # purple rectangle = linker DNA 'L'
FIRE_ORANGE = "#F0A246"        # MSP q>0.05 tier (paper "orange"/yellow)
FIRE_RED = "#E01B24"           # MSP 0.01<q<=0.05 tier (bright red)
FIRE_DKRED = "#7D3219"         # MSP q<=0.01 tier (dark red-brown)
MCPG_FILL = "#7D3219"          # brown m5CpG histogram
FIRE_BAR = "#1A1A1A"           # black FIRE score bars
COV_NUC = "#B0B1B5"            # grey   = nucleosome coverage (stack top)
COV_LNK = "#8A78B6"            # purple = linker coverage (stack middle)
COV_FIRE = "#C0392B"           # red    = FIRE coverage (stack bottom)
BANNER_BG = "#D2E8D2"          # pale green "Fiber-seq" banner
GENE_EXON = "#2F4F6F"          # dark slate = gene-structure exon boxes
LABEL_BOX = "#D8ECD8"          # pale green sample label box
HILITE = "#C9E4F2"             # light blue highlight of the zoomed fiber
HILITE_EC = "#7FB3D0"
GAP_ANN = "#7A7A7A"
VALID_FORMATS = ("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg",
                 "jpeg", "pgf", "raw", "rgba")

logger = logging.getLogger("fiberbrowser")


def setup_logging() -> None:
    """Initialize module logging (spec: bin scripts use logging, not print)."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")

# ---------------- SAM parsing ----------------

def read_to_ref(pos0: int, cigar: str) -> Dict[int, int]:
    ref_of_read: Dict[int, int] = {}
    ref, read = pos0, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            for i in range(n):
                ref_of_read[read + i] = ref + i
            read += n
            ref += n
        elif op in "DN":
            ref += n
        else:
            read += n
    return ref_of_read

def aligned_ref_segments(pos0: int, cigar: str) -> List[Tuple[int, int]]:
    """Merged (ref_start, ref_end) of aligned blocks (M/=/X); gaps = D/N."""
    segs: List[Tuple[int, int]] = []
    ref = pos0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            if segs and segs[-1][1] == ref:
                segs[-1] = (segs[-1][0], segs[-1][1] + n)
            else:
                segs.append((ref, ref + n))
            ref += n
        elif op in "DN":
            ref += n
    return segs

def mm_positions(entry: str, seq: str) -> List[int]:
    """MM tag entry -> read positions of modified bases."""
    m = re.match(r"([ACGTN])([+-])([a-z]+)(\.?)", entry)
    if not m:
        return []
    base = m.group(1)
    nums = entry[m.end():].lstrip(",?")
    try:
        deltas = [int(x) for x in nums.split(",") if x != ""]
    except ValueError:
        return []
    idx = [i for i, c in enumerate(seq) if c.upper() == base]
    out, cur = [], -1
    for d in deltas:
        cur = d if not out else cur + 1 + d
        if 0 <= cur < len(idx):
            out.append(idx[cur])
        else:
            break
    return out

def clip_blocks_to_segs(blocks: List[Tuple[int, int]],
                        segs: List[Tuple[int, int]]
                        ) -> List[Tuple[int, int]]:
    """Split ref blocks at aligned-segment boundaries so nothing bridges a
    deletion gap (a Ma block spanning a D in read coords would otherwise be
    stretched across the deleted reference stretch)."""
    out: List[Tuple[int, int]] = []
    for a, b in blocks:
        for sa, sb in segs:
            lo, hi = max(a, sa), min(b, sb)
            if hi > lo:
                out.append((lo, hi))
    return out

def clip_blocks_q(blocks: List[Tuple[int, int, int]],
                  segs: List[Tuple[int, int]]
                  ) -> List[Tuple[int, int, int]]:
    out: List[Tuple[int, int, int]] = []
    for a, b, q in blocks:
        for sa, sb in segs:
            lo, hi = max(a, sa), min(b, sb)
            if hi > lo:
                out.append((lo, hi, q))
    return out

def parse_blocks_q(ma_value: str, key: str, rmap: Dict[int, int],
                   aq: List[int], aq_off: List[int]
                   ) -> List[Tuple[int, int, int]]:
    """Ma `key` blocks -> (ref_start, ref_end, quality).  `aq_off` is a
    one-element counter indexing into the Aq array (qualities are grouped
    per annotation in Ma-string order; in these files only `fire` has
    qualities, and the Aq length matches the fire block count exactly)."""
    out: List[Tuple[int, int, int]] = []
    last_read = max(rmap) if rmap else 0
    for seg in ma_value.split(";"):
        m = re.match(rf"{re.escape(key)}\.([A-Za-z0-9]*):", seg)
        if not m:
            continue
        for blk in seg[m.end():].split(","):
            if "-" not in blk:
                continue
            s, l = blk.split("-")[:2]
            rs, rl = int(s), int(l)
            q = aq[aq_off[0]] if aq_off[0] < len(aq) else 255
            aq_off[0] += 1
            a = rmap.get(rs)
            b = rmap.get(min(rs + rl - 1, last_read))
            if a is None or b is None:
                continue
            out.append((min(a, b), max(a, b), q))
    return out

def parse_blocks(ma_value: str, key: str, rmap: Dict[int, int]
                 ) -> List[Tuple[int, int]]:
    """Ma `key` blocks (no quality) -> reference-coordinate (start, end).

    `ma_value` is the Ma tag body; `key` names the section (nuc/msp);
    `rmap` maps read to reference positions.  Blocks whose endpoints do
    not map (outside the alignment) are dropped.
    """
    blocks: List[Tuple[int, int]] = []
    last_read = max(rmap) if rmap else 0
    for seg in ma_value.split(";"):
        m = re.match(rf"{re.escape(key)}\.([A-Za-z0-9]*):", seg)
        if not m:
            continue
        for blk in seg[m.end():].split(","):
            if "-" not in blk:
                continue
            s, l = blk.split("-")[:2]
            rs, rl = int(s), int(l)
            a = rmap.get(rs)
            b = rmap.get(min(rs + rl - 1, last_read))
            if a is None or b is None:
                continue
            blocks.append((min(a, b), max(a, b)))
    return blocks

def fire_tier_color(aq: int) -> str:
    """Aq = precision*255 (u8); q = 1 - Aq/255.  Paper tiers by FDR q."""
    if aq >= 252:      # precision cap 0.99 -> u8 252 (q<=0.01)
        return FIRE_DKRED
    if aq >= 243:      # 0.01 < q <= 0.05
        return FIRE_RED
    return FIRE_ORANGE  # 0.05 < q <= 0.10

def _sam_records(path: str):
    """Yield (name, flag, ref_name, pos0, cigar, seq, tags) from a text SAM.

    Tag values keep their SAM strings; Aq arrives comma-joined so both
    record sources share one downstream parse.
    """
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            tags: Dict[str, str] = {}
            for t in f[11:]:
                p = t.split(":", 2)
                if len(p) == 3:
                    tags[p[0]] = p[2]
            yield (f[0], int(f[1]), f[2], int(f[3]) - 1, f[5], f[9], tags)


_CIGAR_OPS = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"}


def _bam_records(path: str, chrom: str, start: int, end: int):
    """Yield the same 7-tuples from BAM/CRAM via pysam.

    Fetches only the locus when an index is present, else scans the file.
    Aq arrays are re-serialized to the text-SAM form so one parser serves
    both sources.
    """
    import pysam                       # lazy: text-SAM runs need no pysam
    bam = pysam.AlignmentFile(path, "rb")
    try:
        try:
            it = bam.fetch(chrom, start, end)
        except ValueError:             # missing index -> scan everything
            it = bam.fetch(until_eof=True)
        for rec in it:
            if rec.is_unmapped or rec.reference_name != chrom:
                continue
            tags: Dict[str, str] = {}
            for key in ("Ma", "MM"):
                if rec.has_tag(key):
                    tags[key] = str(rec.get_tag(key))
            if rec.has_tag("Aq"):
                tags["Aq"] = ",".join(str(int(v)) for v in rec.get_tag("Aq"))
            cigar = "".join(f"{n}{_CIGAR_OPS[op]}" for op, n in
                            (rec.cigartuples or []))
            yield (rec.query_name or "", rec.flag, chrom,
                   rec.reference_start, cigar,
                   rec.query_sequence or "", tags)
    finally:
        bam.close()


def _iter_records(path: str, chrom: str, start: int, end: int):
    """Dispatch by suffix: .bam/.cram via pysam, anything else text SAM."""
    if path.lower().endswith((".bam", ".cram")):
        yield from _bam_records(path, chrom, start, end)
    else:
        yield from _sam_records(path)


def parse_sam(path: str, chrom: str, start: int, end: int) -> List[dict]:
    """Primary fibers overlapping the locus -> molecule records.

    Accepts text SAM or indexed BAM/CRAM.  Secondary/supplementary
    alignments (flag 0x900) are skipped so split segments never masquerade
    as whole molecules.
    """
    fibers: List[dict] = []
    for name, flag, ref_name, pos0, cigar, seq, tags in _iter_records(
            path, chrom, start, end):
        if ref_name != chrom or flag & 0x900:
            continue
        if not cigar or cigar == "*":
            continue
        segs = aligned_ref_segments(pos0, cigar)
        if not segs or segs[-1][1] <= start or segs[0][0] >= end:
            continue
        rmap = read_to_ref(pos0, cigar)
        ma = tags.get("Ma", "")
        aq = [int(x) for x in tags.get("Aq", "").replace("B", "").split(",")
              if x.strip().isdigit()]
        cpg: List[int] = []
        for part in tags.get("MM", "").split(";"):
            # cytosine modifications: C+m/C+h on read Cs, G-m/G-h on read
            # Gs (minus strand CpG); the m5CpG track pools these calls
            if part[:2] in ("C+", "G-"):
                cpg += [rmap[p] for p in mm_positions(part, seq)
                        if p in rmap and start <= rmap[p] < end]
        nuc = clip_blocks_to_segs(parse_blocks(ma, "nuc", rmap), segs)
        msp = clip_blocks_to_segs(parse_blocks(ma, "msp", rmap), segs)
        fire_q = clip_blocks_q(parse_blocks_q(ma, "fire", rmap, aq, [0]),
                               segs)
        fibers.append({"name": name, "segs": segs, "nuc": nuc, "msp": msp,
                       "fire_q": fire_q, "cpg": cpg})
    return fibers

# ---------------- aggregation ----------------

def split_msp_by_fire(msp: List[Tuple[int, int]],
                      fire: List[Tuple[int, int, int]]
                      ) -> Tuple[List[Tuple[int, int]],
                                 List[Tuple[int, int, int]]]:
    """MSP intervals -> (linker purple parts, FIRE parts with q).

    FIRE blocks sit inside MSP blocks: the FIRE sub-intervals are regulatory
    elements (RE/P); the remaining MSP sub-intervals are linker DNA (L).
    """
    linkers: List[Tuple[int, int]] = []
    fires: List[Tuple[int, int, int]] = []
    for a, b in msp:
        cuts = sorted((max(a, fa), min(b, fb), q) for fa, fb, q in fire
                      if fa < b and a < fb)
        cur = a
        for ca, cb, q in cuts:
            if ca > cur:
                linkers.append((cur, ca))
            fires.append((ca, cb, q))
            cur = max(cur, cb)
        if cur < b:
            linkers.append((cur, b))
    return linkers, fires

# alignment gaps (CIGAR D/N) up to this many bp are bridged when drawing the
# fiber line: the molecule is continuous DNA, and at zoom scale a ~40 bp
# alignment gap renders as a visible break in the element chain.  Larger
# gaps (real structural deletions, e.g. the validated 377-1868 bp SVs) stay
# open.
BRIDGE_MAX_GAP = 150

def bridge_segs(segs: List[Tuple[int, int]],
                max_gap: int = BRIDGE_MAX_GAP) -> List[Tuple[int, int]]:
    """Merge alignment blocks separated by small deletions/skips (D/N).

    Used for drawing the connecting line and the 'n' regions so the fiber
    reads as one continuous chain; the raw segs remain untouched for
    coverage accounting.
    """
    out: List[Tuple[int, int]] = []
    for a, b in segs:
        if out and a - out[-1][1] <= max_gap:
            out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out

# moving-average window for the between-group difference tracks used by
# --zoom-center-mode
DIFF_SMOOTH_BP = 200

def diff_center(samples, group_labels, start: int, end: int, kind: str) -> int:
    """Genomic position of the maximal absolute between-group difference.

    kind: "diff-fire" (FIRE score tracks) or "diff-cpg" (m5CpG density).
    Exactly two sample groups are required.  Positions without aligned
    coverage in either group are ignored (a one-sided coverage gap must not
    drive the choice); the per-base difference is smoothed over
    DIFF_SMOOTH_BP before taking the maximum.
    """
    uniq = sorted(set(group_labels))
    if len(uniq) != 2:
        raise ValueError(
            f"{kind} zoom-centre needs exactly two sample groups, "
            f"got {uniq}")

    def pooled(g: str) -> List[dict]:
        return [f for (_lbl, fl), gr in zip(samples, group_labels)
                if gr == g for f in fl]

    fa, fb = pooled(uniq[0]), pooled(uniq[1])
    if kind == "diff-fire":
        sa = fire_score_track(fa, start, end)
        sb = fire_score_track(fb, start, end)
    else:
        sa = per_base(fa, "cpg", start, end)
        sb = per_base(fb, "cpg", start, end)
    ca = per_base(fa, "cov", start, end)
    cb = per_base(fb, "cov", start, end)
    sm = smooth([abs(a - b) for a, b in zip(sa, sb)], DIFF_SMOOTH_BP)
    best_i, best_v = None, -1.0
    for i in range(end - start):
        if ca[i] <= 0 or cb[i] <= 0:
            continue
        if sm[i] > best_v:
            best_i, best_v = i, sm[i]
    if best_i is None:
        logger.warning(
            f"{kind}: no position has coverage in both groups; "
            "falling back to the locus midpoint")
        return (start + end) // 2
    return start + best_i

def per_base(fibers: List[dict], kind: str, start: int, end: int
             ) -> List[float]:
    """Per-base signal of `kind` over [start, end).

    kind: "cpg" (m5CpG mark counts), "cov" (aligned depth, M/=/X only),
    "nuc" / "lnk" / "fire" (nucleosome / linker / FIRE-element coverage).
    Returns a float list of length end-start.
    """
    arr = [0.0] * (end - start)
    for rec in fibers:
        if kind == "cpg":
            for x in rec["cpg"]:
                if start <= x < end:
                    arr[x - start] += 1.0
        elif kind == "cov":            # true aligned depth (M/=/X blocks)
            for a, b in rec["segs"]:
                for i in range(max(a, start) - start, min(b, end) - start):
                    arr[i] += 1.0
        elif kind == "nuc":
            for a, b in rec["nuc"]:
                for i in range(max(a, start) - start, min(b, end) - start):
                    arr[i] += 1.0
        elif kind == "lnk":
            for a, b in split_msp_by_fire(rec["msp"], rec["fire_q"])[0]:
                for i in range(max(a, start) - start, min(b, end) - start):
                    arr[i] += 1.0
        elif kind == "fire":
            for a, b, _q in rec["fire_q"]:
                for i in range(max(a, start) - start, min(b, end) - start):
                    arr[i] += 1.0
    return arr

def smooth(arr: List[float], k: int = 20) -> List[float]:
    """Centered moving average of window k, edge-clamped at both ends."""
    n = len(arr)
    pre = [0.0]
    for v in arr:
        pre.append(pre[-1] + v)
    return [ (pre[min(i + k // 2 + 1, n)] - pre[max(0, i - k // 2)])
             / (min(i + k // 2 + 1, n) - max(0, i - k // 2))
             for i in range(n) ]

def fire_score_track(fibers: List[dict], start: int, end: int
                     ) -> List[float]:
    """FIRE score S_g = -(50/C_g) * sum(log10(q_i)), q_i = 1 - Aq/255.
    Positions with <3 overlapping FIRE elements are unscored (-1 in the FIRE
    pipeline) and plotted as 0."""
    cov = per_base(fibers, "cov", start, end)
    qs: List[List[float]] = [[] for _ in range(end - start)]
    for rec in fibers:
        for a, b, aq in rec["fire_q"]:
            q = 1.0 - aq / 255.0
            for i in range(max(a, start) - start, min(b, end) - start):
                qs[i].append(q)
    out = [0.0] * (end - start)
    for i, lst in enumerate(qs):
        if len(lst) >= 3 and cov[i] > 0:
            out[i] = -(50.0 / cov[i]) * sum(math.log10(q) for q in lst)
    return out

# ---------------- gene structure (optional GTF) ----------------

def parse_gtf_genes(gtf_path: str, chrom: str, start: int, end: int
                    ) -> List[dict]:
    """Gene models overlapping [start, end) from a GTF/GTF.gz.

    One pass over the annotation; keeps only `gene`/`exon` features of
    `chrom` inside the locus.  Exons are merged per gene_id (union of all
    transcripts) so the drawn model stays compact.  Returns models sorted by
    start: {'name', 'strand', 'start', 'end', 'exons': [(a, b), ...]}.
    """
    import gzip
    opener = gzip.open if gtf_path.endswith(".gz") else open
    genes: Dict[str, dict] = {}
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != chrom or f[2] not in ("gene", "exon"):
                continue
            a, b = int(f[3]), int(f[4])
            if b <= start or a >= end:
                continue
            attrs = f[8]
            m = re.search(r'gene_id "([^"]+)"', attrs)
            gid = m.group(1) if m else f"{chrom}:{a}-{b}"
            nm = re.search(r'gene_name "([^"]+)"', attrs)
            g = genes.setdefault(gid, {"name": nm.group(1) if nm and
                                       nm.group(1) else gid,
                                       "strand": f[6],
                                       "start": a, "end": b, "exons": []})
            g["start"] = min(g["start"], a)
            g["end"] = max(g["end"], b)
            if f[2] == "exon":
                g["exons"].append((a, b))
    for g in genes.values():
        merged: List[Tuple[int, int]] = []
        for a, b in sorted(set(g["exons"])):
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        g["exons"] = merged
    return sorted(genes.values(), key=lambda g: (g["start"], g["end"]))

def draw_gene_track(ax, genes: List[dict], start: int, end: int, span: int
                    ) -> None:
    """Gene-structure row(s): exon boxes on a strand-arrowed backbone line,
    one stacked row per gene model, name above each model."""
    ax.set_xlim(0, span)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    n_rows = max(1, len(genes))
    ax.set_ylim(-0.5, n_rows - 0.5)
    for idx, g in enumerate(genes):
        y = n_rows - 1 - idx                # first gene on top
        gs, ge = max(g["start"], start), min(g["end"], end)
        ax.plot([gs - start, ge - start], [y, y], lw=0.9, color=NUC_LINE,
                solid_capstyle="butt", zorder=1)
        for a, b in g["exons"]:
            lo, hi = max(a, start), min(b, end)
            if hi > lo:
                ax.add_patch(Rectangle((lo - start, y - 0.17), hi - lo, 0.34,
                                       facecolor=GENE_EXON, ec="none",
                                       zorder=2))
        step = max(span // 15, 200)         # strand chevrons along the line
        marker = ">" if g["strand"] != "-" else "<"
        x = gs + step // 2
        while x < ge - step // 4:
            ax.plot([x - start], [y], marker=marker, ms=2.6, color=NUC_LINE,
                    zorder=3)
            x += step
        ax.text((gs + ge) / 2 - start, y + 0.22, g["name"], ha="center",
                va="bottom", fontsize=6.5)
    ax.set_ylabel("gene", fontsize=7, rotation=0, ha="right", va="center")

# ---------------- plot ----------------

def plot_locus(locus, samples, out_base, image_formats=None,
               max_fibers=12, dpi=200, zoom=False, zoom_bp=2100, gtf=None,
               zoom_center=None, zoom_center_mode=None, group_labels=None):
    """Render the paper-style locus figure (Fig 2H scheme) to disk.

    locus: (chrom, start, end, name); samples: [(label, fibers)] in display
    order; out_base: output prefix (writes out_base + "." + fmt per format);
    image_formats: list of matplotlib formats, None -> ["png"]; max_fibers:
    fiber rows drawn per sample (most aligned bases first); zoom: add the
    paper-style selection + zoom row; zoom_bp: width of the zoomed region.
    """
    chrom, start, end, name = locus
    formats = image_formats or ["png"]     # spec: None -> ["png"]
    span = end - start
    n = len(samples)
    genes = parse_gtf_genes(gtf, chrom, start, end) if gtf else []
    n_gene_rows = max(1, len(genes))
    gene_ratio = 0.5 + 0.35 * n_gene_rows
    # layout top to bottom: m5CpG tracks of all samples (stacked at the top),
    # then per-sample (FIRE score, coverage, single fibers), then the optional
    # gene-structure track (--gtf), then the optional zoom row (--zoom).
    if zoom:
        n_rows = n + 3 * n + (1 if gtf else 0) + 1
        ratios = ([0.8] * n + [0.7, 0.7, 3.0] * n
                  + ([gene_ratio] if gtf else []) + [1.2])
    else:
        n_rows = n + 3 * n + (1 if gtf else 0)
        ratios = ([0.8] * n + [0.7, 0.7, 3.2] * n
                  + ([gene_ratio] if gtf else []))
    fig, axes = plt.subplots(
        n_rows, 1,
        figsize=(11, 1.2 + n * 1.5 + 0.12 * min(max_fibers, 90)
                 + (0.35 * n_gene_rows if gtf else 0)),
        gridspec_kw={"hspace": 0.10, "height_ratios": ratios})

    data_axes = axes[:-1] if zoom else axes
    for ax in data_axes:
        ax.set_xlim(0, span)
        ax.set_xticks([])   # no coordinate axes (paper style)

    row = 0
    # --- m5CpG tracks of all samples, stacked at the top of the section ---
    # paper style: discrete brown spikes (rug), not a filled silhouette
    for sname, fibers in samples:
        ax = axes[row]; row += 1
        mc = per_base(fibers, "cpg", start, end)
        pos = [x for x in range(span) if mc[x] > 0]
        ax.vlines(pos, 0, [mc[x] for x in pos], color=MCPG_FILL, lw=0.7)
        ax.set_ylim(0, max(max(mc) * 1.1, 1))
        ax.set_yticks([])
        ax.set_ylabel(f"{sname}\nm5CpG", fontsize=6.5, rotation=0,
                      ha="right", va="center")
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)

    # --- per-sample blocks: FIRE score, coverage, single fibers ---
    last_fiber_ax = None
    stacks = []   # (sname, fiber axes, fiber list) for the optional zoom
    for sname, fibers in samples:
        # FIRE score: black filled blocks of the true step signal (no
        # smoothing: the real S_g is piecewise constant over scored
        # intervals; the paper's sloped silhouettes come from its own data
        # varying per base, not from a rendering trick)
        ax = axes[row]; row += 1
        fs = fire_score_track(fibers, start, end)
        ax.fill_between(range(span), fs, step="mid", color=FIRE_BAR, lw=0)
        ax.set_ylim(0, max(max(fs) * 1.15, 5))
        ax.set_yticks([])
        ax.set_ylabel("FIRE\nscore", fontsize=7, rotation=0, ha="right",
                      va="center")
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)

        # fiber coverage: opaque stacked element coverages (paper pixels:
        # pure colors inside each band, ~1px boundaries — red bottom,
        # purple middle, grey top; the grey band is nucleosome coverage)
        ax = axes[row]; row += 1
        nuc = per_base(fibers, "nuc", start, end)
        lnk = per_base(fibers, "lnk", start, end)
        fir = per_base(fibers, "fire", start, end)
        # draw the widest band first so inner layers stay visible
        ax.fill_between(range(span),
                        [f + l + u for f, l, u in zip(fir, lnk, nuc)],
                        color=COV_NUC, lw=0)
        ax.fill_between(range(span), [f + l for f, l in zip(fir, lnk)],
                        color=COV_LNK, lw=0)
        ax.fill_between(range(span), fir, color=COV_FIRE, lw=0)
        tot = [f + l + u for f, l, u in zip(fir, lnk, nuc)]
        ax.set_ylim(0, max(max(tot) * 1.1, 1.0))
        ax.set_yticks([])
        ax.set_ylabel("coverage", fontsize=7, rotation=0, ha="right",
                      va="center")
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)

        # annotate interior zero-coverage gaps >=500bp (deletions)
        gap_start = None
        for i, v in enumerate(tot):
            if v < 0.5:
                if gap_start is None:
                    gap_start = i
            else:
                if gap_start is not None and i - gap_start >= 500:
                    mid = (gap_start + i) / 2
                    ax.annotate(
                        f"deletion {((i - gap_start) / 1000):.1f} kb",
                        xy=(mid, max(max(tot) * 0.15, 0.3)),
                        fontsize=6, color=GAP_ANN, ha="center", va="bottom")
                gap_start = None

        # --- single fibers: black line, purple linkers, tiered FIRE MSPs ---
        ax = axes[row]; row += 1
        last_fiber_ax = ax

        def inwin(rec):
            return sum(min(b, end) - max(a, start)
                       for a, b in rec["segs"] if a < end and b > start)
        fib_list = sorted(fibers, key=inwin, reverse=True)[:max_fibers]
        stacks.append((sname, ax, fib_list))

        for i, rec in enumerate(fib_list):
            for a, b in bridge_segs(rec["segs"]):
                lo, hi = max(a, start), min(b, end)
                if hi <= lo:
                    continue
                ax.plot([lo - start, hi - start], [i, i], lw=0.9,
                        color=NUC_LINE, zorder=1, solid_capstyle="butt")
            linkers, fires = split_msp_by_fire(rec["msp"], rec["fire_q"])
            for a, b in linkers:           # purple rectangles = linker DNA
                if b < start or a > end:
                    continue
                ax.add_patch(Rectangle(
                    (max(a, start) - start, i - 0.25),
                    min(b, end) - max(a, start), 0.50,
                    facecolor=LINKER, ec="none", zorder=2))
            for a, b, q in fires:          # FIRE MSPs by q tier
                if b < start or a > end:
                    continue
                ax.add_patch(Rectangle(
                    (max(a, start) - start, i - 0.28),
                    min(b, end) - max(a, start), 0.56,
                    facecolor=fire_tier_color(q), ec="none", zorder=3))
        ax.set_ylim(-0.6, len(fib_list) - 0.2)
        ax.set_yticks([])
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
        # (no ylabel here: the green label box already names the stack —
        # a "fibers" ylabel collided with it and only its tail showed)

        # big bracket + green label box at the left of the fiber stack
        ax.plot([-0.012, -0.012], [0.0, 1.0], transform=ax.transAxes,
                color=NUC_LINE, lw=1.0, clip_on=False, zorder=5)
        ax.plot([-0.012, -0.004], [0.0, 0.0], transform=ax.transAxes,
                color=NUC_LINE, lw=1.0, clip_on=False, zorder=5)
        ax.plot([-0.012, -0.004], [1.0, 1.0], transform=ax.transAxes,
                color=NUC_LINE, lw=1.0, clip_on=False, zorder=5)
        ax.text(-0.022, 0.5, f"{sname}\nSingle Fibers", transform=ax.transAxes,
                fontsize=6.5, ha="right", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.35", facecolor=LABEL_BOX,
                          edgecolor="#7A9A7A", lw=0.6))

    # --- optional gene-structure track (--gtf): one stacked row per gene
    # model overlapping the locus; the zoom's grey key shadow below extends
    # across this track so the selection lines up with the gene features.
    if gtf:
        ax = axes[row]; row += 1
        draw_gene_track(ax, genes, start, end, span)

    # --- optional zoom (--zoom): ONE selection for the whole figure, paper
    # style - a light-grey shadow (the key region, ~1/3 of the zoomed width)
    # crosses all tracks, the representative fiber row carries a light-blue
    # box exactly the size of the zoomed region, dashed lines run from the
    # shadow's two endpoints at the stack bottom to the shadow's two edges
    # inside the zoom row (shadow -> shadow), and the shadow continues
    # inside the zoom row.
    if zoom:
        best = None
        for sname_s, axf, flist in stacks:
            for ridx, rec in enumerate(flist):
                nf = sum(1 for a, b, _q in rec["fire_q"]
                         if a < end and b > start)
                iw = sum(min(b, end) - max(a, start)
                         for a, b in rec["segs"] if a < end and b > start)
                if best is None or (nf, iw) > best[0]:
                    best = ((nf, iw), sname_s, axf, ridx, rec)
        _key, sname, axf, ridx, rec = best
        if zoom_center is not None:
            zmid = min(max(int(zoom_center), start), end)
            if not any(a <= zmid < b for _lbl, fl in samples for rec in fl
                       for a, b in rec["segs"]):
                logger.warning(
                    f"--zoom-center {zoom_center} has no read coverage in any "
                    f"sample (deletion/assembly gap?); the zoom row will be "
                    f"empty — consider omitting --zoom-center for automatic "
                    f"FIRE-based selection")
        elif zoom_center_mode:
            zmid = diff_center(samples, group_labels or [], start, end,
                               zoom_center_mode)
            logger.info(f"zoom-centre ({zoom_center_mode}): "
                        f"{chrom}:{zmid:,}")
        else:
            fires_in = [(a, b, q) for a, b, q in rec["fire_q"]
                        if a < end and b > start]
            if fires_in:
                a0, b0, _q0 = min(fires_in, key=lambda t: t[2])
                zmid = (max(a0, start) + min(b0, end)) // 2
            else:
                zmid = (start + end) // 2
        zw = min(zoom_bp, span)          # zoomed region = the blue box
        z0 = max(start, min(zmid - zw // 2, end - zw))
        z1 = z0 + zw
        sw = max(zw // 3, 50)            # key shadow = 1/3 of the zoom
        s0 = max(z0, zmid - sw // 2)
        s1 = min(z1, s0 + sw)

        # grey key shadow: ONE continuous translucent band over the whole
        # data-axes stack, in figure coordinates.  Per-track axvspan copies
        # broke into row fragments (white hspace gaps) and hid behind the
        # opaque element blocks, so the band read as dashes, not a band.
        invf = fig.transFigure.inverted()
        p_top = invf.transform(data_axes[0].transData.transform(
            (s0 - start, data_axes[0].get_ylim()[1])))
        p_bot = invf.transform(data_axes[-1].transData.transform(
            (s1 - start, data_axes[-1].get_ylim()[0])))
        fig.add_artist(Rectangle(
            (p_top[0], p_bot[1]), p_bot[0] - p_top[0], p_top[1] - p_bot[1],
            facecolor="#B9B9B9", alpha=0.30, lw=0, zorder=5,
            clip_on=False))
        axf.add_patch(Rectangle(         # blue box = the zoom region
            (z0 - start, ridx - 0.55), z1 - z0, 1.1,
            facecolor=HILITE, ec=HILITE_EC, lw=0.7, zorder=0.4))

        # zoom axes: span most of the data-axes width so the row reads as a
        # true magnification (px/bp = track scale * span/zoom_bp; the old
        # half-width inset ate the magnification, so shallow zooms looked
        # unzoomed).  Placed so the shadow's magnified rect sits under the
        # track shadow when the margins allow, keeping the callout short.
        pos = axes[-1].get_position()
        pos_d = data_axes[-1].get_position()
        zoff, zlen_w = z0, z1 - z0
        lane = max(zlen_w // 10, 120)    # narrow lane for the 150 bp bar
        inv0 = fig.transFigure.inverted()
        cs = inv0.transform(data_axes[-1].transData.transform(
            ((s0 + s1) / 2 - start, data_axes[-1].get_ylim()[0])))[0]
        zwidth = 0.85 * pos_d.width
        lo = pos_d.x0 + 0.03 * pos_d.width
        hi = pos_d.x0 + 0.97 * pos_d.width - zwidth
        xc = (zlen_w / 2) / (zlen_w + lane + 260)   # shadow-rect centre frac
        zx0 = min(max(cs - xc * zwidth, lo), hi)
        axes[-1].set_position([zx0, pos.y0, zwidth, pos.height])
        ax = axes[-1]
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlim(0, zlen_w + lane + 260)   # right lane for the 150 bp bar
        ax.set_ylim(-0.8, 1.5)
        ax.add_patch(Rectangle((s0 - z0, -0.8), s1 - s0, 2.3,
                               facecolor="#B9B9B9", alpha=0.25, lw=0,
                               zorder=5, clip_on=False))
        for a, b in bridge_segs(rec["segs"]):
            lo, hi = max(a, z0), min(b, z1)
            if hi > lo:
                ax.plot([lo - zoff, hi - zoff], [0, 0], lw=1.4,
                        color=NUC_LINE, zorder=1, solid_capstyle="butt")
        linkers, fires = split_msp_by_fire(rec["msp"], rec["fire_q"])
        covered: List[Tuple[int, int]] = []
        label_cands: List[Tuple[str, int, int]] = []   # type, lo, hi (abs)
        for a, b in linkers:               # purple rectangles = L
            lo, hi = max(a, z0), min(b, z1)
            if hi <= lo:
                continue
            ax.add_patch(Rectangle((lo - zoff, -0.22), hi - lo, 0.44,
                                   facecolor=LINKER, ec="none", zorder=2))
            covered.append((lo, hi))
            label_cands.append(("L", lo, hi))
        for a, b, q in fires:              # FIRE MSPs = RE (tier color)
            lo, hi = max(a, z0), min(b, z1)
            if hi <= lo:
                continue
            ax.add_patch(Rectangle((lo - zoff, -0.25), hi - lo, 0.50,
                                   facecolor=fire_tier_color(q),
                                   ec="none", zorder=3))
            covered.append((lo, hi))
            label_cands.append(("RE", lo, hi))
        for a, b in bridge_segs(rec["segs"]):   # 'n' over nucleosome gaps
            lo, hi = max(a, z0), min(b, z1)
            if hi <= lo:
                continue
            cur = lo
            for ca, cb in sorted(covered):
                if cb <= cur or ca >= hi:
                    continue
                if ca > cur:
                    label_cands.append(("n", cur, ca))
                cur = max(cur, cb)
            if hi > cur:
                label_cands.append(("n", cur, hi))
        # one label per structure type, on its widest instance: repeating
        # 'L'/'RE'/'n' over every segment cluttered the row and overprinted
        # wherever structures sit close together.  Colliding labels are
        # staggered vertically instead of overlapping.
        inv_px = ax.transData.inverted()
        x_per_px = inv_px.transform((1, 0))[0] - inv_px.transform((0, 0))[0]
        placed: List[Tuple[float, float, int]] = []   # x0, x1, level
        for typ in ("L", "RE", "n"):
            lst = [c for c in label_cands if c[0] == typ]
            if not lst:
                continue
            _, lo, hi = max(lst, key=lambda c: c[2] - c[1])
            cx = (lo + hi) / 2 - zoff
            w = (len(typ) * 0.62 * 7.5 * dpi / 72 + 4) * x_per_px
            cx = min(max(cx, w / 2), zlen_w - w / 2)
            x0, x1 = cx - w / 2, cx + w / 2
            margin = 4 * x_per_px
            level = 0
            while any(lv == level and
                      not (x1 + margin < p0 or x0 - margin > p1)
                      for p0, p1, lv in placed):
                level += 1
            placed.append((x0, x1, level))
            ax.text(cx, 0.42 + 0.34 * level, typ, ha="center", va="bottom",
                    fontsize=7.5, zorder=6)
        # 150 bp scale bar on the SAME level as the fiber line, in its own
        # lane clearly separated from the zoom content
        ax.plot([zlen_w + lane + 20, zlen_w + lane + 170], [0, 0], lw=2,
                color=NUC_LINE)
        ax.text(zlen_w + lane + 95, 0.10, "150 bp", ha="center", va="bottom",
                fontsize=7)
        ax.set_ylabel(f"{sname}\nzoom", fontsize=7, rotation=0, ha="right",
                      va="center")

        # dashed connectors: the grey key shadow's two endpoints at the bottom
        # of the data-track stack -> the same shadow's two edges inside the
        # zoom row (its magnified continuation).  Both anchors sit on the
        # shadow, so the callout never crosses the fiber rows.
        inv = fig.transFigure.inverted()
        axb = data_axes[-1]
        yb = axb.get_ylim()[0]
        pL = inv.transform(axb.transData.transform((s0 - start, yb)))
        pR = inv.transform(axb.transData.transform((s1 - start, yb)))
        ztop = ax.get_ylim()[1]
        zL = inv.transform(ax.transData.transform((s0 - z0, ztop)))
        zR = inv.transform(ax.transData.transform((s1 - z0, ztop)))
        for px, py, zx, zy in ((pL[0], pL[1], zL[0], zL[1]),
                               (pR[0], pR[1], zR[0], zR[1])):
            fig.add_artist(Line2D([px, zx], [py, zy], color="#9A9A9A",
                                  lw=0.8, ls=(0, (4, 3)), zorder=6))

    # (no genomic scale bar here: the paper's browser shots carry none;
    # the 150 bp bars inside the optional zoom rows follow the paper)

    axes[0].set_title(f"{name}  ({chrom}:{start:,}-{end:,})",
                      fontsize=11, loc="center", pad=26)
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=LINKER, label="L linker DNA (MSP)"),
        Rectangle((0, 0), 1, 1, facecolor=FIRE_ORANGE,
                  label="MSP FIRE q>0.05"),
        Rectangle((0, 0), 1, 1, facecolor=FIRE_RED,
                  label="MSP FIRE 0.01<q<=0.05"),
        Rectangle((0, 0), 1, 1, facecolor=FIRE_DKRED,
                  label="MSP FIRE q<=0.01 (RE/P)"),
    ]
    # legend: ONE row, centered below the title, inside the track span
    pos_d0 = data_axes[0].get_position()
    cx = pos_d0.x0 + pos_d0.width / 2
    fig.legend(handles=handles, frameon=False, fontsize=7, ncol=4,
               loc="lower center",
               bbox_to_anchor=(cx, pos_d0.y0 + pos_d0.height + 0.002))
    out = Path(out_base)
    out.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{out_base}.{fmt}", dpi=dpi, bbox_inches="tight",
                    facecolor="white")
        logger.info(f"saved {out_base}.{fmt}")
    plt.close(fig)

def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry: parse arguments, load fibers per sample, render every
    requested figure variant.  Returns 0 on success.

    The sample arguments consume the config's native dicts (JSON), so the
    workflow layer can pass its config values through without inventing an
    intermediate string format.
    """
    p = argparse.ArgumentParser(
        description="Paper-style fiber-seq locus browser (Fig2H scheme).")
    p.add_argument("--locus", required=True, help="chrom:start-end|name")
    p.add_argument("--sample-sams", required=True,
                   help="JSON dict sample_id -> SAM/BAM path (the config-"
                        "native form, passed straight through)")
    p.add_argument("--sample-group-dict", default="{}",
                   help="JSON dict sample_id -> group label shown in the "
                        "figure labels, e.g. '{\"s1\": \"aged\"}'")
    p.add_argument("--out-dir", required=True,
                   help="output directory; writes <name>.<fmt> and "
                        "<name>_zoom.<fmt> per requested variant")
    p.add_argument("--format", dest="formats", action="append", default=None,
                   help="image format (repeatable; png pdf svg ps eps tif "
                        "tiff jpg jpeg pgf raw rgba)")
    p.add_argument("--variant", action="append", default=None,
                   choices=["plain", "zoom"],
                   help="which figure variants to write (repeatable)")
    p.add_argument("--max-fibers", type=int, default=12)
    p.add_argument("--zoom-bp", type=int, default=2100,
                   help="width of the zoomed region in bp (the light-blue "
                        "box; the grey key shadow is 1/3 of it)")
    p.add_argument("--gtf", default=None,
                   help="optional GTF/GTF.gz: adds a gene-structure track "
                        "(exon boxes, strand arrows, names) at the bottom of "
                        "every figure; the zoom's grey key shadow extends "
                        "across it")
    p.add_argument("--zoom-center", type=int, default=None,
                   help="genomic coordinate the blue zoom box (and its grey "
                        "key shadow) is centered on; default: the best "
                        "FIRE call of the representative fiber (locus "
                        "midpoint when the locus has no FIRE); out-of-locus "
                        "values are clamped to the locus")
    p.add_argument("--zoom-center-mode", default=None,
                   choices=["diff-fire", "diff-cpg"],
                   help="automatic zoom-centre rule: 'diff-fire' centres on "
                        "the position of maximal absolute between-group "
                        "FIRE-score difference, 'diff-cpg' on the maximal "
                        "m5CpG difference (exactly two groups in "
                        "--sample-group-dict are required). An explicit "
                        "--zoom-center takes precedence")
    p.add_argument("--sample-order", default=None,
                   help="comma-separated sample_ids fixing the display order "
                        "top to bottom (e.g. 's_youth,s_aged'); samples not "
                        "listed keep their input order below the listed ones")
    p.add_argument("--dpi", type=int, default=200)
    args = p.parse_args(argv)
    m = re.match(r"(\S+):(\d+)-(\d+)\|(.+)", args.locus)
    if not m:
        p.error("--locus must be chrom:start-end|name")
    locus = (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
    formats = args.formats or ["png"]
    bad = [f for f in formats if f not in VALID_FORMATS]
    if bad:
        p.error(f"unsupported formats {bad}")
    setup_logging()
    sams = json.loads(args.sample_sams)
    groups = json.loads(args.sample_group_dict)
    samples_map: Dict[str, tuple] = {}
    for sid, sam_path in sams.items():
        label = f"{sid} ({groups[sid]})" if sid in groups else sid
        fibers = parse_sam(sam_path, *locus[:3])
        nfire = sum(len(r["fire_q"]) for r in fibers)
        logger.info(f"{label}: {len(fibers)} primary fibers, {nfire} FIRE calls")
        samples_map[sid] = (label, fibers)
    order = ([s.strip() for s in args.sample_order.split(",") if s.strip()]
             if args.sample_order else [])
    unknown = [s for s in order if s not in samples_map]
    if unknown:
        p.error(f"--sample-order ids not in --sample-sams: {unknown}")
    ordered = order + [s for s in sams if s not in order]
    samples = [samples_map[s] for s in ordered]
    group_labels = [groups.get(s, s) for s in ordered]
    variants = args.variant or ["plain"]
    out_dir = Path(args.out_dir)
    name = locus[3]
    for variant in variants:
        out_base = str(out_dir / (name if variant == "plain" else f"{name}_zoom"))
        plot_locus(locus, samples, out_base, formats, args.max_fibers, args.dpi,
                   zoom=(variant == "zoom"), zoom_bp=args.zoom_bp,
                   gtf=args.gtf, zoom_center=args.zoom_center,
                   zoom_center_mode=args.zoom_center_mode,
                   group_labels=group_labels)
    return 0

if __name__ == "__main__":
    sys.exit(main())
