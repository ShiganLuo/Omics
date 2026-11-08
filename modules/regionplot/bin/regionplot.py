#!/usr/bin/env python3
"""Render region-focused multi-sample track stacking plots with pyGenomeTracks.

Given a regions file (BED rows or a ``[NAME=]chr:start-end`` list, formats
may be mixed), one or more track groups (label + per-sample bigWig paths)
and an optional group map (sample -> biological group), writes a
pyGenomeTracks ``.ini`` per region panel and calls ``pyGenomeTracks`` to
produce PNG/SVG stacking plots: one wiggle track line per sample, samples of
the same biological group drawn in the same color, a group color legend on
top, and an optional gene-model track from GTF at the bottom.

The script is self-contained (module bin convention): logging helper is
copied in, no ``src.common`` imports.

Usage:
    python regionplot.py \
        --regions regions.bed \
        --track-tsv tracks.tsv \
        --group-map group_map.tsv \
        --gtf genes.gtf --gtf-track \
        --outdir plots/panel \
        --panel-name MYC_locus \
        --formats png svg \
        --track-height 1.5 --width 12
"""
from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# pyGenomeTracks plotting runs headless inside containers
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# ----------------------------------------------------------------------------
# Logging (behavior-compatible copy of src/common/util/LogUtil.setup_logger;
# module bin scripts must stay self-contained)
# ----------------------------------------------------------------------------


def setup_logger(logger_name: str, log_file: Optional[str] = None) -> logging.Logger:
    """Return a logger with the repo's standard format.

    Parameters
    ----------
    logger_name : str
        Logger name (usually the rule name).
    log_file : str, optional
        Also write log records to this path (created with parents).

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")
    if not logger.handlers:
        sh = logging.StreamHandler()
        sh.setFormatter(formatter)
        logger.addHandler(sh)
    if log_file:
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    return logger


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    """Logger per repo convention; console stream redirected to stdout."""
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("regionplot", log_file=str(log_file))
    else:
        lg = setup_logger("regionplot")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


# ----------------------------------------------------------------------------
# Input parsing
# ----------------------------------------------------------------------------

def _sanitize_panel_name(name: str) -> str:
    """Make a region name safe as a path element (panel directory name)."""
    for ch in ("/", " ", ":", "="):
        name = name.replace(ch, "_")
    return name


_REGION_RE = re.compile(r"([^:]+):([\d,]+)-([\d,]+)")


def _parse_region_string(s: str) -> Tuple[str, int, int, str]:
    """Parse one ``[NAME=]chrom:start-end`` region string.

    Parameters
    ----------
    s : str
        Region string, e.g. ``MYC_locus=chr8:127735000-127745000`` or
        ``chr8:127735000-127745000`` (name defaults to the string with
        ``:``/``-``/``,`` replaced by ``_``). Thousands separators
        (``1,277,350``) are tolerated and stripped.

    Returns
    -------
    (chrom, start, end, name)
        Coordinates passed to pyGenomeTracks verbatim (same 0-based
        half-open window a BED row ``chrom start end`` would plot).

    Raises
    ------
    ValueError
        On a malformed string or end <= start.
    """
    name = s
    region = s
    if "=" in s:
        name, region = s.split("=", 1)
    match = _REGION_RE.fullmatch(region.strip())
    if match is None:
        raise ValueError(
            f"malformed region {s!r}, expected [NAME=]chrom:start-end")
    chrom = match.group(1)
    start = int(match.group(2).replace(",", ""))
    end = int(match.group(3).replace(",", ""))
    if end <= start:
        raise ValueError(f"region end <= start: {s!r}")
    if name == s:
        # default name mirrors the BED default chrom_start_end shape
        name = f"{chrom}_{start}_{end}"
    return chrom, start, end, _sanitize_panel_name(name)


def parse_regions(path: Path) -> List[Tuple[str, int, int, str]]:
    """Read regions from a BED-style or ``chr:start-end`` list file.

    Each non-comment row is auto-detected: tab-separated rows with >= 3
    columns parse as BED (chrom/start/end[/name]); anything else parses
    as one ``[NAME=]chrom:start-end`` string. Formats may be mixed.

    Parameters
    ----------
    path : Path
        Regions file (BED or chr:start-end list; ``#``/``track``/
        ``browser`` lines skipped).

    Returns
    -------
    list of (chrom, start, end, panel_name)
        panel_name is the BED col4 / ``NAME=`` prefix when present else
        ``chrom_start_end``. Duplicate names get ``.1``/``.2`` suffixes
        so panel wildcards stay unique.

    Raises
    ------
    ValueError
        On a malformed row or zero regions.
    """
    regions: List[Tuple[str, int, int, str]] = []
    seen: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            if "\t" in line:
                f = line.split("\t")
                if len(f) < 3:
                    raise ValueError(f"regions BED row needs >=3 columns: {line!r}")
                chrom, start, end = f[0], int(f[1]), int(f[2])
                if end <= start:
                    raise ValueError(f"region end <= start: {line!r}")
                name = (f[3].strip() if len(f) >= 4 and f[3].strip()
                        else f"{chrom}_{start}_{end}")
                name = _sanitize_panel_name(name)
            else:
                chrom, start, end, name = _parse_region_string(line.strip())
            if name in seen:
                seen[name] += 1
                name = f"{name}.{seen[name]}"
            else:
                seen[name] = 0
            regions.append((chrom, start, end, name))
    if not regions:
        raise ValueError(f"no regions found in {path}")
    return regions


def parse_track_tsv(path: Path) -> List[Dict[str, str]]:
    """Read track group definitions.

    Parameters
    ----------
    path : Path
        TSV written by the rule with header
        ``label<TAB>sample_id<TAB>bw_path<TAB>order``; rows are expected
        pre-sorted by the rule (group-major, then sample order).

    Returns
    -------
    list of dict
        One dict per row with keys label/sample_id/bw_path/order.
    """
    rows: List[Dict[str, str]] = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:3] != ["label", "sample_id", "bw_path"]:
            raise ValueError(f"unexpected track TSV header: {header}")
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rows.append({
                "label": f[0],
                "sample_id": f[1],
                "bw_path": f[2],
                "order": f[3] if len(f) > 3 else "",
            })
    return rows


def parse_group_map(path: Optional[Path]) -> Dict[str, str]:
    """Read optional sample -> group map (header sample_id<TAB>group)."""
    groups: Dict[str, str] = {}
    if path is None or not Path(path).exists():
        return groups
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].strip():
                groups[f[0]] = f[1].strip()
    return groups


# ----------------------------------------------------------------------------
# Group colors and legend
# ----------------------------------------------------------------------------

# tab10 base palette; group members get lightness-ramped shades of the
# group color so within-group samples stay visually related but separable.
_BASE_COLORS = [
    "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
    "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f",
]


def _shade(color: str, factor: float) -> str:
    """Lighten (factor>0) or darken (factor<0) a hex color.

    Parameters
    ----------
    color : str
        Hex color like ``#1f77b4``.
    factor : float
        Shift amount in [-1, 1]; 0 returns the color unchanged.

    Returns
    -------
    str
        Shifted hex color.
    """

    def clip(v: int) -> int:
        return max(0, min(255, v))

    r = int(color[1:3], 16)
    g = int(color[3:5], 16)
    b = int(color[5:7], 16)
    if factor >= 0:
        r, g, b = (clip(int(v + (255 - v) * factor)) for v in (r, g, b))
    else:
        r, g, b = (clip(int(v * (1 + factor))) for v in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"


def assign_group_colors(rows: List[Dict[str, str]], group_map: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Assign one color per biological group and per (group, sample).

    Parameters
    ----------
    rows : list of dict
        Track rows (label/sample_id/bw_path/order).
    group_map : dict
        sample_id -> group name; missing samples get no explicit color.

    Returns
    -------
    (group_colors, track_colors)
        group_colors: group -> hex (for the legend); track_colors:
        ``(group, sample_id)`` -> hex (group base shaded per member index).
    """
    # preserve first-appearance order of groups (rows are pre-sorted)
    group_order: List[str] = []
    member_index: Dict[Tuple[str, str], int] = {}
    member_count: Dict[str, int] = {}
    for row in rows:
        grp = group_map.get(row["sample_id"])
        if grp is None:
            continue
        if grp not in member_count:
            member_count[grp] = 0
            group_order.append(grp)
        key = (grp, row["sample_id"])
        if key not in member_index:
            member_index[key] = member_count[grp]
            member_count[grp] += 1
    group_colors: Dict[str, str] = {
        grp: _BASE_COLORS[i % len(_BASE_COLORS)] for i, grp in enumerate(group_order)
    }
    track_colors: Dict[str, str] = {}
    for (grp, sid), idx in member_index.items():
        # spread members within a group between -15% and +30% lightness
        n = max(member_count[grp] - 1, 1)
        factor = -0.15 + 0.45 * idx / n
        track_colors[f"{grp}\t{sid}"] = _shade(group_colors[grp], factor)
    return group_colors, track_colors


def write_legend(group_colors: Dict[str, str], out_png: Path, out_svg: Optional[Path]) -> None:
    """Write the group color legend as a standalone figure.

    The legend is composed on top of the pyGenomeTracks figure by the
    caller (SVG/PNG passthrough); a standalone copy is delivered for
    reuse in figures.
    """
    if not group_colors:
        return
    fig, ax = plt.subplots(figsize=(6, 0.4 + 0.35 * len(group_colors)))
    ax.axis("off")
    handles = [
        Line2D([0], [0], marker="s", color="none",
                   markerfacecolor=color, markersize=10, label=grp)
        for grp, color in group_colors.items()
    ]
    ax.legend(handles=handles, loc="center left", frameon=False,
              title="Groups", title_fontsize=10, fontsize=9)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    if out_svg:
        fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# pyGenomeTracks ini + rendering
# ----------------------------------------------------------------------------

def build_ini(rows: List[Dict[str, str]], group_map: Dict[str, str],
              track_colors: Dict[str, str], gtf: Optional[Path],
              track_height: float, width_in: float, out_ini: Path) -> List[Tuple[str, str, str]]:
    """Write the pyGenomeTracks tracks.ini for one region panel.

    Parameters
    ----------
    rows, group_map, track_colors, gtf, region, track_height, width_in
        As described in :func:`render_panel`.
    out_ini : Path
        Destination ini path.

    Returns
    -------
    list of (sample_id, track_group, bw_path)
        The audit rows also written to track_sources.tsv by the caller.
    """
    # group-major ordering: preserve the rule's pre-sorted row order
    lines: List[str] = ["[spacer]\nheight = 0.05\n"]
    audit: List[Tuple[str, str, str]] = []
    for row in rows:
        sid = row["sample_id"]
        grp = group_map.get(sid)
        label = f"{sid} ({grp})" if grp else sid
        color = track_colors.get(f"{grp}\t{sid}", "#1f77b4") if grp else "#1f77b4"
        lines.append(f"[{sid}:{row['label']}]\n")
        lines.append(f"file = {row['bw_path']}\n")
        lines.append(f"title = {label}\n")
        lines.append(f"color = {color}\n")
        lines.append("height = {:.2f}\n".format(track_height))
        lines.append("display = lines\n")
        lines.append("min_value = 0\n")
        lines.append("overlap = 0\n")
        lines.append("\n")
        audit.append((sid, row["label"], row["bw_path"]))
    if gtf:
        lines.append("[genes]\n")
        lines.append(f"file = {gtf}\n")
        lines.append("height = 2.0\n")
        lines.append("title = Genes\n")
        lines.append("file_type = gtf\n")
        lines.append("merge_transcripts = true\n")
        lines.append("prefered_name = gene_name\n")
        lines.append("\n")
    lines.append("[x-axis]\n")
    lines.append("fontsize = 10\n")
    out_ini.parent.mkdir(parents=True, exist_ok=True)
    out_ini.write_text("".join(lines))
    return audit


def render_panel(outdir: Path, region: Tuple[str, int, int, str], rows: List[Dict[str, str]],
                 group_map: Dict[str, str], gtf: Optional[Path], formats: List[str],
                 track_height: float, width_in: float) -> List[Path]:
    """Render one region panel.

    Parameters
    ----------
    outdir : Path
        Panel output directory; receives ``<panel>.<fmt>`` images,
        ``tracks.ini``, ``track_sources.tsv`` and ``group_map.tsv``.
    region : (chrom, start, end, panel_name)
    rows : track rows (pre-sorted, group-major).
    group_map, gtf, formats, track_height, width_in
        As in the CLI.

    Returns
    -------
    list of Path
        Every file written (all declared rule outputs must appear here).
    """
    chrom, start, end, panel = region
    written: List[Path] = []
    group_colors, track_colors = assign_group_colors(rows, group_map)

    ini = outdir / "tracks.ini"
    audit = build_ini(rows, group_map, track_colors, gtf,
                      track_height, width_in, ini)
    written.append(ini)

    # audit table: which bw fed which track line
    tsv = outdir / "track_sources.tsv"
    with open(tsv, "w") as fh:
        fh.write("sample_id\ttrack_group\tbw_path\n")
        for sid, grp, bw in audit:
            fh.write(f"{sid}\t{grp}\t{bw}\n")
    written.append(tsv)

    if group_colors:
        gm = outdir / "group_map.tsv"
        with open(gm, "w") as fh:
            fh.write("sample_id\tgroup\n")
            for sid in sorted(group_map):
                fh.write(f"{sid}\t{group_map[sid]}\n")
        written.append(gm)
        legend_png = outdir / "group_legend.png"
        write_legend(group_colors, legend_png, outdir / "group_legend.svg" if "svg" in formats else None)
        written.append(legend_png)
        if "svg" in formats:
            written.append(outdir / "group_legend.svg")

    region_str = f"{chrom}:{start}-{end}"
    exe = shutil.which("pyGenomeTracks") or "pyGenomeTracks"
    for fmt in formats:
        img = outdir / f"{panel}.{fmt}"
        cmd = [exe, "--tracks", str(ini), "--region", region_str,
               "--outFileName", str(img),
               "--width", str(width_in),
               "--dpi", "150",
               "--trackLabelFraction", "0.22"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(
                f"pyGenomeTracks failed for {region_str}: {res.stderr.strip()[:2000]}")
        if not img.exists():
            raise RuntimeError(f"pyGenomeTracks produced no output for {region_str}: {img}")
        written.append(img)
        logger.info(f"rendered {img.name} for {region_str} ({len(rows)} tracks)")
    return written


def run(args: argparse.Namespace) -> int:
    regions = parse_regions(Path(args.regions))
    if args.panel:
        regions = [r for r in regions if r[3] == args.panel]
        if not regions:
            logger.error(f"panel {args.panel!r} not found in {args.regions}")
            return 1
    rows = parse_track_tsv(Path(args.track_tsv))
    group_map = parse_group_map(Path(args.group_map) if args.group_map else None)

    # verify all bw inputs exist up front: pgT errors are cryptic otherwise
    missing = [r["bw_path"] for r in rows if not Path(r["bw_path"]).exists()]
    if missing:
        logger.error(f"missing bigWig inputs: {missing}")
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    for region in regions:
        render_panel(args.outdir, region, rows, group_map,
                     Path(args.gtf) if args.gtf and args.gtf_track else None,
                     list(args.formats), args.track_height, args.width)
    logger.info(f"rendered {len(regions)} panel(s) into {args.outdir}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="regionplot.py",
        description="Region-focused multi-sample track stacking plots via pyGenomeTracks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-r", "--regions", type=Path, required=True,
                   help="regions file: BED rows (chrom/start/end[/name]) or "
                        "[NAME=]chrom:start-end list lines; formats may be mixed")
    p.add_argument("--panel", type=str, default=None,
                   help="render only this panel (region name); default renders all regions")
    p.add_argument("-t", "--track-tsv", type=Path, required=True,
                   help="track table: label<TAB>sample_id<TAB>bw_path")
    p.add_argument("-g", "--group-map", type=Path, default=None,
                   help="sample->group map TSV (sample_id<TAB>group); enables group coloring + legend")
    p.add_argument("--gtf", type=str, default=None,
                   help="gene annotation GTF for the gene-model track")
    p.add_argument("--gtf-track", action="store_true",
                   help="add the GTF gene-model track (requires --gtf)")
    p.add_argument("-o", "--outdir", type=Path, required=True,
                   help="output directory for this panel")
    p.add_argument("--formats", nargs="+", default=["png", "svg"],
                   help="output image formats")
    p.add_argument("--track-height", type=float, default=1.5,
                   help="track height in inches per sample track")
    p.add_argument("--width", type=float, default=12, help="figure width in inches")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in [args.regions, args.track_tsv] + ([args.group_map] if args.group_map else []):
        if f is not None and not Path(f).exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.gtf_track and not args.gtf:
        logger.error("--gtf-track requires --gtf")
        return 1
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
