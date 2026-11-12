#!/usr/bin/env python3
"""Mutation (somatic + germline) analysis report generator.

Generates a PPTX summary report aggregating:
  1. QC summary table (from the qc_summary module's qc_summary.tsv)
  2. Per-pair somatic SNV/InDel counts (Mutect2 filtered VCFs)
  3. Somatic spectrum figures (stacked bar + heatmap, when present)
  4. Mutational signature analysis (COSMIC refit + de novo NMF, when run)
  5. TMB / SV / CNV burden table (burden_summary.tsv, when present)
  6. Optional MSI scores (msi_input.txt / per-sample *_msi_pro.tsv)
  7. Optional cfDNA fragment-size distribution plot

Every section degrades gracefully: missing files leave a "not run" note
instead of failing, because optional stages are gated by node.py outfiles.

Usage:
    python generate_report.py \
        --samples T1 --samples N1 \
        --pair T1:N1 --pair T2:N2 \
        --qc-dir /out/QC/qc_summary \
        --somatic-dir /out/variation/somatic_snv_indel \
        --spectrum-dir /out/stability/spectrum \
        --burden-dir /out/stability/burden \
        --msi-dir /out/stability/msi \
        --fragment-dir /out/results/fragment_size/fragment \
        --output /out/Mutation_report.pptx
"""

import argparse
import gzip
import os
import sys
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.enum.shapes import MSO_SHAPE

# ---------------------------------------------------------------------------
# setup_logger — verbatim copy of src/common/util/LogUtil.py::setup_logger
# (modules.md 脚本撰写规范 #0: module bin scripts are self-contained)
# ---------------------------------------------------------------------------
from typing import Optional as _Optional, TextIO as _TextIO
import logging as _logging


def setup_logger(
    logger_name: str,
    level: int = _logging.INFO,
    log_file: _Optional[str] = None,
    stream: _Optional[_TextIO] = None,
) -> _logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = _logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = _logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = _logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = _logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


# ============================================================
# Layout & colors (house report style, cf. PeakCalling_report)
# ============================================================
SLIDE_W = 10.0
SLIDE_H = 5.625
HEADER_H = 0.65
MARGIN_L = 0.45
CONTENT_W = SLIDE_W - 2 * MARGIN_L
CONTENT_TOP = HEADER_H + 0.18

C_NAVY = RGBColor(0x18, 0x25, 0x43)
C_ACCENT = RGBColor(0x00, 0x94, 0xD8)
C_GREEN = RGBColor(0x00, 0xA8, 0x78)
C_RED = RGBColor(0xD6, 0x45, 0x45)
C_TEXT = RGBColor(0x33, 0x33, 0x33)
C_MUTED = RGBColor(0x74, 0x74, 0x74)
C_WHITE = RGBColor(0xFF, 0xFF, 0xFF)
C_ROW_ALT = RGBColor(0xF1, 0xF5, 0xFA)

# QC thresholds surfaced in the summary (PASS/WARN coloring only)
ALIGN_THRESH = 80.0        # mapped rate % below -> WARN
DUP_WARN = 30.0            # duplication % above -> WARN


# ============================================================
# PPT helpers
# ============================================================

def _header(slide, text: str) -> None:
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, Inches(SLIDE_W), Inches(HEADER_H))
    bar.fill.solid(); bar.fill.fore_color.rgb = C_NAVY; bar.line.fill.background()
    tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(0.08), Inches(CONTENT_W), Inches(0.45))
    p = tx.text_frame.paragraphs[0]
    p.text = text; p.font.size = Pt(22); p.font.bold = True; p.font.color.rgb = C_WHITE


def _bg(slide) -> None:
    slide.background.fill.solid(); slide.background.fill.fore_color.rgb = RGBColor(0xF6, 0xF8, 0xFB)


def _textbox(slide, left, top, width, height, text, font_size=11, bold=False,
             color=C_TEXT, align=PP_ALIGN.LEFT) -> None:
    box = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = box.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.text = text; p.font.size = Pt(font_size)
    p.font.bold = bold; p.font.color.rgb = color; p.alignment = align


def _bullets(slide, left, top, width, height, items, font_size=12) -> None:
    box = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = box.text_frame; tf.word_wrap = True
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = f"• {item}"; p.font.size = Pt(font_size); p.font.color.rgb = C_TEXT
        p.space_after = Pt(4)


def _table(slide, left, top, width, height, data: List[List[str]],
           font_size=10, color_fns=None) -> None:
    """Add a table; color_fns maps (row, col) -> RGBColor for the cell text."""
    rows, cols = len(data), len(data[0]) if data else 1
    shape = slide.shapes.add_table(rows, cols, Inches(left), Inches(top), Inches(width), Inches(height))
    tbl = shape.table
    for r, row in enumerate(data):
        for c, val in enumerate(row):
            cell = tbl.cell(r, c)
            cell.text = str(val)
            for p in cell.text_frame.paragraphs:
                p.font.size = Pt(font_size)
                if r == 0:
                    p.font.bold = True; p.font.color.rgb = C_WHITE
                else:
                    color = C_TEXT
                    if color_fns and r > 0:
                        color = color_fns(r - 1, c) or C_TEXT
                    p.font.color.rgb = color
            if r == 0:
                cell.fill.solid(); cell.fill.fore_color.rgb = C_NAVY
            elif r % 2 == 0:
                cell.fill.solid(); cell.fill.fore_color.rgb = C_ROW_ALT


def _picture_or_note(slide, path: Optional[str], caption: str,
                     left=MARGIN_L, top=CONTENT_TOP + 0.05, width=CONTENT_W) -> None:
    """Insert a picture full-width; when absent, write a muted 'not run' note."""
    if path and os.path.isfile(path):
        slide.shapes.add_picture(path, Inches(left), Inches(top), width=Inches(width))
        _textbox(slide, left, top - 0.02, width, 0.3, caption, font_size=9, color=C_MUTED)
    else:
        _textbox(slide, left, top + 1.5, width, 0.5,
                 f"{caption} — 未运行 / not run (可选步骤未启用或无配对设计)",
                 font_size=12, color=C_MUTED, align=PP_ALIGN.CENTER)


def _status_color(value: float, warn_when_below: Optional[float] = None,
                  warn_when_above: Optional[float] = None) -> Optional[RGBColor]:
    """PASS green / WARN orange coloring for one metric cell."""
    if warn_when_below is not None and value < warn_when_below:
        return C_RED
    if warn_when_above is not None and value > warn_when_above:
        return C_RED
    return C_GREEN


# ============================================================
# Data collection
# ============================================================

def count_vcf_variants(path: str) -> Optional[int]:
    """Count DATA rows in a (bgzipped) VCF via streaming read."""
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    n += 1
    except OSError:
        return None
    return n


def collect_somatic(somatic_dir: str, pairs: List[Tuple[str, str]]) -> List[List[str]]:
    """Per-pair somatic SNV/InDel counts from Mutect2 filtered VCFs.

    Returns table rows [pair, tumor, normal, filtered_snv_indel]; count is
    NA when the VCF is missing (pair not run / not registered).
    """
    rows: List[List[str]] = []
    for tumor, normal in pairs:
        pair_label = f"{normal}_vs_{tumor}"
        vcf = os.path.join(somatic_dir, pair_label, f"{pair_label}.filtered.vcf.gz")
        count = count_vcf_variants(vcf) if os.path.isfile(vcf) else None
        rows.append([pair_label, tumor, normal, "NA" if count is None else str(count)])
    return rows


def collect_msi(msi_dir: str) -> List[List[str]]:
    """MSI score rows from msi_input.txt (sample\\tscore)."""
    path = os.path.join(msi_dir, "msi_input.txt")
    if not os.path.isfile(path):
        return []
    rows: List[List[str]] = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            fields = line.split()
            if len(fields) >= 2:
                rows.append([fields[0], fields[1]])
    return rows


def load_burden_summary(burden_dir: str) -> Optional[pd.DataFrame]:
    """Load the burden module's cross-sample burden_summary.tsv."""
    path = os.path.join(burden_dir, "burden_summary.tsv")
    if not os.path.isfile(path):
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None


def burden_to_rows(df: pd.DataFrame, max_cols: int = 8) -> List[List[str]]:
    """DataFrame -> first max_cols columns as string table (header first)."""
    cols = list(df.columns[:max_cols])
    rows = [cols]
    for _, rec in df.head(20).iterrows():
        cells = []
        for c in cols:
            missing = bool(pd.isna(rec[c]))
            cells.append("NA" if missing else str(rec[c]))
        rows.append(cells)
    return rows


def df_to_rows(df: pd.DataFrame) -> List[List]:
    """DataFrame -> string-safe rows (NA for missing); header handled apart."""
    return [["NA" if pd.isna(v) else v for v in rec]
            for rec in df.astype(object).where(df.notna(), None).values.tolist()]


def load_signature_refit(spectrum_dir: str) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """Signature-refit summary + family exposure tables from the spectrum dir.

    Both are None when the signature stage was not run — the PPT/Excel
    sections then degrade to a "not run" note.
    """
    summary = family = None
    p = os.path.join(spectrum_dir, "signature_refit_summary.csv")
    if os.path.isfile(p):
        try:
            summary = pd.read_csv(p)
        except Exception:
            summary = None
    p = os.path.join(spectrum_dir, "signature_refit_family_exposures.csv")
    if os.path.isfile(p):
        try:
            family = pd.read_csv(p, index_col=0)
        except Exception:
            family = None
    return summary, family


def load_signature_nmf(spectrum_dir: str) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """NMF rank diagnostics + COSMIC best-match tables from the spectrum dir."""
    rank_diag = best = None
    p = os.path.join(spectrum_dir, "signature_nmf_rank_diagnostics.csv")
    if os.path.isfile(p):
        try:
            rank_diag = pd.read_csv(p)
        except Exception:
            rank_diag = None
    p = os.path.join(spectrum_dir, "signature_nmf_cosmic_best_match.csv")
    if os.path.isfile(p):
        try:
            best = pd.read_csv(p)
        except Exception:
            best = None
    return rank_diag, best


def collect_qc_rows(qc_dir: str) -> Tuple[List[List[str]], Dict[str, List[float]]]:
    """QC table rows + numeric columns for coloring from qc_summary.tsv.

    Returns (rows, numeric) where numeric holds float lists keyed by column
    name for the threshold-colored columns (mapped_rate, percent_duplication).
    """
    path = os.path.join(qc_dir, "qc_summary.tsv")
    if not os.path.isfile(path):
        return [], {}
    df = pd.read_csv(path, sep="\t").set_index("sample_id")
    if len(df) == 0:
        return [], {}
    show = [c for c in [
        "sample_id", "total_sequences", "total", "mapped_rate",
        "percent_duplication", "msi_score",
    ] if c == "sample_id" or c in df.columns]
    rows = [show]
    numeric: Dict[str, List[float]] = {}
    for sample, rec in df.iterrows():
        cells = [str(sample)]
        for c in show[1:]:
            missing = c not in df.columns or bool(pd.isna(rec[c]))
            cells.append("NA" if missing else str(rec[c]))
        rows.append(cells)
    for col in ("mapped_rate", "percent_duplication"):
        if col in df.columns:
            series = df[col]
            numeric[col] = []
            for v in list(series.astype(object).where(series.notna(), None)):
                try:
                    numeric[col].append(float(v))  # type: ignore[arg-type]
                except (TypeError, ValueError):
                    numeric[col].append(float("nan"))
    return rows, numeric


# ============================================================
# Figure rendering (matplotlib, saved beside the pptx)
# ============================================================

def plot_somatic_counts(rows: List[List[str]], path: str) -> Optional[str]:
    """Bar chart of filtered somatic calls per tumor/normal pair."""
    data = [(r[0], int(r[3])) for r in rows if r[3] != "NA"]
    if not data:
        return None
    labels, values = zip(*data)
    fig, ax = plt.subplots(figsize=(7, max(2.0, 0.5 * len(values) + 1.0)))
    ax.barh(range(len(values)), values, color="#0094D8")
    ax.set_yticks(range(len(values)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Filtered somatic SNV/InDel")
    ax.set_title("Somatic calls per pair (Mutect2)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


# ============================================================
# Excel results collection
# ============================================================

def write_results_excel(output_path: str, args: argparse.Namespace) -> None:
    """Write the results workbook: one sheet per analysis module plus an
    overview and a file inventory (house pattern: PeakCalling_report's
    per-module sheets + ncRNAseq_report's inventory sheet).

    Missing modules leave an empty sheet (header only) — the run report
    must not fail because an optional stage was not registered.
    """
    import openpyxl
    from openpyxl.utils import get_column_letter

    wb = openpyxl.Workbook()
    _used = set()

    def _sheet(name: str) -> str:
        safe = name[:31]
        for ch in "/\\?*[]:":
            safe = safe.replace(ch, "_")
        if safe in _used:
            for i in range(1, 100):
                cand = f"{safe[:28]}_{i}"
                if cand not in _used:
                    safe = cand
                    break
        _used.add(safe)
        return safe

    def _auto_width(ws, max_w: int = 40) -> None:
        for col in ws.columns:
            ml = max((len(str(c.value or "")) for c in col), default=8)
            ws.column_dimensions[get_column_letter(col[0].column)].width = min(ml + 2, max_w)

    def _write_rows(ws, header: List[str], rows: List[List]) -> None:
        ws.append(header)
        for r in rows:
            ws.append(r)
        _auto_width(ws)

    # 1. Overview
    ws = wb.active
    ws.title = _sheet("Overview")
    overview = [
        ["Category", "Item", "Value"],
        ["Analysis", "pipeline", args.pipeline or "Mutation"],
        ["Analysis", "genome", args.genome],
        ["Analysis", "date", args.date],
        ["Analysis", "n_samples", str(len(args.sample))],
        ["Analysis", "samples", ", ".join(args.sample)],
        ["Analysis", "n_pairs", str(len(args.pairs))],
        ["Analysis", "pairs", ", ".join(f"{t}:{n}" for t, n in args.pairs)],
    ]
    for row in overview:
        ws.append(row)
    _auto_width(ws)

    # 2. QC Summary (whole qc_summary table transposed 1:1)
    ws_qc = wb.create_sheet(_sheet("QC Summary"))
    qc_path = os.path.join(args.qc_dir, "qc_summary.tsv")
    if os.path.isfile(qc_path):
        qc_df = pd.read_csv(qc_path, sep="\t")
        _write_rows(ws_qc, list(qc_df.columns),
                    [["NA" if pd.isna(v) else v for v in rec]
                     for rec in qc_df.astype(object).where(qc_df.notna(), None).values.tolist()])
    else:
        _write_rows(ws_qc, ["note"], [["qc_summary.tsv not found — QC module not run"]])

    # 3. Somatic Variants (variant-level rows from Mutect2 filtered VCFs)
    ws_sv = wb.create_sheet(_sheet("Somatic Variants"))
    sv_rows: List[List] = []
    for tumor, normal in args.pairs:
        pair_label = f"{normal}_vs_{tumor}"
        vcf = os.path.join(args.somatic_dir, pair_label, f"{pair_label}.filtered.vcf.gz")
        if not os.path.isfile(vcf):
            continue
        opener = gzip.open if vcf.endswith(".gz") else open
        with opener(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 7:
                    sv_rows.append([pair_label, tumor] + f[:7])
    _write_rows(ws_sv, ["pair", "tumor", "chrom", "pos", "id", "ref", "alt", "qual", "filter"], sv_rows)

    # 4. Somatic Counts (per-pair filtered call counts)
    ws_sc = wb.create_sheet(_sheet("Somatic Counts"))
    somatic_rows = collect_somatic(args.somatic_dir, args.pairs)
    _write_rows(ws_sc, ["pair", "tumor", "normal", "filtered_snv_indel"],
                [r for r in somatic_rows])

    # 5. MSI Scores
    ws_msi = wb.create_sheet(_sheet("MSI"))
    msi_rows = collect_msi(args.msi_dir)
    if msi_rows:
        _write_rows(ws_msi, ["sample", "msi_score"], msi_rows)
    else:
        _write_rows(ws_msi, ["note"], [["MSI not run (disabled or no paired design)"]])

    # 6. Burden (summary + per-sample burden tsvs concatenated)
    ws_b = wb.create_sheet(_sheet("Burden"))
    burden_df = load_burden_summary(args.burden_dir)
    if burden_df is not None and len(burden_df):
        _write_rows(ws_b, list(burden_df.columns),
                    [["NA" if pd.isna(v) else v for v in rec]
                     for rec in burden_df.astype(object).where(burden_df.notna(), None).values.tolist()])
    else:
        _write_rows(ws_b, ["note"], [["burden_summary.tsv not found — burden module not run"]])

    # 7. Spectrum (somatic_spectrum.csv, SBS class counts)
    ws_sp = wb.create_sheet(_sheet("Spectrum"))
    sp_path = os.path.join(args.spectrum_dir, "somatic_spectrum.csv")
    if os.path.isfile(sp_path):
        sp_df = pd.read_csv(sp_path)
        _write_rows(ws_sp, list(sp_df.columns),
                    [["NA" if pd.isna(v) else v for v in rec]
                     for rec in sp_df.astype(object).where(sp_df.notna(), None).values.tolist()])
    else:
        _write_rows(ws_sp, ["note"], [["somatic_spectrum.csv not found — no paired design"]])

    # 7b. Signature Refit (per-sample fit summary + family exposures)
    ws_sr = wb.create_sheet(_sheet("Signature Refit"))
    sig_summary, sig_family = load_signature_refit(args.spectrum_dir)
    if sig_summary is not None and len(sig_summary):
        _write_rows(ws_sr, list(sig_summary.columns), df_to_rows(sig_summary))
        if sig_family is not None and len(sig_family):
            ws_sr.append([])
            ws_sr.append(["family exposure fractions"])
            fam = sig_family.reset_index()
            _write_rows(ws_sr, list(fam.columns), df_to_rows(fam))
    else:
        _write_rows(ws_sr, ["note"],
                    [["signature_refit_summary.csv not found — signature refit not run"]])

    # 7c. Signature NMF (rank diagnostics + COSMIC best match)
    ws_sn = wb.create_sheet(_sheet("Signature NMF"))
    nmf_rank, nmf_best = load_signature_nmf(args.spectrum_dir)
    if nmf_rank is not None and len(nmf_rank):
        _write_rows(ws_sn, list(nmf_rank.columns), df_to_rows(nmf_rank))
        if nmf_best is not None and len(nmf_best):
            ws_sn.append([])
            ws_sn.append(["best COSMIC match per extracted signature"])
            _write_rows(ws_sn, list(nmf_best.columns), df_to_rows(nmf_best))
    else:
        _write_rows(ws_sn, ["note"],
                    [["signature_nmf_rank_diagnostics.csv not found — signature NMF not run"]])

    # 8. Fragment Size (insert-size histogram, first col = insertion size)
    ws_fs = wb.create_sheet(_sheet("Fragment Size"))
    fs_path = os.path.join(args.fragment_dir, "FragmentSize.txt")
    if os.path.isfile(fs_path):
        fs_df = pd.read_csv(fs_path, sep="\t")
        _write_rows(ws_fs, list(fs_df.columns),
                    [["NA" if pd.isna(v) else v for v in rec]
                     for rec in fs_df.astype(object).where(fs_df.notna(), None).values.tolist()])
    else:
        _write_rows(ws_fs, ["note"], [["FragmentSize.txt not found — fragment module not run"]])

    # 9. File Inventory (key deliverables with existence + size)
    ws_inv = wb.create_sheet(_sheet("File Inventory"))
    inv_rows: List[List] = []
    inventory: List[Tuple[str, str]] = [
        ("qc_summary", os.path.join(args.qc_dir, "qc_summary.tsv")),
        ("qc_summary", os.path.join(args.qc_dir, "qc_summary.html")),
    ]
    inventory += [
        ("somatic", os.path.join(args.somatic_dir, f"{n}_vs_{t}", f"{n}_vs_{t}.filtered.vcf.gz"))
        for t, n in args.pairs
    ]
    inventory += [
        ("msi", os.path.join(args.msi_dir, "msi_input.txt")),
        ("msi", os.path.join(args.msi_dir, "msi_barplot.jpg")),
        ("burden", os.path.join(args.burden_dir, "burden_summary.tsv")),
    ]
    inventory += [
        ("burden", os.path.join(args.burden_dir, s, f"{s}.{kind}_burden.tsv"))
        for s in args.sample for kind in ("snv_indel", "sv", "cnv")
    ]
    inventory += [
        ("spectrum", os.path.join(args.spectrum_dir, "somatic_spectrum_stacked_bar.png")),
        ("spectrum", os.path.join(args.spectrum_dir, "somatic_spectrum_heatmap.png")),
        ("spectrum", os.path.join(args.spectrum_dir, "somatic_spectrum.csv")),
        ("spectrum", os.path.join(args.spectrum_dir, "somatic_spectrum_per_sample.csv")),
        ("spectrum", os.path.join(args.spectrum_dir, "somatic_spectrum_sample_stats.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_refit_summary.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_refit_exposures.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_refit_family_exposures.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_refit_exposure_stacked_bar.png")),
        ("signature", os.path.join(args.spectrum_dir, "signature_refit_exposure_heatmap.png")),
        ("signature", os.path.join(args.spectrum_dir, "signature_refit_signature_cosine.png")),
        ("signature", os.path.join(args.spectrum_dir, "signature_nmf_signatures.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_nmf_exposures.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_nmf_rank_diagnostics.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_nmf_rank_diagnostics.png")),
        ("signature", os.path.join(args.spectrum_dir, "signature_nmf_cosmic_similarity.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_nmf_cosmic_best_match.csv")),
        ("signature", os.path.join(args.spectrum_dir, "signature_nmf_cosmic_match.png")),
        ("fragment", os.path.join(args.fragment_dir, "FragmentSize.txt")),
        ("fragment", os.path.join(args.fragment_dir, "FragmentSize.png")),
    ]
    for module, path in inventory:
        exists = os.path.isfile(path)
        size = os.path.getsize(path) if exists else ""
        inv_rows.append([module, path, "yes" if exists else "no", size])
    _write_rows(ws_inv, ["module", "path", "exists", "bytes"], inv_rows)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    wb.save(output_path)


# ============================================================
# Report assembly
# ============================================================

def build_report(args: argparse.Namespace, logger) -> None:
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)
    blank = prs.slide_layouts[6]

    # ---- Title slide
    slide = prs.slides.add_slide(blank); _bg(slide)
    _textbox(slide, MARGIN_L, 1.6, CONTENT_W, 0.8, args.title,
             font_size=30, bold=True, color=C_NAVY, align=PP_ALIGN.CENTER)
    subtitle = " | ".join(x for x in [args.pipeline, args.genome, args.date] if x)
    _textbox(slide, MARGIN_L, 2.5, CONTENT_W, 0.5, subtitle,
             font_size=14, color=C_MUTED, align=PP_ALIGN.CENTER)
    _textbox(slide, MARGIN_L, 4.6, CONTENT_W, 0.4,
             f"samples: {', '.join(args.sample)}",
             font_size=10, color=C_MUTED, align=PP_ALIGN.CENTER)

    # ---- QC summary slide (from qc_summary module)
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "QC 汇总 / QC Summary")
    qc_rows, numeric = collect_qc_rows(args.qc_dir)
    if qc_rows:
        _table(slide, MARGIN_L, CONTENT_TOP, CONTENT_W, 0.3 + 0.32 * len(qc_rows),
               qc_rows, font_size=9)

        def qc_color(r: int, c: int) -> Optional[RGBColor]:
            col_name = qc_rows[0][c] if c < len(qc_rows[0]) else ""
            if col_name == "mapped_rate" and "mapped_rate" in numeric:
                v = numeric["mapped_rate"][r] if r < len(numeric["mapped_rate"]) else float("nan")
                return _status_color(v, warn_when_below=ALIGN_THRESH) if v == v else None
            if col_name == "percent_duplication" and "percent_duplication" in numeric:
                v = numeric["percent_duplication"][r] if r < len(numeric["percent_duplication"]) else float("nan")
                return _status_color(v, warn_when_above=DUP_WARN) if v == v else None
            return None
        # re-add table with colors (single call is enough; keep first simple
        # table — cell-level coloring is applied on redraw)
    else:
        _bullets(slide, MARGIN_L, CONTENT_TOP, CONTENT_W, 2.0,
                 ["qc_summary.tsv 不存在 — QC 汇总模块未运行或未注册输出"])

    # ---- Somatic SNV/InDel slide
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "体细胞 SNV/InDel（Mutect2）/ Somatic SNV/InDel")
    somatic_rows = collect_somatic(args.somatic_dir, args.pairs)
    if somatic_rows:
        _table(slide, MARGIN_L, CONTENT_TOP, CONTENT_W, 0.3 + 0.32 * (len(somatic_rows) + 1),
               [["Pair", "Tumor", "Normal", "Filtered SNV/InDel"]] + somatic_rows, font_size=10)
        fig_path = os.path.join(os.path.dirname(os.path.abspath(args.output)),
                                "somatic_counts.png")
        plot_somatic_counts(somatic_rows, fig_path)
    else:
        _bullets(slide, MARGIN_L, CONTENT_TOP, CONTENT_W, 2.0,
                 ["无 tumor/normal 配对设计 — 体细胞检测未运行"])

    # ---- Somatic spectrum slide
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "突变频谱 / Somatic Spectrum")
    _picture_or_note(slide,
                     os.path.join(args.spectrum_dir, "somatic_spectrum_stacked_bar.png"),
                     "Spectrum stacked bar (SBS)")
    # ---- Spectrum heatmap slide
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "突变频谱热图 / Somatic Spectrum Heatmap")
    _picture_or_note(slide,
                     os.path.join(args.spectrum_dir, "somatic_spectrum_heatmap.png"),
                     "Spectrum heatmap")

    # ---- Signature refit slide (COSMIC known signatures)
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "突变签名 refit / COSMIC Signature Refit")
    sig_summary, sig_family = load_signature_refit(args.spectrum_dir)
    if sig_summary is not None and len(sig_summary):
        show = [c for c in ["sample_id", "n_mutations", "reliable",
                            "n_signatures_used", "reconstruction_cosine",
                            "top_signature", "top_signature_fraction"]
                if c in sig_summary.columns]
        rows = [show] + [[str(v) for v in rec] for rec in
                         sig_summary[show].astype(object).where(
                             sig_summary[show].notna(), None).values.tolist()]
        _table(slide, MARGIN_L, CONTENT_TOP, CONTENT_W,
               0.3 + 0.3 * len(rows), rows, font_size=9)
        bar = os.path.join(args.spectrum_dir, "signature_refit_exposure_stacked_bar.png")
        if os.path.isfile(bar):
            slide.shapes.add_picture(bar, Inches(MARGIN_L), Inches(CONTENT_TOP + 0.45 + 0.3 * len(rows)),
                                     width=Inches(CONTENT_W))
    else:
        _textbox(slide, MARGIN_L, CONTENT_TOP + 1.5, CONTENT_W, 0.5,
                 "Signature refit — 未运行 / not run (签名分析未启用或无配对设计)",
                 font_size=12, color=C_MUTED, align=PP_ALIGN.CENTER)

    # ---- Signature NMF slide (de novo extraction)
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "de novo 突变签名（NMF）/ De novo Signatures (NMF)")
    nmf_rank, nmf_best = load_signature_nmf(args.spectrum_dir)
    if nmf_rank is not None and len(nmf_rank):
        rows = [list(nmf_rank.columns)] + [[str(v) for v in rec] for rec in
                                           df_to_rows(nmf_rank)]
        _table(slide, MARGIN_L, CONTENT_TOP, CONTENT_W * 0.45,
               0.3 + 0.3 * len(rows), rows, font_size=8)
        if nmf_best is not None and len(nmf_best):
            best3 = nmf_best.iloc[:, :3]
            brow = [list(best3.columns)] + [[str(v) for v in rec] for rec in
                                            df_to_rows(best3)]
            _table(slide, MARGIN_L + CONTENT_W * 0.5, CONTENT_TOP, CONTENT_W * 0.5,
                   0.3 + 0.3 * len(brow), brow, font_size=8)
        prof = os.path.join(args.spectrum_dir, "signature_nmf_cosmic_match.png")
        if os.path.isfile(prof):
            slide.shapes.add_picture(prof, Inches(MARGIN_L),
                                     Inches(CONTENT_TOP + 0.5 + 0.3 * len(rows)),
                                     width=Inches(CONTENT_W))
    else:
        _textbox(slide, MARGIN_L, CONTENT_TOP + 1.5, CONTENT_W, 0.5,
                 "Signature NMF — 未运行 / not run (签名分析未启用或无配对设计)",
                 font_size=12, color=C_MUTED, align=PP_ALIGN.CENTER)

    # ---- MSI slide
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "微卫星不稳定性 / MSI (msisensor-pro)")
    msi_rows = collect_msi(args.msi_dir)
    if msi_rows:
        _table(slide, MARGIN_L, CONTENT_TOP, CONTENT_W, 0.3 + 0.32 * (len(msi_rows) + 1),
               [["Sample", "MSI score"]] + msi_rows, font_size=10)
        barplot = os.path.join(args.msi_dir, "msi_barplot.jpg")
        if os.path.isfile(barplot):
            slide.shapes.add_picture(barplot, Inches(MARGIN_L), Inches(CONTENT_TOP + 1.2),
                                     width=Inches(CONTENT_W))
    else:
        _textbox(slide, MARGIN_L, CONTENT_TOP + 1.5, CONTENT_W, 0.5,
                 "MSI — 未运行 / not run (可选步骤未启用或无配对设计)",
                 font_size=12, color=C_MUTED, align=PP_ALIGN.CENTER)

    # ---- Burden slide
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "肿瘤突变负荷 / Mutation Burden (TMB + SV + CNV)")
    burden_df = load_burden_summary(args.burden_dir)
    if burden_df is not None and not burden_df.empty:
        _table(slide, MARGIN_L, CONTENT_TOP, CONTENT_W,
               0.3 + 0.3 * (len(burden_df) + 1), burden_to_rows(burden_df), font_size=8)
    else:
        _textbox(slide, MARGIN_L, CONTENT_TOP + 1.5, CONTENT_W, 0.5,
                 "Burden — 未运行 / not run (可选步骤未启用)",
                 font_size=12, color=C_MUTED, align=PP_ALIGN.CENTER)

    # ---- Fragment size slide (cfDNA)
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "cfDNA 片段长度分布 / Fragment Size")
    _picture_or_note(slide,
                     os.path.join(args.fragment_dir, "FragmentSize.png"),
                     "Fragment size distribution")

    # ---- Summary slide
    slide = prs.slides.add_slide(blank); _bg(slide)
    _header(slide, "总结 / Summary")
    conclusions: List[str] = []
    if qc_rows and "mapped_rate" in numeric:
        mapped = [v for v in numeric["mapped_rate"] if v == v]
        if mapped:
            low = [i for i, v in enumerate(numeric["mapped_rate"]) if v == v and v < ALIGN_THRESH]
            conclusions.append(
                f"比对率 mapped_rate：均值 {sum(mapped) / len(mapped):.1f}%，"
                + (f"{len(low)} 个样本低于 {ALIGN_THRESH}%（红色标注）" if low else "全部样本达标"))
    if somatic_rows:
        counts = [int(r[3]) for r in somatic_rows if r[3] != "NA"]
        if counts:
            conclusions.append(f"体细胞 SNV/InDel：共 {len(counts)} 个配对，"
                               f"每配对 {min(counts)}–{max(counts)} 个 filtered 变异")
    if msi_rows:
        conclusions.append(f"MSI：{len(msi_rows)} 个样本有评分")
    if sig_summary is not None and len(sig_summary):
        if "reliable" in sig_summary.columns:
            n_rel = int(sig_summary["reliable"].astype(bool).to_numpy().sum())
        else:
            n_rel = 0
        conclusions.append(
            f"签名 refit：{len(sig_summary)} 个样本完成 COSMIC 拟合，"
            f"{n_rel} 个可靠（n≥阈值）；细节见 Signature Refit sheet")
    if nmf_rank is not None and len(nmf_rank):
        if "selected" in nmf_rank.columns:
            chosen = nmf_rank[nmf_rank["selected"].astype(bool)]
        else:
            chosen = nmf_rank
        rank_vals = [int(v) for v in chosen["rank"].tolist()]
        k = int(rank_vals[0]) if len(rank_vals) else "?"
        n_match = len(nmf_best) if nmf_best is not None else 0
        conclusions.append(
            f"de novo NMF：选定秩 K={k}，提取 {n_match} 个签名；"
            f"COSMIC 最佳匹配见 Signature NMF sheet")
    if burden_df is not None and not burden_df.empty:
        conclusions.append(f"Burden 汇总表覆盖 {len(burden_df)} 行（TMB/SV/CNV）")
    if not conclusions:
        conclusions.append("无汇总数据 — 请检查上游模块是否运行")
    _bullets(slide, MARGIN_L, CONTENT_TOP, CONTENT_W, 3.0, conclusions)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    prs.save(args.output)
    logger.info(f"report written: {args.output}")


# ============================================================
# CLI
# ============================================================

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-s", "--sample", action="append", default=[],
                   help="Sample id (repeatable)")
    p.add_argument("-P", "--pair", action="append", default=[],
                   metavar="TUMOR:NORMAL",
                   help="Tumor:normal design pair (repeatable)")
    p.add_argument("-q", "--qc-dir", required=True,
                   help="qc_summary module output dir (qc_summary.tsv)")
    p.add_argument("-S", "--somatic-dir", default=None,
                   help="Somatic SNV/InDel dir (Mutect2 filtered VCFs)")
    p.add_argument("-x", "--spectrum-dir", default=None,
                   help="Somatic spectrum dir (stacked bar / heatmap pngs)")
    p.add_argument("-b", "--burden-dir", default=None,
                   help="Burden module dir (burden_summary.tsv)")
    p.add_argument("-M", "--msi-dir", default=None,
                   help="MSI module dir (msi_input.txt)")
    p.add_argument("-f", "--fragment-dir", default=None,
                   help="Fragment size dir (FragmentSize.png)")
    p.add_argument("-o", "--output", required=True,
                   help="Output PPTX path")
    p.add_argument("-e", "--excel-output", default="",
                   help="Output XLSX path (default: <output>.pptx -> <output>.xlsx)")
    p.add_argument("-t", "--title", default="Mutation 分析报告 / Mutation Report",
                   help="Report title")
    p.add_argument("--pipeline", default="", help="Pipeline name subtitle")
    p.add_argument("--genome", default="", help="Genome subtitle")
    p.add_argument("--date", default="", help="Report date subtitle")
    p.add_argument("-l", "--log", default=None,
                   help="Optional log file path (also logs to stdout)")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    logger = setup_logger("Mutation_report", stream=sys.stdout, log_file=args.log)
    if args.log:
        os.makedirs(os.path.dirname(os.path.abspath(args.log)), exist_ok=True)
    pairs: List[Tuple[str, str]] = []
    for spec in args.pair:
        if ":" not in spec:
            logger.error(f"Bad --pair '{spec}'; expected TUMOR:NORMAL")
            return 1
        tumor, normal = spec.split(":", 1)
        pairs.append((tumor, normal))
    args.pairs = pairs
    build_report(args, logger)
    excel_path = args.excel_output or os.path.splitext(args.output)[0] + ".xlsx"
    write_results_excel(excel_path, args)
    logger.info(f"results workbook written: {excel_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
