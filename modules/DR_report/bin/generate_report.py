#!/usr/bin/env python3
"""Generate the DR workflow joint analysis report (PPTX + XLSX).

The DR workflow pairs a DNA (WES somatic) chain and an RNA-seq chain on the
same biological samples and adds DR-unique joint analyses: genomic
instability burden (SNV/InDel/SV/CNV), DRclass D×R evidence classification
(DNA-anchor / RNA-anchor labels, penetrance) and TransEffect transcriptional
effect analysis. This report consumes the finished outputs of all chains and
renders one PPTX slide deck plus one XLSX file inventory. It never
recomputes upstream analyses.

Section outputs are best-effort: an optional stage that was not enabled
degrades to a "not run" note on its slide and an empty sheet in the
workbook, never an error.

Usage:
    python3 generate_report.py \
        --analysis-dir /data/.../DR \
        --dna-dir /data/.../DR/DNA --rna-dir /data/.../DR/RNA \
        --drclass-dir /data/.../DR/DRclass --teffect-dir /data/.../DR/TransEffect \
        --burden-dir /data/.../DR/DNA/stability/burden \
        --output DR_report.pptx --file-inventory DR_report_files.xlsx \
        --title "DNA×RNA 联合分析报告" --genome GRCm39 --genome-label GRCm39 \
        --date 2026-11-25 --lang zh --img-dir figures --format png \
        --pair TUMOR:NORMAL --contrast A_vs_B
"""

import argparse
import gzip
import os
import sys
import tempfile
import time
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

# setup_logger copied verbatim from src/common/util/LogUtil (module bin
# scripts are self-contained per modules.md 脚本撰写规范 #0).
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


logger = setup_logger("DR_report", stream=sys.stdout)


PlotFormat = str

# matplotlib savefig's full format set (LogUtil.PlotFormat); anything else
# fails fast at parse time instead of dying inside savefig (modules.md #9)
VALID_FORMATS = ("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg",
                 "jpeg", "pgf", "raw", "rgba")

# ── 上游产出图 raster-first 变体解析(modules.md #9 formats 规范;与
# Mutation_report/RNAseq_report 一致)──
_PICK_EXTS = ("png", "jpg", "jpeg", "tif", "tiff",
              "pdf", "svg", "ps", "eps", "pgf", "raw", "rgba")
_RASTER_EXTS = ("png", "jpg", "jpeg", "tif", "tiff")


def _pick_variant(base: str) -> str:
    """Raster-first existing variant of a consumed producer figure.

    *base* is the conventional path (with extension); the scan probes
    ``<stem>.<fmt>`` in raster-first order and returns the first existing
    file, falling back to *base* when no variant exists.
    """
    stem, ext = os.path.splitext(base)
    if ext.lstrip(".").lower() not in _PICK_EXTS:
        stem = base
    for fmt in _PICK_EXTS:
        cand = "{}.{}".format(stem, fmt)
        if os.path.isfile(cand):
            return cand
    return base


def _embeddable(path: Optional[str]) -> bool:
    """True when *path* exists as a raster PPTX can embed (python-pptx
    embeds rasters only, so never hand vector variants to add_picture)."""
    return (bool(path) and os.path.isfile(path)
            and os.path.splitext(path)[1].lower() in _RASTER_EXTS)


# ============================================================
# I18N (zh default / en)
# ============================================================

I18N = {
    "zh": {
        "title_default": "DNA×RNA 联合分析报告 / DR Joint Analysis Report",
        "overview": "流程总览",
        "dna_qc": "DNA 质控概览",
        "dna_variants": "DNA 变异概览(体细胞 SNV/InDel + MSI)",
        "spectrum": "突变谱与突变签名",
        "fragment": "片段长度分布",
        "burden": "基因组稳定性负荷(Burden)",
        "rna_qc": "RNA 质控概览",
        "te_chimeric": "TE 与嵌合转录本",
        "de": "差异表达概览",
        "func": "功能富集概览",
        "fusion": "融合基因概览",
        "drclass": "DRclass D×R 联合证据分类",
        "teffect": "TransEffect 转录效应分析",
        "conclusion": "结论",
        "not_run": "未运行 / not run (可选步骤未启用或无配对设计)",
        "sample_table": "样本",
        "chain_dna": "DNA 链(体细胞变异 / 稳定性)",
        "chain_rna": "RNA 链(TE / 差异 / 功能 / 融合)",
        "chain_dr": "DR 独特分析(burden / DRclass / TransEffect)",
        "up": "上调", "down": "下调", "contrast": "对比",
        "label": "标签", "penetrance": "穿透率", "sample": "样本",
    },
    "en": {
        "title_default": "DR Joint Analysis Report",
        "overview": "Workflow Overview",
        "dna_qc": "DNA QC Overview",
        "dna_variants": "DNA Variant Overview (somatic SNV/InDel + MSI)",
        "spectrum": "Somatic Spectrum & Mutational Signatures",
        "fragment": "Fragment Size Distribution",
        "burden": "Genomic Instability Burden",
        "rna_qc": "RNA QC Overview",
        "te_chimeric": "TE & Chimeric Transcripts",
        "de": "Differential Expression Overview",
        "func": "Functional Enrichment Overview",
        "fusion": "Fusion Gene Overview",
        "drclass": "DRclass D×R Joint Evidence Classification",
        "teffect": "TransEffect Transcriptional Effect Analysis",
        "conclusion": "Conclusion",
        "not_run": "not run (optional stage disabled or no design pair)",
        "sample_table": "Sample",
        "chain_dna": "DNA chain (somatic variants / stability)",
        "chain_rna": "RNA chain (TE / DE / function / fusion)",
        "chain_dr": "DR-unique analyses (burden / DRclass / TransEffect)",
        "up": "Up", "down": "Down", "contrast": "Contrast",
        "label": "Label", "penetrance": "Penetrance", "sample": "Sample",
    },
}


def _t(key: str, lang: str) -> str:
    return I18N.get(lang, I18N["zh"]).get(key, I18N["zh"].get(key, key))


def _setup_cjk_font() -> None:
    """Best-effort CJK font registration for matplotlib labels."""
    from matplotlib import font_manager
    for name in ("WenQuanYi Zen Hei", "WenQuanYi Micro Hei", "Noto Sans CJK SC",
                 "Source Han Sans SC", "SimHei"):
        try:
            path = font_manager.findfont(name, fallback_to_default=False)
            if path and os.path.isfile(path):
                font_manager.fontManager.addfont(path)
                plt.rcParams["font.sans-serif"] = [name] + plt.rcParams["font.sans-serif"]
                plt.rcParams["axes.unicode_minus"] = False
                return
        except Exception:
            continue


# ============================================================
# Layout & colors (house report style, cf. Mutation_report)
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
# PPT helpers (copied from modules/Mutation_report/bin/generate_report.py)
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


def _picture_or_note(slide, path: Optional[str], caption: str, lang: str,
                     left=MARGIN_L, top=CONTENT_TOP + 0.05, width=CONTENT_W) -> None:
    """Insert a picture full-width; when absent, write a muted 'not run' note."""
    if path and _embeddable(path):
        slide.shapes.add_picture(path, Inches(left), Inches(top), width=Inches(width))
        _textbox(slide, left, top - 0.02, width, 0.3, caption, font_size=9, color=C_MUTED)
    else:
        _textbox(slide, left, top + 1.5, width, 0.5,
                 f"{caption} — {_t('not_run', lang)}",
                 font_size=12, color=C_MUTED, align=PP_ALIGN.CENTER)


def _add_picture(slide, path: str, left: float, top: float,
                 max_w: float, max_h: float) -> None:
    """Insert a picture scaled to fit (max_w x max_h), keeping aspect ratio."""
    from PIL import Image as PILImage  # python-pptx 依赖 Pillow,惰性导入
    img = PILImage.open(path)
    aspect = img.size[0] / max(img.size[1], 1)
    width = max_w
    height = width / aspect
    if height > max_h:
        height = max_h
        width = height * aspect
    slide.shapes.add_picture(path, Inches(left), Inches(top),
                             Inches(width), Inches(height))


# Bottom margin floor for content (report-module-creation table sizing):
# y + table height must stay below this line.
_SLIDE_BOTTOM = SLIDE_H - 0.22


# ============================================================
# Figure rendering (matplotlib; embedded raster + figures/ 导出)
# ============================================================

# 报告自有图除嵌入 PPTX 外,还作为独立交付物写到 <report_dir>/figures/
# <name>.<fmt>(每种请求格式一个文件);manifest.tsv 逐行列出
# name<TAB>filename,让 DAG 在不逐一声明条件图的情况下保持诚实。
FIG_DIR = ""  # set in main(): <img_dir>
_MANIFEST_PATH = ""  # set in main(): <img_dir>/manifest.tsv
_FIG_MANIFEST: List[Tuple[str, str]] = []


def save_fig(fig, name: str, image_formats: Optional[List[PlotFormat]] = None,
             dpi: int = 200) -> str:
    """Save one report figure as ``figures/<name>.<fmt>`` for each requested
    format (standalone deliverable) and return the raster copy embedded in
    the PPTX (the png variant when produced, else a tempfile).

    image_formats: Save formats. Default ``["png"]``.
    """
    if image_formats is None:
        image_formats = ["png"]
    bad = [f for f in image_formats if f not in VALID_FORMATS]
    if bad:
        raise ValueError(
            f"Unsupported image format(s): {', '.join(bad)}; "
            f"valid: {', '.join(VALID_FORMATS)}")
    os.makedirs(FIG_DIR, exist_ok=True)
    for fmt in image_formats:
        p = os.path.join(FIG_DIR, f"{name}.{fmt}")
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        _FIG_MANIFEST.append((name, f"{name}.{fmt}"))
    if "png" in image_formats:
        return os.path.join(FIG_DIR, f"{name}.png")
    fd, tmp = tempfile.mkstemp(prefix=f"{name}_", suffix=".png")
    os.close(fd)
    fig.savefig(tmp, dpi=dpi, bbox_inches="tight")
    return tmp


def write_fig_manifest() -> str:
    """Write figures/manifest.tsv (one ``name<TAB>filename`` line per
    exported figure file). Always written so the DAG stays honest even
    when no figure was rendered."""
    os.makedirs(FIG_DIR, exist_ok=True)
    with open(_MANIFEST_PATH, "w", encoding="utf-8") as fh:
        for name, filename in _FIG_MANIFEST:
            fh.write(f"{name}\t{filename}\n")
    return _MANIFEST_PATH


# ============================================================
# Data collection (DNA chain; collectors copied from Mutation_report)
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


def read_table(path: str, sep: str = "\t") -> Optional[pd.DataFrame]:
    """Read a TSV/CSV into a DataFrame; None when absent or unparseable."""
    if not path or not os.path.isfile(path):
        return None
    try:
        return pd.read_csv(path, sep=sep)
    except Exception:
        return None


def load_burden_summary(burden_dir: str) -> Optional[pd.DataFrame]:
    """Load the burden module's cross-sample burden_summary.tsv."""
    return read_table(os.path.join(burden_dir, "burden_summary.tsv"), sep="\t")


def df_to_rows(df: pd.DataFrame) -> List[List]:
    """DataFrame -> string-safe rows (NA for missing); header handled apart."""
    return [["NA" if pd.isna(v) else v for v in rec]
            for rec in df.astype(object).where(df.notna(), None).values.tolist()]


def df_to_table_rows(df: pd.DataFrame, max_cols: int = 8,
                     max_rows: int = 20) -> List[List[str]]:
    """DataFrame -> [header, *rows] string table (capped, NA for missing)."""
    cols = list(df.columns[:max_cols])
    rows: List[List[str]] = [[str(c) for c in cols]]
    for _, rec in df.head(max_rows).iterrows():
        cells = []
        for c in cols:
            v = rec[c]
            cells.append("NA" if pd.isna(v) else str(v))
        rows.append(cells)
    return rows


def collect_qc_rows(qc_path: str) -> Tuple[List[List[str]], Dict[str, List[float]]]:
    """QC table rows + numeric columns for coloring from qc_summary.tsv.

    Known informative columns are preferred when present; otherwise the
    first columns are shown so chain-specific qc_summary layouts render.
    """
    if not os.path.isfile(qc_path):
        return [], {}
    try:
        df = pd.read_csv(qc_path, sep="\t")
    except Exception:
        return [], {}
    if len(df) == 0:
        return [], {}
    if "sample_id" in df.columns:
        df = df.set_index("sample_id")
    known = ["total_sequences", "total", "mapped_rate", "percent_duplication", "msi_score"]
    show = [c for c in known if c in df.columns]
    if not show:
        show = list(df.columns[:7])
    rows = [[_hdr(c) for c in ["sample_id"] + show]]
    for sample, rec in df.iterrows():
        cells = [str(sample)]
        for c in show:
            missing = bool(pd.isna(rec[c]))
            cells.append("NA" if missing else str(rec[c]))
        rows.append(cells)
    numeric: Dict[str, List[float]] = {}
    for col in ("mapped_rate", "percent_duplication"):
        if col in df.columns:
            series = df[col]
            numeric[col] = []
            for v in list(series.astype(object).where(series.notna(), None)):
                try:
                    numeric[col].append(float(v))
                except (TypeError, ValueError):
                    numeric[col].append(float("nan"))
    return rows, numeric


def _hdr(name: str) -> str:
    """Column header display name (sample_id -> Sample)."""
    return "sample" if name == "sample_id" else str(name)


def collect_sexcheck_rows(path: str, max_cols: int = 8,
                          max_rows: int = 10) -> List[List[str]]:
    """sexcheck.tsv as a string table (generic columns)."""
    df = read_table(path, sep="\t")
    if df is None or df.empty:
        return []
    return df_to_table_rows(df, max_cols=max_cols, max_rows=max_rows)


# ============================================================
# DR-unique data collection (burden / DRclass / TransEffect)
# ============================================================

DRCLASS_LABELS = ["D1R1", "D1R2", "D1R3", "D1R4", "D1R6", "D1R5",
                  "D2R1", "D3R1", "D6R1", "D5R1", "D0R1", "R1-ANCHOR-REJ"]


def load_drclass_summary(drclass_dir: str) -> Optional[pd.DataFrame]:
    """DRclass_summary.tsv (per-sample label counts + penetrance + QC)."""
    return read_table(os.path.join(drclass_dir, "DRclass_summary.tsv"), sep="\t")


def load_teffect_tables(teffect_dir: str) -> Dict[str, pd.DataFrame]:
    """All TransEffect tables present under teffect_dir, keyed by stem."""
    out: Dict[str, pd.DataFrame] = {}
    for stem in ("rna_support_sites", "vaf_ratios", "variant_gene_delta",
                 "cnv_dosage_summary", "cnv_segments", "program_overlap"):
        df = read_table(os.path.join(teffect_dir, f"{stem}.tsv"), sep="\t")
        if df is not None and not df.empty:
            out[stem] = df
    return out


# ============================================================
# RNA chain data collection
# ============================================================

def count_lines(path: str) -> int:
    """Count non-empty lines of a text file (0 when absent)."""
    if not path or not os.path.isfile(path):
        return 0
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n


def collect_de_summary(rna_dir: str, genome_label: str,
                       contrasts: List[str]) -> List[List[str]]:
    """Per-contrast up/down gene counts from the function-stage gene lists.

    Rows: [contrast, up, down]; both counts NA when the lists are absent.
    """
    rows: List[List[str]] = []
    for contrast in contrasts:
        fdir = os.path.join(rna_dir, "function", genome_label, contrast)
        up = count_lines(os.path.join(fdir, "up_genes.txt"))
        down = count_lines(os.path.join(fdir, "down_genes.txt"))
        rows.append([contrast,
                     str(up) if up else ("0" if os.path.isdir(fdir) else "NA"),
                     str(down) if down else ("0" if os.path.isdir(fdir) else "NA")])
    return rows


def collect_func_top_terms(func_dir: str, csv_name: str, top: int = 3) -> List[str]:
    """Top term names of one enrichment CSV (first column), for bullets."""
    df = read_table(os.path.join(func_dir, csv_name), sep=",")
    if df is None or df.empty or df.shape[1] == 0:
        return []
    return [str(v) for v in df.iloc[:top, 0].tolist() if not pd.isna(v)]


# ============================================================
# Report-own charts
# ============================================================

def plot_somatic_counts(rows: List[List[str]],
                        image_formats: Optional[List[PlotFormat]] = None,
                        dpi: int = 200) -> Optional[str]:
    """Bar chart of filtered somatic calls per tumor/normal pair."""
    if not image_formats:
        return None
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
    raster = save_fig(fig, "somatic_calls_per_pair", image_formats, dpi=dpi)
    plt.close(fig)
    return raster


def plot_burden_overview(df: pd.DataFrame,
                         image_formats: Optional[List[PlotFormat]] = None,
                         dpi: int = 200) -> Optional[str]:
    """Three-panel per-sample bars: SNV/InDel, SV, CNV altered segments."""
    if df is None or df.empty or not image_formats:
        return None
    panels = [("snv_count", "SNV count"), ("sv_total", "SV total"),
              ("cnv_total_segments", "CNV segments")]
    present = [(c, t) for c, t in panels if c in df.columns]
    if not present:
        return None
    fig, axes = plt.subplots(1, len(present),
                             figsize=(3.2 * len(present), 3.0),
                             squeeze=False)
    for ax, (col, title) in zip(axes[0], present):
        vals = pd.to_numeric(df[col], errors="coerce").fillna(0)
        ax.bar(range(len(vals)), vals.tolist(), color="#0094D8")
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels([str(s) for s in df.get("sample", df.index)],
                           rotation=45, ha="right", fontsize=7)
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    raster = save_fig(fig, "burden_overview", image_formats, dpi=dpi)
    plt.close(fig)
    return raster


def plot_drclass_labels(df: pd.DataFrame,
                        image_formats: Optional[List[PlotFormat]] = None,
                        dpi: int = 200) -> Optional[str]:
    """Stacked bar of D×R label counts per DRclass sample."""
    if df is None or df.empty or not image_formats:
        return None
    labels = [c for c in DRCLASS_LABELS if c in df.columns]
    if not labels:
        return None
    samples = [str(s) for s in df.iloc[:, 0].tolist()]
    bottom = [0.0] * len(samples)
    fig, ax = plt.subplots(figsize=(7, max(2.5, 0.6 * len(samples) + 1.5)))
    cmap = plt.get_cmap("tab20")
    for i, label in enumerate(labels):
        vals = pd.to_numeric(df[label], errors="coerce").fillna(0).tolist()
        ax.bar(range(len(samples)), vals, bottom=bottom, label=label,
               color=cmap(i % 20))
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Loci")
    ax.set_title("DRclass DxR label composition per sample")
    ax.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    raster = save_fig(fig, "drclass_label_composition", image_formats, dpi=dpi)
    plt.close(fig)
    return raster


def plot_de_updown(rows: List[List[str]],
                   image_formats: Optional[List[PlotFormat]] = None,
                   dpi: int = 200) -> Optional[str]:
    """Grouped bars of up/down gene counts per contrast."""
    if not image_formats:
        return None
    data = [(r[0], r[1], r[2]) for r in rows
            if r[1] not in ("NA", "") or r[2] not in ("NA", "")]
    if not data:
        return None
    labels = [d[0] for d in data]
    ups = [float(d[1]) if d[1] not in ("NA", "") else 0.0 for d in data]
    downs = [float(d[2]) if d[2] not in ("NA", "") else 0.0 for d in data]
    x = range(len(labels))
    fig, ax = plt.subplots(figsize=(7, max(2.5, 0.4 * len(labels) + 1.5)))
    ax.bar([i - 0.2 for i in x], ups, width=0.4, label="up", color="#D64545")
    ax.bar([i + 0.2 for i in x], downs, width=0.4, label="down", color="#0094D8")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Genes")
    ax.set_title("DE genes per contrast")
    ax.legend(fontsize=8)
    fig.tight_layout()
    raster = save_fig(fig, "de_updown_per_contrast", image_formats, dpi=dpi)
    plt.close(fig)
    return raster


# ============================================================
# Excel workbook
# ============================================================

_MAX_ROWS = 5000


def _safe_sheet_name(name: str) -> str:
    for ch in "[]:*?/\\ ":
        name = name.replace(ch, "_")
    return name[:31]


class _SheetWriter:
    """Multi-sheet workbook writer with unique sheet names and NaN-safe
    column widths (pandas 3.x keeps NaN through astype(str), so an
    all-NaN column must not reach int())."""

    def __init__(self, output_path: str):
        import openpyxl
        self._openpyxl = openpyxl
        self.wb = openpyxl.Workbook()
        self.wb.remove(self.wb.active)
        self._used: set = set()
        self.output_path = output_path

    def unique_name(self, name: str) -> str:
        base = _safe_sheet_name(name)
        if base not in self._used:
            self._used.add(base)
            return base
        for i in range(1, 100):
            candidate = base[: 31 - len(str(i))] + str(i)
            if candidate not in self._used:
                self._used.add(candidate)
                return candidate
        return base

    def write_df(self, name: str, df: Optional[pd.DataFrame],
                 first_col_note: str = "") -> None:
        ws = self.wb.create_sheet(self.unique_name(name))
        if first_col_note:
            ws.append([first_col_note])
        if df is None or df.empty:
            ws.append(["(no data / 未运行)"])
            return
        if len(df) > _MAX_ROWS:
            df = df.head(_MAX_ROWS)
        ws.append([str(c) for c in df.columns])
        for _, row in df.iterrows():
            ws.append(["" if pd.isna(row[c]) else str(row[c]) for c in df.columns])
        self._auto_width(ws, df)

    def write_rows(self, name: str, rows: List[List[str]]) -> None:
        ws = self.wb.create_sheet(self.unique_name(name))
        if not rows:
            ws.append(["(no data / 未运行)"])
            return
        for row in rows[: _MAX_ROWS + 1]:
            ws.append([str(v) for v in row])
        self._auto_width(ws)

    def _auto_width(self, ws, df: Optional[pd.DataFrame] = None) -> None:
        openpyxl = self._openpyxl
        for i, col_cells in enumerate(ws.iter_cols(min_row=1, max_row=min(ws.max_row, 200)), 1):
            max_len = 0
            for cell in col_cells:
                val = "" if cell.value is None else str(cell.value)
                max_len = max(max_len, len(val))
            ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = min(max_len + 2, 50)

    def save(self) -> None:
        self.wb.save(self.output_path)


# ============================================================
# Slides
# ============================================================

def _qc_sexcheck_slide(slide, root: str, lang: str) -> None:
    """DNA/RNA QC slide: qc_summary table + sexcheck table, hard-clipped to
    the slide bottom margin (report-module-creation table sizing rule:
    y + N * row_h must stay below 5.405")."""
    qc_rows, _numeric = collect_qc_rows(os.path.join(root, "QC", "qc_summary", "qc_summary.tsv"))
    sex_rows = collect_sexcheck_rows(os.path.join(root, "QC", "sexcheck", "sexcheck.tsv"))
    y = CONTENT_TOP + 0.1
    if qc_rows:
        h = min(0.3 * len(qc_rows), 1.9)
        _table(slide, MARGIN_L, y, CONTENT_W, h, qc_rows, font_size=9)
        y += h + 0.25
    if sex_rows:
        _textbox(slide, MARGIN_L, y, CONTENT_W, 0.3, "sexcheck", font_size=11, bold=True)
        y += 0.32
        h = max(0.6, min(0.28 * len(sex_rows), _SLIDE_BOTTOM - y))
        _table(slide, MARGIN_L, y, CONTENT_W, h, sex_rows, font_size=9)
    if not qc_rows and not sex_rows:
        _textbox(slide, MARGIN_L, CONTENT_TOP + 1.5, CONTENT_W, 0.5,
                 _t("not_run", lang), font_size=12, color=C_MUTED, align=PP_ALIGN.CENTER)


def build_report(args: argparse.Namespace) -> None:
    lang = args.lang
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)
    blank = prs.slide_layouts[6]

    def new_slide(title: str):
        s = prs.slides.add_slide(blank)
        _bg(s)
        _header(s, title)
        return s

    pairs: List[Tuple[str, str]] = args.pairs
    genome = args.genome or args.genome_label

    # ---- S1 title ----
    s = prs.slides.add_slide(blank)
    _bg(s)
    bar = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, Inches(1.6), Inches(SLIDE_W), Inches(1.5))
    bar.fill.solid(); bar.fill.fore_color.rgb = C_NAVY; bar.line.fill.background()
    _textbox(s, 0.6, 1.85, SLIDE_W - 1.2, 0.7, args.title or _t("title_default", lang),
             font_size=28, bold=True, color=C_WHITE, align=PP_ALIGN.CENTER)
    sub = " / ".join(x for x in (args.subtitle, genome, args.date) if x)
    _textbox(s, 0.6, 2.55, SLIDE_W - 1.2, 0.4, sub, font_size=14, color=C_WHITE,
             align=PP_ALIGN.CENTER)
    if args.pipeline:
        _textbox(s, 0.8, 3.4, SLIDE_W - 1.6, 1.2, args.pipeline, font_size=11,
                 color=C_MUTED, align=PP_ALIGN.CENTER)

    # ---- S2 overview ----
    s = new_slide(_t("overview", lang))
    cards = [
        (_t("chain_dna", lang),
         f"{len(args.dna_samples) or len(args.samples)} samples; "
         f"{len(pairs)} somatic pair(s)"),
        (_t("chain_rna", lang),
         f"{len(args.rna_samples) or len(args.samples)} samples; "
         f"{len(args.contrasts)} contrast(s)"),
        (_t("chain_dr", lang), "DRclass / TransEffect / burden"),
    ]
    for i, (name, desc) in enumerate(cards):
        x = MARGIN_L + i * (CONTENT_W / 3)
        card = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x),
                                  Inches(CONTENT_TOP + 0.3),
                                  Inches(CONTENT_W / 3 - 0.2), Inches(1.8))
        card.fill.solid(); card.fill.fore_color.rgb = C_WHITE
        card.line.color.rgb = C_ACCENT
        _textbox(s, x + 0.15, CONTENT_TOP + 0.5, CONTENT_W / 3 - 0.5, 0.6, name,
                 font_size=13, bold=True, color=C_NAVY)
        _textbox(s, x + 0.15, CONTENT_TOP + 1.15, CONTENT_W / 3 - 0.5, 0.8, desc,
                 font_size=10, color=C_MUTED)
    if args.pipeline:
        _textbox(s, MARGIN_L, CONTENT_TOP + 2.4, CONTENT_W, 1.6, args.pipeline,
                 font_size=10, color=C_TEXT)

    # ---- S3 DNA QC ----
    s = new_slide(f"{_t('dna_qc', lang)} (DNA)")
    _qc_sexcheck_slide(s, args.dna_dir, lang)

    # ---- S4 DNA variants (somatic + MSI) ----
    s = new_slide(_t("dna_variants", lang))
    somatic_rows = collect_somatic(os.path.join(args.dna_dir, "variation", "somatic_snv_indel"), pairs)
    if somatic_rows:
        header = [["pair", "tumor", "normal", "filtered_snv_indel"]] + somatic_rows
        _table(s, MARGIN_L, CONTENT_TOP + 0.1, CONTENT_W * 0.55,
               min(2.2, 0.3 * len(header)), header, font_size=9)
        img = plot_somatic_counts(somatic_rows, args.formats)
        if img:
            _add_picture(s, img, MARGIN_L + CONTENT_W * 0.58, CONTENT_TOP + 0.1,
                         CONTENT_W * 0.4, 2.4)
    msi_rows = collect_msi(os.path.join(args.dna_dir, "stability", "msi"))
    if msi_rows:
        top = CONTENT_TOP + 2.5
        _textbox(s, MARGIN_L, top, CONTENT_W, 0.3, "MSI score", font_size=11, bold=True)
        _table(s, MARGIN_L, top + 0.3, CONTENT_W * 0.5,
               min(1.8, 0.3 * (len(msi_rows) + 1)),
               [["sample", "msi_score"]] + msi_rows, font_size=9)
    if not somatic_rows and not msi_rows:
        _textbox(s, MARGIN_L, CONTENT_TOP + 1.5, CONTENT_W, 0.5,
                 _t("not_run", lang), font_size=12, color=C_MUTED, align=PP_ALIGN.CENTER)

    # ---- S5 spectrum / signatures ----
    s = new_slide(_t("spectrum", lang))
    spec_dir = os.path.join(args.dna_dir, "stability", "spectrum")
    figs = [_pick_variant(os.path.join(spec_dir, n)) for n in
            ("somatic_spectrum_stacked_bar.png", "somatic_spectrum_heatmap.png",
             "signature_refit_exposure_stacked_bar.png", "signature_refit_exposure_heatmap.png")]
    shown = 0
    for i, p in enumerate(figs):
        if _embeddable(p):
            _add_picture(s, p, MARGIN_L + (i % 2) * (CONTENT_W / 2),
                         CONTENT_TOP + 0.1 + (i // 2) * 2.15, CONTENT_W / 2 - 0.15, 2.0)
            shown += 1
    if not shown:
        _picture_or_note(s, None, _t("spectrum", lang), lang)

    # ---- S6 fragment size ----
    s = new_slide(_t("fragment", lang))
    frag = _pick_variant(os.path.join(args.dna_dir, "results", "fragment_size",
                                      "fragment", "FragmentSize.png"))
    _picture_or_note(s, frag, _t("fragment", lang), lang)

    # ---- S7 burden (DR-unique) ----
    s = new_slide(f"★ {_t('burden', lang)}")
    bdf = load_burden_summary(args.burden_dir)
    if bdf is not None and not bdf.empty:
        _table(s, MARGIN_L, CONTENT_TOP + 0.1, CONTENT_W * 0.98,
               min(2.3, 0.28 * (min(len(bdf), 12) + 1)), df_to_table_rows(bdf, max_cols=10),
               font_size=8)
        img = plot_burden_overview(bdf, args.formats)
        if img:
            _add_picture(s, img, MARGIN_L + CONTENT_W * 0.1,
                         CONTENT_TOP + 2.5, CONTENT_W * 0.8, 1.9)
    else:
        _picture_or_note(s, None, _t("burden", lang), lang)

    # ---- S8 RNA QC ----
    s = new_slide(f"{_t('rna_qc', lang)} (RNA)")
    _qc_sexcheck_slide(s, args.rna_dir, lang)

    # ---- S9 TE chimeric ----
    s = new_slide(_t("te_chimeric", lang))
    te_dir = os.path.join(args.rna_dir, "transcripts", args.genome_label, "TE_chimeric")
    te_df = read_table(os.path.join(te_dir, "TE_chimeric_sample_summary.tsv"), sep="\t")
    if te_df is not None and not te_df.empty:
        _table(s, MARGIN_L, CONTENT_TOP + 0.1, CONTENT_W * 0.6,
               min(2.4, 0.28 * (min(len(te_df), 12) + 1)),
               df_to_table_rows(te_df, max_cols=8), font_size=8)
        img = _pick_variant(os.path.join(te_dir, "TE_chimeric_te_type_top.png"))
        if _embeddable(img):
            _add_picture(s, img, MARGIN_L + CONTENT_W * 0.63, CONTENT_TOP + 0.1,
                         CONTENT_W * 0.35, 2.4)
    else:
        _picture_or_note(s, None, _t("te_chimeric", lang), lang)

    # ---- S10 DE ----
    s = new_slide(_t("de", lang))
    de_rows = collect_de_summary(args.rna_dir, args.genome_label, args.contrasts)
    if de_rows:
        _table(s, MARGIN_L, CONTENT_TOP + 0.1, CONTENT_W * 0.5,
               min(2.4, 0.28 * (len(de_rows) + 1)),
               [[_t("contrast", lang), _t("up", lang), _t("down", lang)]] + de_rows,
               font_size=9)
        img = plot_de_updown(de_rows, args.formats)
        if img:
            _add_picture(s, img, MARGIN_L + CONTENT_W * 0.52, CONTENT_TOP + 0.1,
                         CONTENT_W * 0.46, 2.6)
    else:
        _picture_or_note(s, None, _t("de", lang), lang)

    # ---- S11 function ----
    s = new_slide(_t("func", lang))
    items: List[str] = []
    for contrast in args.contrasts:
        fdir = os.path.join(args.rna_dir, "function", args.genome_label, contrast)
        top_up = collect_func_top_terms(fdir, "go_up.csv")
        top_kegg = collect_func_top_terms(fdir, "kegg_up.csv")
        line = f"{contrast}: "
        line += "GO↑ " + ", ".join(top_up) if top_up else ""
        line += (" | " if top_up and top_kegg else "") + ("KEGG↑ " + ", ".join(top_kegg) if top_kegg else "")
        if line != f"{contrast}: ":
            items.append(line)
        if len(items) >= 8:
            break
    if items:
        _bullets(s, MARGIN_L, CONTENT_TOP + 0.1, CONTENT_W, 2.6, items, font_size=10)
        img = _pick_variant(os.path.join(args.rna_dir, "function", args.genome_label,
                                         "gsva", "gsva_heatmap.png"))
        if _embeddable(img):
            _add_picture(s, img, MARGIN_L + CONTENT_W * 0.25, CONTENT_TOP + 2.8,
                         CONTENT_W * 0.5, 1.5)
    else:
        _picture_or_note(s, None, _t("func", lang), lang)

    # ---- S12 fusion ----
    s = new_slide(_t("fusion", lang))
    fu_dir = os.path.join(args.rna_dir, "fusion", args.genome_label, "arriba_report")
    fu_df = read_table(os.path.join(fu_dir, "per_sample_summary.tsv"), sep="\t")
    if fu_df is not None and not fu_df.empty:
        _table(s, MARGIN_L, CONTENT_TOP + 0.1, CONTENT_W * 0.9,
               min(2.3, 0.28 * (min(len(fu_df), 12) + 1)),
               df_to_table_rows(fu_df, max_cols=9), font_size=8)
        img = _pick_variant(os.path.join(fu_dir, "figures", "fig1_per_sample_counts.png"))
        if _embeddable(img):
            _add_picture(s, img, MARGIN_L + CONTENT_W * 0.2, CONTENT_TOP + 2.5,
                         CONTENT_W * 0.6, 1.9)
    else:
        _picture_or_note(s, None, _t("fusion", lang), lang)

    # ---- S13 DRclass (DR-unique) ----
    s = new_slide(f"★ {_t('drclass', lang)}")
    ddf = load_drclass_summary(args.drclass_dir)
    if ddf is not None and not ddf.empty:
        _table(s, MARGIN_L, CONTENT_TOP + 0.1, CONTENT_W * 0.98,
               min(2.2, 0.26 * (min(len(ddf), 10) + 1)),
               df_to_table_rows(ddf, max_cols=10), font_size=7)
        img = plot_drclass_labels(ddf, args.formats)
        if img:
            _add_picture(s, img, MARGIN_L + CONTENT_W * 0.12, CONTENT_TOP + 2.4,
                         CONTENT_W * 0.76, 2.0)
    else:
        _picture_or_note(s, None, _t("drclass", lang), lang)

    # ---- S14 TransEffect (DR-unique) ----
    s = new_slide(f"★ {_t('teffect', lang)}")
    teffect = load_teffect_tables(args.teffect_dir)
    if teffect:
        vaf = teffect.get("vaf_ratios")
        if vaf is not None:
            _table(s, MARGIN_L, CONTENT_TOP + 0.1, CONTENT_W * 0.55,
                   min(2.2, 0.28 * (min(len(vaf), 8) + 1)),
                   df_to_table_rows(vaf, max_cols=6), font_size=8)
        fig_paths = [_pick_variant(os.path.join(args.teffect_dir, n)) for n in
                     ("fig_rna_support_vaf.png", "fig_variant_gene_delta.png",
                      "fig_cnv_overview.png", "fig_program_overlap.png")]
        shown = 0
        for p in fig_paths:
            if _embeddable(p):
                _add_picture(s, p, MARGIN_L + CONTENT_W * 0.58,
                             CONTENT_TOP + 0.1 + shown * 1.15, CONTENT_W * 0.4, 1.1)
                shown += 1
        if not shown and vaf is None:
            _picture_or_note(s, None, _t("teffect", lang), lang)
    else:
        _picture_or_note(s, None, _t("teffect", lang), lang)

    # ---- S15 conclusion ----
    s = new_slide(_t("conclusion", lang))
    n_pairs = sum(1 for r in collect_somatic(
        os.path.join(args.dna_dir, "variation", "somatic_snv_indel"), pairs) if r[3] != "NA")
    items = [
        f"DNA: {len(args.dna_samples) or len(args.samples)} samples, "
        f"{n_pairs} somatic pair(s) with calls",
        f"RNA: {len(args.rna_samples) or len(args.samples)} samples, "
        f"{len(args.contrasts)} contrast(s)",
    ]
    bdf = load_burden_summary(args.burden_dir)
    items.append(f"burden: {len(bdf)} sample(s) aggregated" if bdf is not None and not bdf.empty
                 else "burden: not run")
    ddf = load_drclass_summary(args.drclass_dir)
    items.append(f"DRclass: {len(ddf)} sample(s) classified" if ddf is not None and not ddf.empty
                 else "DRclass: not run")
    items.append(f"TransEffect: {len(load_teffect_tables(args.teffect_dir))} table(s)"
                 if os.path.isdir(args.teffect_dir) else "TransEffect: not run")
    _bullets(s, MARGIN_L, CONTENT_TOP + 0.3, CONTENT_W, 3.2, items, font_size=13)

    prs.save(args.output)


# ============================================================
# Excel inventory
# ============================================================

def write_results_excel(output_path: str, args: argparse.Namespace) -> None:
    """Write the results workbook: one sheet per analysis output plus a
    file inventory. Missing stages leave a '(no data)' sheet — the report
    must not fail because an optional stage was not registered."""
    w = _SheetWriter(output_path)
    gl = args.genome_label

    # DNA chain
    w.write_df("dna_qc", read_table(os.path.join(args.dna_dir, "QC", "qc_summary", "qc_summary.tsv")))
    w.write_df("dna_sexcheck", read_table(os.path.join(args.dna_dir, "QC", "sexcheck", "sexcheck.tsv")))
    w.write_rows("somatic_calls",
                 [["pair", "tumor", "normal", "filtered_snv_indel"]] + collect_somatic(
                     os.path.join(args.dna_dir, "variation", "somatic_snv_indel"), args.pairs))
    w.write_rows("msi_scores", [["sample", "msi_score"]] + collect_msi(
        os.path.join(args.dna_dir, "stability", "msi")))
    spec_dir = os.path.join(args.dna_dir, "stability", "spectrum")
    w.write_df("spectrum_sample_stats", read_table(os.path.join(spec_dir, "somatic_spectrum_sample_stats.csv"), ","))
    w.write_df("signature_refit_summary", read_table(os.path.join(spec_dir, "signature_refit_summary.csv"), ","))
    w.write_df("signature_nmf_best", read_table(os.path.join(spec_dir, "signature_nmf_cosmic_best_match.csv"), ","))
    w.write_df("burden_summary", load_burden_summary(args.burden_dir))

    # DR-unique
    w.write_df("drclass_summary", load_drclass_summary(args.drclass_dir))
    if os.path.isdir(args.drclass_dir):
        for sample in sorted(os.listdir(args.drclass_dir)):
            for kind in ("dna_anchored", "rna_anchored"):
                p = os.path.join(args.drclass_dir, sample, f"{sample}.{kind}.tsv")
                if os.path.isfile(p):
                    w.write_df(f"drclass_{sample}_{kind}", read_table(p))
    for stem, df in load_teffect_tables(args.teffect_dir).items():
        w.write_df(f"teffect_{stem}", df)

    # RNA chain
    w.write_df("rna_qc", read_table(os.path.join(args.rna_dir, "QC", "qc_summary", "qc_summary.tsv")))
    w.write_df("rna_sexcheck", read_table(os.path.join(args.rna_dir, "QC", "sexcheck", "sexcheck.tsv")))
    te_dir = os.path.join(args.rna_dir, "transcripts", gl, "TE_chimeric")
    w.write_df("te_chimeric_sample", read_table(os.path.join(te_dir, "TE_chimeric_sample_summary.tsv")))
    w.write_df("te_chimeric_group", read_table(os.path.join(te_dir, "TE_chimeric_group_summary.tsv")))
    tecount = os.path.join(args.rna_dir, "counts", "TEtranscripts", gl, f"{gl}_TEcount.tsv")
    if os.path.isfile(tecount) and os.path.getsize(tecount) < 10 * 1024 * 1024:
        w.write_df("tecount_matrix", read_table(tecount))
    fc = os.path.join(args.rna_dir, "counts", "featureCounts", gl, f"{gl}_featureCounts.tsv")
    if os.path.isfile(fc) and os.path.getsize(fc) < 10 * 1024 * 1024:
        w.write_df("featurecounts_matrix", read_table(fc))
    w.write_df("gsva_scores", read_table(os.path.join(args.rna_dir, "function", gl, "gsva", "gsva_scores.tsv")))
    for contrast in args.contrasts:
        de_dir = os.path.join(args.rna_dir, "diff_expression", gl, contrast, "TEcount")
        w.write_df(f"de_{contrast}_gene_name",
                   read_table(os.path.join(de_dir, f"{contrast}.TEcount_Gene.name.tsv")))
        fdir = os.path.join(args.rna_dir, "function", gl, contrast)
        for csv_name, sheet in (("go_up.csv", "go_up"), ("go_down.csv", "go_down"),
                                ("kegg_up.csv", "kegg_up"), ("kegg_down.csv", "kegg_down"),
                                (os.path.join("GSEA", "TEcount_Gene_GSEA.csv"), "gsea")):
            df = read_table(os.path.join(fdir, csv_name), ",")
            w.write_df(f"func_{contrast}_{sheet}", df,
                       first_col_note=contrast)
    fu_dir = os.path.join(args.rna_dir, "fusion", gl, "arriba_report")
    for name, sheet in (("per_sample_summary.tsv", "fusion_summary"),
                        ("recurrent_fusions.tsv", "recurrent_fusions"),
                        ("high_medium_confidence_fusions.tsv", "high_medium_fusions"),
                        ("inframe_fusions.tsv", "inframe_fusions")):
        w.write_df(sheet, read_table(os.path.join(fu_dir, name)))

    # file inventory
    inv: List[List[str]] = [["path"]]
    for root in (args.dna_dir, args.rna_dir, args.drclass_dir, args.teffect_dir):
        if root and os.path.isdir(root):
            for dirpath, _dirnames, filenames in os.walk(root):
                for fn in sorted(filenames):
                    if fn.endswith((".tsv", ".csv", ".txt")):
                        inv.append([os.path.join(dirpath, fn)])
    w.write_rows("file_inventory", inv)
    w.save()


# ============================================================
# CLI
# ============================================================

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-a", "--analysis-dir", required=True,
                   help="Workflow output root (contains DNA/ RNA/ DRclass/ TransEffect/)")
    p.add_argument("-d", "--dna-dir", default="", help="DNA chain output root")
    p.add_argument("-r", "--rna-dir", default="", help="RNA chain output root")
    p.add_argument("-k", "--drclass-dir", default="", help="DRclass module output dir")
    p.add_argument("-t", "--teffect-dir", default="", help="TransEffect module output dir")
    p.add_argument("-b", "--burden-dir", default="", help="burden module output dir")
    p.add_argument("-o", "--output", required=True, help="Output PPTX path")
    p.add_argument("-e", "--file-inventory", default="", dest="file_inventory",
                   help="Output XLSX path (default: <output stem>_files.xlsx)")
    p.add_argument("-i", "--img-dir", default="", dest="img_dir",
                   help="Report figure export dir (default: <output dir>/figures)")
    p.add_argument("-s", "--sample", action="append", default=[], dest="samples",
                   metavar="ID", help="All sample ids (repeatable)")
    p.add_argument("-D", "--dna-sample", action="append", default=[], dest="dna_samples",
                   metavar="ID", help="DNA sample id (repeatable)")
    p.add_argument("-R", "--rna-sample", action="append", default=[], dest="rna_samples",
                   metavar="ID", help="RNA sample id (repeatable)")
    p.add_argument("-p", "--paired-sample", action="append", default=[],
                   dest="paired_samples", metavar="ID", help="Paired-end sample (repeatable)")
    p.add_argument("-u", "--single-sample", action="append", default=[],
                   dest="single_samples", metavar="ID", help="Single-end sample (repeatable)")
    p.add_argument("-c", "--contrast", action="append", default=[], dest="contrasts",
                   metavar="A_vs_B", help="DE contrast (repeatable)")
    p.add_argument("-P", "--pair", action="append", default=[], dest="pair_specs",
                   metavar="TUMOR:NORMAL",
                   help="Tumor:normal somatic design pair (repeatable)")
    p.add_argument("-T", "--title", default="DNA×RNA 联合分析报告 / DR Joint Analysis Report",
                   help="Report title")
    p.add_argument("-S", "--subtitle", default="", help="Report subtitle")
    p.add_argument("-L", "--pipeline", default="", help="Pipeline name subtitle")
    p.add_argument("-g", "--genome", default="", help="Genome subtitle")
    p.add_argument("-G", "--genome-label", default="", dest="genome_label",
                   help="Genome dir label of the RNA tree (e.g. GRCm39)")
    p.add_argument("-y", "--date", default="", help="Report date subtitle")
    p.add_argument("-l", "--lang", default="zh", choices=["zh", "en"],
                   help="Report language (default: zh)")
    p.add_argument("-m", "--log", default=None,
                   help="Optional log file path (also logs to stdout)")
    p.add_argument("-f", "--format", action="append", dest="formats", metavar="FMT",
                   help="Image format (png/pdf/svg/...). Repeatable. Default: png. "
                        "Own figures are exported as <img-dir>/<name>.<fmt>.")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    global logger
    if args.log:
        os.makedirs(os.path.dirname(os.path.abspath(args.log)), exist_ok=True)
    logger = setup_logger("DR_report", stream=sys.stdout, log_file=args.log)

    # resolve section roots (explicit flag wins, else derive from analysis-dir)
    if not args.dna_dir:
        args.dna_dir = os.path.join(args.analysis_dir, "DNA")
    if not args.rna_dir:
        args.rna_dir = os.path.join(args.analysis_dir, "RNA")
    if not args.drclass_dir:
        args.drclass_dir = os.path.join(args.analysis_dir, "DRclass")
    if not args.teffect_dir:
        args.teffect_dir = os.path.join(args.analysis_dir, "TransEffect")
    if not args.burden_dir:
        args.burden_dir = os.path.join(args.dna_dir, "stability", "burden")
    if not args.genome_label:
        args.genome_label = args.genome or "all"
    if not args.formats:
        args.formats = ["png"]
    if not args.file_inventory:
        args.file_inventory = os.path.splitext(args.output)[0] + "_files.xlsx"

    pairs: List[Tuple[str, str]] = []
    for spec in args.pair_specs:
        if ":" not in spec:
            logger.error(f"Bad --pair '{spec}'; expected TUMOR:NORMAL")
            return 1
        tumor, normal = spec.split(":", 1)
        pairs.append((tumor, normal))
    args.pairs = pairs

    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)
    global FIG_DIR, _MANIFEST_PATH
    FIG_DIR = args.img_dir or os.path.join(out_dir, "figures")
    _MANIFEST_PATH = os.path.join(FIG_DIR, "manifest.tsv")
    os.makedirs(FIG_DIR, exist_ok=True)

    _setup_cjk_font()
    logger.info(f"Building DR report: dna={args.dna_dir} rna={args.rna_dir} "
                f"drclass={args.drclass_dir} teffect={args.teffect_dir} "
                f"burden={args.burden_dir} genome_label={args.genome_label}")
    build_report(args)
    manifest_path = write_fig_manifest()
    logger.info(f"figure manifest written: {manifest_path}")
    write_results_excel(args.file_inventory, args)
    logger.info(f"results workbook written: {args.file_inventory}")
    logger.info(f"DR report written: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
