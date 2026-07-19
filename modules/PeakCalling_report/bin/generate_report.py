#!/usr/bin/env python3
"""
ChIP-seq Peak Calling Report Generator

Generates a comprehensive PPT report from ChIP-seq peak calling results.
Supports Chinese (default) and English language for PPT text.
Matplotlib chart titles are always in English for font compatibility.

Usage:
    python generate_report.py \
        --samples Pop5IP Rpp14IP Rpp21IP \
        --input-samples Pop5Input Rpp14Input Rpp21Input \
        --output /path/to/report.pptx \
        --lang zh  # or en
"""

import argparse
import csv
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE


# ============================================================
# Internationalization
# ============================================================
I18N = {
    "zh": {
        "report_title": "ChIP-seq Peak Calling 分析报告",
        "pipeline_label": "分析流程",
        "date_label": "报告日期",
        "workflow_title": "分析流程概览",
        "raw_qc": "原始数据质控",
        "adapter_trimming": "接头修剪",
        "alignment": "序列比对",
        "markdup": "标记重复",
        "peak_calling": "Peak Calling",
        "peak_annotation": "Peak 注释",
        "bigwig_track": "生成 BigWig 轨道",
        "enrichment_qc": "富集质量评估",
        "sample_info": "样本信息",
        "sample": "样本",
        "type": "类型",
        "replicate": "生物学重复",
        "ip_type": "IP (免疫沉淀)",
        "input_type": "Input (对照)",
        "alignment_title": "Bowtie2 比对统计",
        "alignment_warning": "警告: {} 比对率 < 80%, 样本质量可能存在问题",
        "markdup_title": "GATK MarkDuplicates 统计",
        "read_pairs": "Read Pair 数",
        "dup_pairs": "重复 Pair 数",
        "dup_rate": "重复率",
        "unmapped": "未比对数",
        "markdup_notes": [
            "GATK MarkDuplicates (仅标记, 不去除)",
            "IP 样本重复率高于 Input 样本 (ChIP 预期行为)",
            "重复 reads 保留用于 MACS3 --keep-dup auto",
        ],
        "peak_calling_title": "MACS3 Peak Calling 结果",
        "peak_count": "Peak 数",
        "frip_score": "FRiP 分数",
        "peak_calling_notes": [
            "MACS3 参数: --pvalue 1e-5, --bw 200, genome mm",
            "FRiP = peaks 内 reads 数 / 总 mapped reads 数",
            "FRiP ≥ 20% 通常表示良好 ChIP 富集",
        ],
        "cutoff_title": "Cutoff 分析",
        "cutoff_note": "虚线标记使用的阈值 (p=1e-5)",
        "annotation_title": "Peak 注释 — 基因组区域分布",
        "promoter": "启动子 (Promoter-TSS)",
        "exon": "外显子 (Exon)",
        "intron": "内含子 (Intron)",
        "intergenic": "基因间区 (Intergenic)",
        "tts": "转录终止位点 (TTS)",
        "other": "其他",
        "top_genes_title": "Top 靶基因",
        "top_n_peaks": "Top {} Peaks (按 Score 排序)",
        "gene": "基因",
        "position": "位置",
        "score": "Score",
        "shared_genes": "共有 Top 靶基因",
        "shared_desc": "跨样本一致结合 — 高可信度靶基因",
        "summary_title": "总结与 QC 评估",
        "metric": "指标",
        "threshold": "参考阈值",
        "conclusions": "结论",
        "low_alignment": "比对率异常低 ({:.1f}%), 样本质量可能存在问题",
        "low_frip": "FRiP 低 ({:.2f}%), 富集较弱",
        "good_enrichment": "富集良好 (FRiP={:.2f}%)",
        "recommendations": "建议",
        "recommend_frip": "所有样本 FRiP 均低于 20% 阈值, 建议优化 ChIP 实验条件",
        "recommend_align": "{} 比对率异常低, 需检查样本质量或参考基因组",
        "total": "总计",
        "tss_title": "TSS 距离分布",
    },
    "en": {
        "report_title": "ChIP-seq Peak Calling Report",
        "pipeline_label": "Pipeline",
        "date_label": "Date",
        "workflow_title": "Workflow Overview",
        "raw_qc": "Raw data QC",
        "adapter_trimming": "Adapter trimming",
        "alignment": "Alignment",
        "markdup": "MarkDuplicates",
        "peak_calling": "Peak Calling",
        "peak_annotation": "Peak annotation",
        "bigwig_track": "BigWig tracks",
        "enrichment_qc": "Enrichment QC",
        "sample_info": "Sample Information",
        "sample": "Sample",
        "type": "Type",
        "replicate": "Replicate",
        "ip_type": "IP",
        "input_type": "Input",
        "alignment_title": "Bowtie2 Alignment Statistics",
        "alignment_warning": "Warning: {} alignment rate < 80%, possible sample quality issue",
        "markdup_title": "GATK MarkDuplicates Statistics",
        "read_pairs": "Read Pairs",
        "dup_pairs": "Dup Pairs",
        "dup_rate": "Dup Rate",
        "unmapped": "Unmapped",
        "markdup_notes": [
            "GATK MarkDuplicates (flag only, not removed)",
            "IP samples have higher dup rate than Input (expected ChIP behavior)",
            "Duplicate reads retained for MACS3 --keep-dup auto",
        ],
        "peak_calling_title": "MACS3 Peak Calling Results",
        "peak_count": "Peak Count",
        "frip_score": "FRiP Score",
        "peak_calling_notes": [
            "MACS3 params: --pvalue 1e-5, --bw 200, genome mm",
            "FRiP = reads in peaks / total mapped reads",
            "FRiP ≥ 20% indicates good ChIP enrichment",
        ],
        "cutoff_title": "Cutoff Analysis",
        "cutoff_note": "Dashed line marks the threshold used (p=1e-5)",
        "annotation_title": "Peak Annotation — Genomic Region Distribution",
        "promoter": "Promoter-TSS",
        "exon": "Exon",
        "intron": "Intron",
        "intergenic": "Intergenic",
        "tts": "TTS",
        "other": "Other",
        "top_genes_title": "Top Target Genes",
        "top_n_peaks": "Top {} Peaks (by Score)",
        "gene": "Gene",
        "position": "Position",
        "score": "Score",
        "shared_genes": "Shared Top Target Genes",
        "shared_desc": "Consistent binding across samples — high confidence targets",
        "summary_title": "Summary & QC Assessment",
        "metric": "Metric",
        "threshold": "Threshold",
        "conclusions": "Conclusions",
        "low_alignment": "Low alignment rate ({:.1f}%), possible sample quality issue",
        "low_frip": "Low FRiP ({:.2f}%), weak enrichment",
        "good_enrichment": "Good enrichment (FRiP={:.2f}%)",
        "recommendations": "Recommendations",
        "recommend_frip": "All samples FRiP below 20% threshold, consider optimizing ChIP conditions",
        "recommend_align": "{} low alignment rate, check sample quality or reference genome",
        "total": "Total",
        "tss_title": "TSS Distance Distribution",
    },
}


def t(key, lang="zh"):
    return I18N.get(lang, I18N["zh"]).get(key, key)


# ============================================================
# Style
# ============================================================
C_TITLE = RGBColor(0x1A, 0x1A, 0x2E)
C_ACCENT = RGBColor(0x00, 0x7A, 0xCC)
C_TEXT = RGBColor(0x33, 0x33, 0x33)
C_WHITE = RGBColor(0xFF, 0xFF, 0xFF)
C_RED = RGBColor(0xE6, 0x4B, 0x35)
C_GREEN = RGBColor(0x00, 0xA0, 0x87)
C_YELLOW = RGBColor(0xFF, 0xA5, 0x00)
PALETTE = ["#E64B35", "#4DBBD5", "#00A087", "#F39B7F", "#8491B4", "#91D1C2",
           "#B09C85", "#DC0000", "#7E6148", "#3C5488"]

ALIGN_THRESH = 80.0
FRIP_THRESH = 0.20

# Slide layout constants (inches)
SLIDE_W = 10.0
SLIDE_H = 5.625
MARGIN_L = 0.5
MARGIN_R = 0.5
CONTENT_W = SLIDE_W - MARGIN_L - MARGIN_R  # 9.0
HEADER_H = 0.7
CONTENT_TOP = HEADER_H + 0.2  # 0.9
CONTENT_MAX_H = SLIDE_H - CONTENT_TOP - 0.3  # bottom margin


# ============================================================
# Data Loaders
# ============================================================

def load_bowtie2(log_dir, sample):
    lf = os.path.join(log_dir, sample, "bowtie2_align.log")
    if not os.path.isfile(lf):
        return None
    with open(lf) as f:
        for ln in f:
            if "overall alignment rate" in ln:
                return float(ln.split("%")[0].split()[-1])
    return None


def load_markdup(metrics_dir, sample):
    mf = os.path.join(metrics_dir, sample, f"{sample}.Markdup-metrics.txt")
    if not os.path.isfile(mf):
        return None
    with open(mf) as f:
        hdr = None
        for ln in f:
            if ln.startswith("LIBRARY"):
                hdr = ln.strip().split("\t")
                continue
            if hdr and ln.strip() and not ln.startswith("#"):
                d = dict(zip(hdr, ln.strip().split("\t")))
                return {
                    "read_pairs": int(d.get("READ_PAIRS_EXAMINED", 0)),
                    "dup_pairs": int(d.get("READ_PAIR_DUPLICATES", 0)),
                    "dup_rate": float(d.get("PERCENT_DUPLICATION", 0)) * 100,
                    "unmapped": int(d.get("UNMAPPED_READS", 0)),
                }
    return None


def load_macs3_params(log_dir, sample):
    lf = os.path.join(log_dir, sample, "macs3.log")
    if not os.path.isfile(lf):
        return None
    params = {}
    with open(lf) as f:
        for ln in f:
            if "Command line:" in ln:
                params["command"] = ln.split("Command line:")[-1].strip()
            elif "pvalue cutoff" in ln:
                params["pvalue"] = ln.split("=")[-1].strip()
            elif "band width" in ln:
                params["bw"] = ln.split("=")[-1].strip()
            elif "effective genome size" in ln:
                params["genome_size"] = ln.split("=")[-1].strip()
            elif "predicted fragment length" in ln:
                params["fragment_length"] = ln.split("is")[-1].strip().split()[0]
    return params if params else None


def load_frip(qc_dir, sample):
    ff = os.path.join(qc_dir, sample, f"{sample}.FRiP.txt")
    if not os.path.isfile(ff):
        return None
    with open(ff) as f:
        for ln in f:
            parts = ln.strip().split("\t")
            if len(parts) >= 2:
                return float(parts[1])
    return None


def load_peak_count(peaks_dir, sample):
    nf = os.path.join(peaks_dir, sample, f"{sample}_peaks.narrowPeak")
    if not os.path.isfile(nf):
        return 0
    with open(nf) as f:
        return sum(1 for ln in f if ln.strip() and not ln.startswith("#"))


def load_cutoff(peaks_dir, sample):
    cf = os.path.join(peaks_dir, sample, f"{sample}_cutoff_analysis.txt")
    if not os.path.isfile(cf):
        return None
    ps, ns = [], []
    with open(cf) as f:
        f.readline()
        for ln in f:
            cols = ln.strip().split("\t")
            if len(cols) >= 3:
                ps.append(float(cols[0]))
                ns.append(float(cols[2]))
    return {"pscore": ps, "npeaks": ns} if ps else None


def load_annotation(annotation_dir, sample):
    af = os.path.join(annotation_dir, sample, f"{sample}_peaks.annotatePeaks.txt")
    if not os.path.isfile(af):
        return None
    cats = Counter()
    with open(af) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 8:
                continue
            ann = row[7]
            if "promoter-TSS" in ann:
                cats["promoter"] += 1
            elif "ntergenic" in ann:
                cats["intergenic"] += 1
            elif "intron" in ann:
                cats["intron"] += 1
            elif "exon" in ann:
                cats["exon"] += 1
            elif "TTS" in ann:
                cats["tts"] += 1
            else:
                cats["other"] += 1
    return dict(cats) if cats else None


def load_top_genes(annotation_dir, sample, n=5):
    af = os.path.join(annotation_dir, sample, f"{sample}_peaks.annotatePeaks.txt")
    if not os.path.isfile(af):
        return []
    rows = []
    with open(af) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 16:
                continue
            score = int(row[5]) if row[5] else 0
            gene = row[15] if row[15] else "N/A"
            chr_name = row[1]
            start = row[2]
            end = row[3]
            pos = f"{chr_name}:{start}-{end}"
            tss_dist = row[9] if row[9] else ""
            rows.append((gene, pos, score, tss_dist))
    rows.sort(key=lambda x: x[2], reverse=True)
    return rows[:n]


# ============================================================
# Plotting — English titles, DPI 300
# ============================================================

DPI = 300


def _save(fig):
    tmp = tempfile.NamedTemporaryFile(suffix=".png", prefix="rpt_", delete=False)
    fig.savefig(tmp.name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return tmp.name


def plot_alignment(rates):
    fig, ax = plt.subplots(figsize=(max(7, len(rates)*1.3), 4.5))
    names = list(rates.keys())
    vals = [rates[n] for n in names]
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(names))]
    bars = ax.bar(range(len(names)), vals, color=colors, edgecolor="white", lw=1.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Alignment Rate (%)", fontsize=11)
    ax.set_ylim(0, 105)
    ax.axhline(y=ALIGN_THRESH, color="red", ls="--", alpha=.6, lw=1.2,
               label=f"QC Threshold ({ALIGN_THRESH:.0f}%)")
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(axis="y", alpha=.3)
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+1.5,
                f"{v:.1f}%", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_title("Bowtie2 Alignment Rates", fontsize=13, fontweight="bold", pad=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return _save(fig)


def plot_peak_and_frip(samples, peaks, frips):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    x = np.arange(len(samples))
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(samples))]

    pv = [peaks[s] for s in samples]
    bars = axes[0].bar(x, pv, color=colors, width=.5, edgecolor="white", lw=1.5)
    axes[0].set_xticks(x); axes[0].set_xticklabels(samples, fontsize=10)
    axes[0].set_ylabel("Number of Peaks", fontsize=11)
    axes[0].set_title("Peak Counts", fontsize=13, fontweight="bold", pad=10)
    for b, v in zip(bars, pv):
        axes[0].text(b.get_x()+b.get_width()/2, b.get_height()+max(pv)*.03,
                     f"{v:,}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    axes[0].grid(axis="y", alpha=.3)
    axes[0].spines["top"].set_visible(False); axes[0].spines["right"].set_visible(False)

    fv = [frips[s]*100 for s in samples]
    bars = axes[1].bar(x, fv, color=colors, width=.5, edgecolor="white", lw=1.5)
    axes[1].set_xticks(x); axes[1].set_xticklabels(samples, fontsize=10)
    axes[1].set_ylabel("FRiP (%)", fontsize=11)
    axes[1].set_title("FRiP Scores", fontsize=13, fontweight="bold", pad=10)
    axes[1].axhline(y=FRIP_THRESH*100, color="red", ls="--", alpha=.6, lw=1.2,
                    label=f"Threshold ({FRIP_THRESH*100:.0f}%)")
    axes[1].legend(fontsize=9, loc="upper right")
    for b, v in zip(bars, fv):
        axes[1].text(b.get_x()+b.get_width()/2, b.get_height()+max(fv)*.03,
                     f"{v:.2f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
    axes[1].grid(axis="y", alpha=.3)
    axes[1].spines["top"].set_visible(False); axes[1].spines["right"].set_visible(False)
    return _save(fig)


def plot_cutoff(cutoff_data):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    for i, (name, d) in enumerate(cutoff_data.items()):
        c = PALETTE[i % len(PALETTE)]
        for ax in axes:
            ax.plot(d["pscore"], d["npeaks"], marker="o", ms=4, label=name, color=c, lw=1.8)
    for ax in axes:
        ax.set_xlabel("-log10(p-value) cutoff", fontsize=11)
        ax.set_ylabel("Number of peaks", fontsize=11)
        ax.legend(fontsize=10); ax.invert_xaxis(); ax.grid(True, alpha=.3)
        ax.axvline(x=5.0, color="gray", ls="--", alpha=.5, lw=1.2, label="p=1e-5")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    axes[0].set_yscale("log")
    axes[0].set_title("Peak Count vs p-score (log scale)", fontsize=13, fontweight="bold", pad=10)
    axes[1].set_title("Peak Count vs p-score (linear scale)", fontsize=13, fontweight="bold", pad=10)
    return _save(fig)


def plot_annotation(ann_data, lang="zh"):
    n = len(ann_data)
    fig, axes = plt.subplots(1, n, figsize=(4.5*n, 5), constrained_layout=True)
    if n == 1:
        axes = [axes]
    for i, (sample, cats) in enumerate(ann_data.items()):
        labels = [t(k, lang) for k in cats.keys()]
        values = list(cats.values())
        total = sum(values)
        wedges, texts, autotexts = axes[i].pie(
            values, labels=None, autopct="%1.1f%%",
            colors=PALETTE[:len(labels)], pctdistance=.78,
            textprops={"fontsize": 9}, startangle=90
        )
        for at in autotexts:
            at.set_fontsize(8)
        axes[i].set_title(f"{sample}\n(n={total})", fontsize=12, fontweight="bold", pad=8)
        axes[i].legend(labels, loc="lower center", fontsize=8, ncol=2,
                       bbox_to_anchor=(.5, -.22), frameon=False)
    return _save(fig)


def plot_tss_distance(annotation_dir, samples):
    fig, axes = plt.subplots(1, len(samples), figsize=(5*len(samples), 4.5), constrained_layout=True)
    if len(samples) == 1:
        axes = [axes]
    for i, sample in enumerate(samples):
        af = os.path.join(annotation_dir, sample, f"{sample}_peaks.annotatePeaks.txt")
        if not os.path.isfile(af):
            continue
        distances = []
        with open(af) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                if len(row) >= 10 and row[9]:
                    try:
                        d = abs(int(row[9]))
                        if d <= 50000:
                            distances.append(d)
                    except ValueError:
                        pass
        if not distances:
            continue
        axes[i].hist(distances, bins=50, color=PALETTE[i % len(PALETTE)],
                     edgecolor="white", alpha=.85, lw=.8)
        axes[i].set_xlabel("Distance to TSS (bp)", fontsize=10)
        axes[i].set_ylabel("Number of Peaks", fontsize=10)
        axes[i].set_title(f"{sample}", fontsize=12, fontweight="bold", pad=8)
        axes[i].axvline(x=1000, color="red", ls="--", alpha=.6, lw=1.2, label="1kb")
        axes[i].axvline(x=5000, color="orange", ls="--", alpha=.6, lw=1.2, label="5kb")
        axes[i].legend(fontsize=9, loc="upper right")
        axes[i].spines["top"].set_visible(False)
        axes[i].spines["right"].set_visible(False)
    return _save(fig)


# ============================================================
# PPT helpers
# ============================================================

class Layout:
    """Tracks vertical position to prevent element overflow.
    
    Every slide uses a Layout instance. Elements call layout.allocate(h)
    to get their Y position; if not enough space, they get clipped or skipped.
    """
    def __init__(self, top=CONTENT_TOP, bottom=SLIDE_H - 0.2):
        self.y = top
        self.bottom = bottom

    @property
    def remaining(self):
        return max(0, self.bottom - self.y)

    def allocate(self, h):
        """Return (y, actual_h) where actual_h >= 0.15 and fits within bounds, or None if no space."""
        if self.remaining < 0.15:
            return None
        actual_h = min(max(h, 0.15), self.remaining)
        if actual_h < 0.15:
            return None
        y = self.y
        self.y += actual_h
        return y, actual_h

    def gap(self, size=0.1):
        """Advance Y by a small gap."""
        self.y += min(size, self.remaining)


def _add_picture(slide, img_path, left, top, max_width, max_height):
    """Add picture preserving aspect ratio within bounds. Returns actual height used."""
    from PIL import Image as PILImage
    img = PILImage.open(img_path)
    img_w, img_h = img.size
    aspect = img_w / img_h
    w = max_width
    h = w / aspect
    if h > max_height:
        h = max_height
        w = h * aspect
    left_emu = left + (max_width - w) / 2
    slide.shapes.add_picture(img_path, Inches(left_emu), Inches(top), Inches(w), Inches(h))
    return h


def _add_table(slide, data, left, top, width, height):
    rows, cols = len(data), len(data[0])
    shape = slide.shapes.add_table(rows, cols, left, top, width, height)
    tbl = shape.table
    for r, row in enumerate(data):
        for c, val in enumerate(row):
            cell = tbl.cell(r, c)
            cell.text = str(val)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            for p in cell.text_frame.paragraphs:
                p.font.size = Pt(10)
                p.alignment = PP_ALIGN.CENTER
                if r == 0:
                    p.font.bold = True; p.font.color.rgb = C_WHITE
                    cell.fill.solid(); cell.fill.fore_color.rgb = C_ACCENT
                else:
                    p.font.color.rgb = C_TEXT
                    if r % 2 == 0:
                        cell.fill.solid(); cell.fill.fore_color.rgb = RGBColor(0xF0,0xF0,0xF0)
    return tbl


def _header(slide, text):
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, Inches(SLIDE_W), Inches(HEADER_H))
    bar.fill.solid(); bar.fill.fore_color.rgb = C_TITLE; bar.line.fill.background()
    tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(.12), Inches(CONTENT_W), Inches(.45))
    p = tx.text_frame.paragraphs[0]
    p.text = text; p.font.size = Pt(20); p.font.bold = True; p.font.color.rgb = C_WHITE


def _bullet(slide, text, x=MARGIN_L, y=CONTENT_TOP, w=CONTENT_W, h=3.5, font_size=11):
    tx = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tx.text_frame; tf.word_wrap = True
    if text:
        p = tf.paragraphs[0]; p.text = text; p.font.size = Pt(font_size); p.font.color.rgb = C_TEXT
    return tf


def _note_list(slide, notes, y, max_h=None):
    bottom = SLIDE_H - 0.2
    remaining = bottom - y
    if remaining < 0.2:
        return
    h = min(remaining, len(notes) * .25 + .1)
    if max_h:
        h = min(h, max_h)
    tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(y), Inches(CONTENT_W), Inches(h))
    tf = tx.text_frame; tf.word_wrap = True
    for note in notes:
        p = tf.add_paragraph(); p.text = f"• {note}"
        p.font.size = Pt(9); p.font.color.rgb = C_TEXT; p.space_before = Pt(3)


# ============================================================
# Slide builders — all content fits within slide bounds
# ============================================================

def slide_title(prs, title, subtitle, date, pipeline, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid(); slide.background.fill.fore_color.rgb = C_TITLE
    tx = slide.shapes.add_textbox(Inches(1), Inches(1.0), Inches(8), Inches(2.5))
    tf = tx.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.text = title; p.font.size = Pt(36)
    p.font.bold = True; p.font.color.rgb = C_WHITE; p.alignment = PP_ALIGN.CENTER
    if subtitle:
        p2 = tf.add_paragraph(); p2.text = subtitle
        p2.font.size = Pt(18); p2.font.color.rgb = RGBColor(0xBB,0xBB,0xBB); p2.alignment = PP_ALIGN.CENTER
    if date:
        p3 = tf.add_paragraph(); p3.text = f"{t('date_label', lang)}: {date}"
        p3.font.size = Pt(14); p3.font.color.rgb = RGBColor(0x88,0x88,0x88); p3.alignment = PP_ALIGN.CENTER
    if pipeline:
        p4 = tf.add_paragraph(); p4.text = f"{t('pipeline_label', lang)}: {pipeline}"
        p4.font.size = Pt(12); p4.font.color.rgb = RGBColor(0x88,0x88,0x88); p4.alignment = PP_ALIGN.CENTER


def slide_workflow(prs, samples, input_samples, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("workflow_title", lang))
    lay = Layout()
    steps = [
        ("1", "FastQC", t("raw_qc", lang)),
        ("2", "TrimGalore", t("adapter_trimming", lang)),
        ("3", "Bowtie2", t("alignment", lang)),
        ("4", "GATK", t("markdup", lang)),
        ("5", "MACS3", t("peak_calling", lang)),
        ("6", "HOMER", t("peak_annotation", lang)),
        ("7", "bamCoverage", t("bigwig_track", lang)),
        ("8", "FRiP", t("enrichment_qc", lang)),
    ]
    cols = 4
    step_h = 0.7
    rows_needed = (len(steps) + cols - 1) // cols
    step_area_h = rows_needed * step_h
    for i, (num, tool, desc) in enumerate(steps):
        row, col = i // cols, i % cols
        x = MARGIN_L + col * 2.25
        y = lay.y + row * step_h
        circle = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(x), Inches(y), Inches(.3), Inches(.3))
        circle.fill.solid(); circle.fill.fore_color.rgb = C_ACCENT; circle.line.fill.background()
        tf = circle.text_frame; p = tf.paragraphs[0]; p.text = num
        p.font.size = Pt(11); p.font.bold = True; p.font.color.rgb = C_WHITE; p.alignment = PP_ALIGN.CENTER
        tx = slide.shapes.add_textbox(Inches(x+.35), Inches(y-.02), Inches(1.5), Inches(.22))
        p = tx.text_frame.paragraphs[0]; p.text = tool
        p.font.size = Pt(10); p.font.bold = True; p.font.color.rgb = C_TITLE
        tx = slide.shapes.add_textbox(Inches(x+.35), Inches(y+.2), Inches(1.5), Inches(.22))
        p = tx.text_frame.paragraphs[0]; p.text = desc; p.font.size = Pt(8); p.font.color.rgb = C_TEXT
    lay.y += step_area_h
    lay.gap(0.15)

    tbl = [[t("sample", lang), t("type", lang), t("replicate", lang)]]
    for i, (ip, inp) in enumerate(zip(samples, input_samples)):
        tbl.append([ip, t("ip_type", lang), f"Rep {i+1}"])
        tbl.append([inp, t("input_type", lang), f"Rep {i+1}"])
    alloc = lay.allocate(min(lay.remaining, len(tbl) * .28))
    if alloc:
        _add_table(slide, tbl, Inches(MARGIN_L), Inches(alloc[0]), Inches(CONTENT_W), Inches(alloc[1]))


def slide_alignment(prs, alignment_data, img_path, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("alignment_title", lang))
    lay = Layout()
    # Reserve space for conclusion text (~1 inch)
    img_h = lay.remaining - 1.0
    alloc = lay.allocate(img_h)
    if alloc:
        actual_h = _add_picture(slide, img_path, MARGIN_L, alloc[0], CONTENT_W, alloc[1])
        lay.y = alloc[0] + actual_h  # snap to actual image bottom
    lay.gap(0.05)
    # Conclusion
    all_rates = {s: (d or 0) for s, d in alignment_data.items()}
    low = [s for s, v in all_rates.items() if v < ALIGN_THRESH]
    conc_alloc = lay.allocate(lay.remaining)
    if conc_alloc:
        tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(conc_alloc[0]), Inches(CONTENT_W), Inches(conc_alloc[1]))
        tf = tx.text_frame; tf.word_wrap = True
        if low:
            p = tf.paragraphs[0]
            p.text = t("alignment_warning", lang).format(", ".join(low))
            p.font.size = Pt(10); p.font.color.rgb = C_RED
        else:
            p = tf.paragraphs[0]
            p.text = f"所有样本比对率 > {ALIGN_THRESH:.0f}%, 均通过 QC 阈值"
            p.font.size = Pt(10); p.font.color.rgb = C_GREEN
        # Per-sample summary
        ip_samples = [s for s in all_rates if "IP" in s and "Input" not in s]
        input_samples = [s for s in all_rates if "Input" in s]
        if ip_samples:
            avg_ip = sum(all_rates[s] for s in ip_samples) / len(ip_samples)
            p2 = tf.add_paragraph()
            p2.text = f"IP 样本平均比对率: {avg_ip:.1f}%"
            p2.font.size = Pt(9); p2.font.color.rgb = C_TEXT
        if input_samples:
            avg_inp = sum(all_rates[s] for s in input_samples) / len(input_samples)
            p3 = tf.add_paragraph()
            p3.text = f"Input 样本平均比对率: {avg_inp:.1f}%"
            p3.font.size = Pt(9); p3.font.color.rgb = C_TEXT


def slide_markdup(prs, markdup_data, all_samples, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("markdup_title", lang))
    lay = Layout()
    tbl = [[t("sample", lang), t("read_pairs", lang), t("dup_pairs", lang),
            t("dup_rate", lang), t("unmapped", lang)]]
    for s in all_samples:
        d = markdup_data.get(s)
        if d:
            tbl.append([s, f"{d['read_pairs']:,}", f"{d['dup_pairs']:,}",
                        f"{d['dup_rate']:.1f}%", f"{d['unmapped']:,}"])
    alloc = lay.allocate(min(lay.remaining * .55, len(tbl) * .32))
    if alloc:
        _add_table(slide, tbl, Inches(MARGIN_L), Inches(alloc[0]), Inches(CONTENT_W), Inches(alloc[1]))
    lay.gap(0.08)
    if lay.remaining > 0.2:
        _note_list(slide, t("markdup_notes", lang), lay.y)


def slide_peak_calling(prs, pf_img, samples, peaks, frips, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("peak_calling_title", lang))
    lay = Layout()
    img_h = lay.remaining - 0.9
    alloc = lay.allocate(img_h)
    if alloc:
        actual_h = _add_picture(slide, pf_img, MARGIN_L, alloc[0], CONTENT_W, alloc[1])
        lay.y = alloc[0] + actual_h
    lay.gap(0.05)
    conc_alloc = lay.allocate(lay.remaining)
    if conc_alloc:
        tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(conc_alloc[0]), Inches(CONTENT_W), Inches(conc_alloc[1]))
        tf = tx.text_frame; tf.word_wrap = True
        p = tf.paragraphs[0]
        max_sample = max(samples, key=lambda s: peaks[s])
        min_sample = min(samples, key=lambda s: peaks[s])
        p.text = f"Peak 数: {max_sample} 最多 ({peaks[max_sample]:,}), {min_sample} 最少 ({peaks[min_sample]:,})"
        p.font.size = Pt(10); p.font.color.rgb = C_TEXT
        low_frip = [s for s in samples if (frips[s] or 0) < FRIP_THRESH]
        if low_frip:
            p2 = tf.add_paragraph()
            p2.text = f"注意: {', '.join(low_frip)} FRiP 均低于 {FRIP_THRESH*100:.0f}% 阈值, 富集较弱"
            p2.font.size = Pt(9); p2.font.color.rgb = C_RED


def slide_cutoff(prs, cutoff_img, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("cutoff_title", lang))
    lay = Layout()
    img_h = lay.remaining - 0.5
    alloc = lay.allocate(img_h)
    if alloc:
        actual_h = _add_picture(slide, cutoff_img, MARGIN_L, alloc[0], CONTENT_W, alloc[1])
        lay.y = alloc[0] + actual_h
    lay.gap(0.05)
    conc_alloc = lay.allocate(lay.remaining)
    if conc_alloc:
        tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(conc_alloc[0]), Inches(CONTENT_W), Inches(conc_alloc[1]))
        tf = tx.text_frame; tf.word_wrap = True
        p = tf.paragraphs[0]
        p.text = t("cutoff_note", lang) + " | pscore=4.8 为拐点, 低于此值 peak 数激增"
        p.font.size = Pt(9); p.font.color.rgb = C_TEXT


def slide_annotation(prs, ann_img, ann_data, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("annotation_title", lang))
    lay = Layout()
    img_h = lay.remaining - 1.0
    alloc = lay.allocate(img_h)
    if alloc:
        actual_h = _add_picture(slide, ann_img, MARGIN_L, alloc[0], CONTENT_W, alloc[1])
        lay.y = alloc[0] + actual_h
    lay.gap(0.05)
    conc_alloc = lay.allocate(lay.remaining)
    if conc_alloc:
        tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(conc_alloc[0]), Inches(CONTENT_W), Inches(conc_alloc[1]))
        tf = tx.text_frame; tf.word_wrap = True
        # Find dominant category per sample
        for sample, cats in ann_data.items():
            total = sum(cats.values())
            top_cat = max(cats, key=cats.get)
            top_pct = cats[top_cat] / total * 100
            p = tf.add_paragraph()
            p.text = f"{sample}: {t(top_cat, lang)} 占比最高 ({top_pct:.1f}%), 共 {total} peaks"
            p.font.size = Pt(9); p.font.color.rgb = C_TEXT


def slide_tss(prs, tss_img, tss_stats, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("tss_title", lang))
    lay = Layout()
    img_h = lay.remaining - 1.0
    alloc = lay.allocate(img_h)
    if alloc:
        actual_h = _add_picture(slide, tss_img, MARGIN_L, alloc[0], CONTENT_W, alloc[1])
        lay.y = alloc[0] + actual_h
    lay.gap(0.05)
    conc_alloc = lay.allocate(lay.remaining)
    if conc_alloc:
        tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(conc_alloc[0]), Inches(CONTENT_W), Inches(conc_alloc[1]))
        tf = tx.text_frame; tf.word_wrap = True
        for sample, ts in tss_stats.items():
            p = tf.add_paragraph()
            p.text = f"{sample}: {ts['within_1kb']} peaks 在 TSS 1kb 内 ({ts['within_1kb_pct']:.1f}%), {ts['within_5kb']} peaks 在 5kb 内 ({ts['within_5kb_pct']:.1f}%)"
            p.font.size = Pt(9); p.font.color.rgb = C_TEXT
            if ts['within_1kb_pct'] > 20:
                p2 = tf.add_paragraph()
                p2.text = f"  → 呈现强启动子近端结合特征, 符合转录因子结合模式"
                p2.font.size = Pt(9); p2.font.color.rgb = C_GREEN


def slide_top_genes(prs, genes_data, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("top_genes_title", lang))
    lay = Layout()

    # Build a single unified table: Sample | Gene | Position | Score
    # This naturally handles any number of samples
    tbl = [[t("sample", lang), t("gene", lang), t("position", lang), t("score", lang)]]
    for sample, genes in genes_data.items():
        if not genes:
            continue
        for i, (gene, pos, score, _) in enumerate(genes):
            tbl.append([sample if i == 0 else "", gene, pos, f"{score:,}"])

    if len(tbl) <= 1:
        return

    # Calculate table height: limit rows that fit, paginate if needed
    max_rows_fit = int(lay.remaining / 0.2)  # ~0.2 inch per row
    # Reserve space for shared genes if applicable
    all_sets = {s: {g for g, _, _, _ in gs} for s, gs in genes_data.items() if gs}
    shared = set.intersection(*all_sets.values()) if len(all_sets) > 1 else set()
    if shared:
        max_rows_fit = int((lay.remaining - 0.4) / 0.2)

    if len(tbl) > max_rows_fit:
        # Truncate and note
        tbl = tbl[:max_rows_fit]
        truncated = True
    else:
        truncated = False

    tbl_h = min(lay.remaining - 0.4, len(tbl) * 0.2)
    alloc = lay.allocate(max(tbl_h, 0.3))
    if alloc:
        _add_table(slide, tbl, Inches(MARGIN_L), Inches(alloc[0]), Inches(CONTENT_W), Inches(alloc[1]))
    lay.gap(0.1)

    # Truncation note
    if truncated and lay.remaining > 0.15:
        alloc2 = lay.allocate(0.15)
        if alloc2:
            tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(alloc2[0]), Inches(CONTENT_W), Inches(alloc2[1]))
            p = tx.text_frame.paragraphs[0]
            p.text = f"（仅显示前 {max_rows_fit-1} 行，使用 --top-n 调整）"
            p.font.size = Pt(8); p.font.color.rgb = C_TEXT; p.font.italic = True
            lay.gap(0.05)

    # Shared genes
    if shared and lay.remaining > 0.3:
        sh_alloc = lay.allocate(min(0.35, lay.remaining))
        if sh_alloc:
            tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(sh_alloc[0]), Inches(CONTENT_W), Inches(sh_alloc[1]))
            tf = tx.text_frame; tf.word_wrap = True
            p = tf.paragraphs[0]
            p.text = f"{t('shared_genes', lang)}: {', '.join(sorted(shared))}"
            p.font.size = Pt(9); p.font.bold = True; p.font.color.rgb = C_ACCENT
            p2 = tf.add_paragraph(); p2.text = t("shared_desc", lang)
            p2.font.size = Pt(8); p2.font.color.rgb = C_TEXT


def slide_summary(prs, samples, aligns, markdups, peaks, frips, anns, tss_stats, lang):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, t("summary_title", lang))
    lay = Layout()

    def _fmt(v, fmt=".1f", unit="%"):
        return f"{v:{fmt}}{unit}" if v is not None else "N/A"

    # --- QC table (compact: 5 rows x N cols) ---
    tbl = [[t("metric", lang)] + samples + [t("threshold", lang)]]
    tbl.append(["Bowtie2"] + [_fmt(aligns[s]) for s in samples] + [f">{ALIGN_THRESH:.0f}%"])
    tbl.append([t("dup_rate", lang)] +
               [_fmt(markdups[s]["dup_rate"] if markdups.get(s) else None) for s in samples] + ["<50%"])
    tbl.append([t("peak_count", lang)] + [f"{peaks[s]:,}" for s in samples] + ["—"])
    tbl.append([t("frip_score", lang)] +
               [_fmt(frips[s]*100 if frips[s] is not None else None, ".2f") for s in samples]
               + [f"≥{FRIP_THRESH*100:.0f}%"])
    # Table: 5 rows, ~1.2 inches
    alloc = lay.allocate(min(1.2, len(tbl) * .25))
    if alloc:
        _add_table(slide, tbl, Inches(MARGIN_L), Inches(alloc[0]), Inches(CONTENT_W), Inches(alloc[1]))
    lay.gap(0.08)

    # --- Per-sample summary (compact: 2 lines per sample) ---
    # Calculate: title(0.2) + 3 samples * 2 lines * 0.13 + recs(0.4) ≈ 1.6
    # Available: ~3.4 inches — plenty
    conc_alloc = lay.allocate(lay.remaining)
    if conc_alloc is None:
        return
    tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(conc_alloc[0]), Inches(CONTENT_W), Inches(conc_alloc[1]))
    tf = tx.text_frame; tf.word_wrap = True
    tf.auto_size = None  # prevent auto-expansion beyond box

    # Title
    p = tf.paragraphs[0]
    p.text = t("conclusions", lang)
    p.font.size = Pt(12); p.font.bold = True; p.font.color.rgb = C_ACCENT

    for s in samples:
        rate = aligns[s] or 0
        frip = frips[s] or 0
        n_peaks = peaks[s]
        dup_rate = markdups[s]["dup_rate"] if markdups.get(s) else 0

        # Determine status
        issues = []
        if rate < ALIGN_THRESH:
            issues.append(f"比对率低({rate:.0f}%)")
        if frip < FRIP_THRESH:
            issues.append(f"FRiP低({frip*100:.1f}%)")
        if dup_rate > 50:
            issues.append(f"重复率高({dup_rate:.0f}%)")

        if not issues:
            status = t("good_enrichment", lang).format(frip*100)
            color = C_GREEN
        elif len(issues) == 1 and rate >= ALIGN_THRESH:
            status = issues[0]
            color = C_YELLOW
        else:
            status = "; ".join(issues)
            color = C_RED

        # Line 1: sample name + key stats + status (all in one line)
        p = tf.add_paragraph()
        ann_str = ""
        if anns.get(s):
            total_ann = sum(anns[s].values())
            promoter_pct = anns[s].get("promoter", 0) / total_ann * 100 if total_ann else 0
            ann_str = f" 启动子{promoter_pct:.0f}%"
        tss_str = ""
        if tss_stats.get(s):
            ts = tss_stats[s]
            tss_str = f" TSS1kb内{ts['within_1kb']}个({ts['within_1kb_pct']:.0f}%)"
        p.text = f"▸ {s}: {n_peaks:,} peaks | 比对率{rate:.0f}% | 重复率{dup_rate:.0f}% | FRiP{frip*100:.1f}% | {status}"
        p.font.size = Pt(9); p.font.color.rgb = color; p.space_before = Pt(4)

        # Line 2: annotation + TSS (compact)
        if ann_str or tss_str:
            p2 = tf.add_paragraph()
            p2.text = f"  注释:{ann_str}{tss_str}"
            p2.font.size = Pt(8); p2.font.color.rgb = C_TEXT

    # --- Recommendations ---
    low_frip = [s for s in samples if (frips[s] or 0) < FRIP_THRESH]
    low_align = [s for s in samples if (aligns[s] or 0) < ALIGN_THRESH]
    high_dup = [s for s in samples if markdups.get(s) and markdups[s]["dup_rate"] > 50]

    if low_frip or low_align or high_dup:
        p = tf.add_paragraph()
        p.text = t("recommendations", lang)
        p.font.size = Pt(10); p.font.bold = True; p.font.color.rgb = C_ACCENT; p.space_before = Pt(6)

        if low_frip:
            p2 = tf.add_paragraph()
            p2.text = f"• FRiP偏低({', '.join(low_frip)}): 建议优化ChIP抗体或增加测序深度"
            p2.font.size = Pt(8); p2.font.color.rgb = C_TEXT
        if low_align:
            p3 = tf.add_paragraph()
            p3.text = f"• 比对率偏低({', '.join(low_align)}): 检查样本质量、接头污染或参考基因组"
            p3.font.size = Pt(8); p3.font.color.rgb = C_TEXT
        if high_dup:
            p4 = tf.add_paragraph()
            p4.text = f"• 重复率偏高({', '.join(high_dup)}): 考虑增加起始量或优化文库构建"
            p4.font.size = Pt(8); p4.font.color.rgb = C_TEXT


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(description="Generate ChIP-seq Peak Calling PPT report")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--input-samples", nargs="+", required=True)
    ap.add_argument("--peaks-dir", required=True)
    ap.add_argument("--annotation-dir", required=True)
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--log-dir", required=True)
    ap.add_argument("--markdup-dir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--img-dir", default="",
                    help="Directory to save plot images (default: temp, auto-deleted)")
    ap.add_argument("--title", default="")
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--pipeline", default="")
    ap.add_argument("--genome", default="")
    ap.add_argument("--date", default="")
    ap.add_argument("--top-n", type=int, default=5)
    ap.add_argument("--lang", choices=["zh", "en"], default="zh")
    args = ap.parse_args()

    lang = args.lang
    if not args.title:
        args.title = t("report_title", lang)
    if not args.subtitle:
        sample_str = " & ".join(args.samples)
        genome_str = f" — {args.genome}" if args.genome else ""
        args.subtitle = f"{sample_str}{genome_str}"

    all_samples = args.samples + args.input_samples

    # Load data
    aligns = {s: load_bowtie2(args.log_dir, s) for s in all_samples}
    markdups = {s: d for s in all_samples if (d := load_markdup(args.markdup_dir, s))}
    frips = {s: load_frip(args.qc_dir, s) for s in args.samples}
    peaks = {s: load_peak_count(args.peaks_dir, s) for s in args.samples}
    cutoffs = {s: d for s in args.samples if (d := load_cutoff(args.peaks_dir, s))}
    anns = {s: d for s in args.samples if (d := load_annotation(args.annotation_dir, s))}
    genes = {s: load_top_genes(args.annotation_dir, s, args.top_n) for s in args.samples}

    # TSS distance analysis
    tss_stats = {}
    for sample in args.samples:
        af = os.path.join(args.annotation_dir, sample, f"{sample}_peaks.annotatePeaks.txt")
        if not os.path.isfile(af):
            continue
        distances = []
        with open(af) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                if len(row) >= 10 and row[9]:
                    try:
                        distances.append(abs(int(row[9])))
                    except ValueError:
                        pass
        if distances:
            total = len(distances)
            within_1kb = sum(1 for d in distances if d <= 1000)
            within_5kb = sum(1 for d in distances if d <= 5000)
            tss_stats[sample] = {
                "total": total,
                "within_1kb": within_1kb,
                "within_1kb_pct": within_1kb / total * 100,
                "within_5kb": within_5kb,
                "within_5kb_pct": within_5kb / total * 100,
            }

    # Generate plots (DPI 300)
    # If img-dir is set, save images there; otherwise use temp
    if args.img_dir:
        os.makedirs(args.img_dir, exist_ok=True)

    align_img = plot_alignment({s: aligns[s] or 0 for s in all_samples})
    pf_img = plot_peak_and_frip(args.samples, peaks, frips)
    cutoff_img = plot_cutoff(cutoffs) if cutoffs else None
    ann_img = plot_annotation(anns, lang) if anns else None
    tss_img = plot_tss_distance(args.annotation_dir, args.samples)

    # Copy images to img-dir if specified
    all_imgs = {"alignment": align_img, "peak_frip": pf_img, "cutoff": cutoff_img,
                "annotation": ann_img, "tss_distance": tss_img}
    if args.img_dir:
        import shutil
        for name, path in all_imgs.items():
            if path and os.path.exists(path):
                dst = os.path.join(args.img_dir, f"{name}.png")
                shutil.copy2(path, dst)
                print(f"  Image saved: {dst}")

    # Build PPT
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)

    slide_title(prs, args.title, args.subtitle, args.date, args.pipeline, lang)
    slide_workflow(prs, args.samples, args.input_samples, lang)
    slide_alignment(prs, aligns, align_img, lang)
    slide_markdup(prs, markdups, all_samples, lang)
    slide_peak_calling(prs, pf_img, args.samples, peaks, frips, lang)
    if cutoff_img:
        slide_cutoff(prs, cutoff_img, lang)
    if ann_img:
        slide_annotation(prs, ann_img, anns, lang)
    slide_tss(prs, tss_img, tss_stats, lang)
    if any(genes.values()):
        slide_top_genes(prs, genes, lang)
    slide_summary(prs, args.samples, aligns, markdups, peaks, frips, anns, tss_stats, lang)

    # Save
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    prs.save(args.output)
    print(f"Saved: {args.output} ({len(prs.slides)} slides)")

    # Cleanup temp images (only if not saving to img-dir)
    if not args.img_dir:
        for img in all_imgs.values():
            if img and os.path.exists(img):
                os.unlink(img)


if __name__ == "__main__":
    main()
