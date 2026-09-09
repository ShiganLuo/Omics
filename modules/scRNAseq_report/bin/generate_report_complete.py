#!/usr/bin/env python3
"""Generate comprehensive scRNA-seq analysis report PPT - Complete version."""
import os
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt
from PIL import Image as PILImage

# Colors
C_NAVY = RGBColor(0x18, 0x25, 0x43)
C_ACCENT = RGBColor(0x00, 0x94, 0xD8)
C_GREEN = RGBColor(0x00, 0xA8, 0x78)
C_RED = RGBColor(0xD6, 0x45, 0x45)
C_ORANGE = RGBColor(0xF0, 0x9A, 0x36)
C_TEXT = RGBColor(0x33, 0x33, 0x33)
C_MUTED = RGBColor(0x74, 0x74, 0x74)
C_WHITE = RGBColor(0xFF, 0xFF, 0xFF)

SLIDE_W = 10.0
SLIDE_H = 5.625
HEADER_H = 0.65
MARGIN_L = 0.45
CONTENT_W = SLIDE_W - MARGIN_L - 0.45

BASE_DIR = Path("/home/luosg/Data/genomeStability/output/luancao/scRNAseq/cell_type_comparison")
HYSTERA_DIR = Path("/home/luosg/Data/genomeStability/output/luancao/scRNAseq/common/5_combine_h5ad/Uterus")
OVARIES_DIR = Path("/home/luosg/Data/genomeStability/output/luancao/scRNAseq/common/5_combine_h5ad/ovaries")


def _header(slide, text):
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, Inches(SLIDE_W), Inches(HEADER_H))
    bar.fill.solid()
    bar.fill.fore_color.rgb = C_NAVY
    bar.line.fill.background()
    tx = slide.shapes.add_textbox(Inches(MARGIN_L), Inches(0.08), Inches(CONTENT_W), Inches(0.45))
    tf = tx.text_frame
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(22)
    p.font.bold = True
    p.font.color.rgb = C_WHITE


def _textbox(slide, left, top, width, height, text, font_size=11, bold=False, color=C_TEXT, align=PP_ALIGN.LEFT):
    box = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.alignment = align
    return box


def _bullets(slide, left, top, width, height, items, font_size=11, color=C_TEXT):
    box = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = box.text_frame
    tf.word_wrap = True
    first = True
    for item in items:
        if first:
            p = tf.paragraphs[0]
            first = False
        else:
            p = tf.add_paragraph()
        p.text = f"• {item}"
        p.font.size = Pt(font_size)
        p.font.color.rgb = color
        p.space_after = Pt(4)
    return box


def _table(slide, left, top, width, height, data, font_size=9):
    rows = len(data)
    cols = len(data[0]) if data else 1
    shape = slide.shapes.add_table(rows, cols, Inches(left), Inches(top), Inches(width), Inches(height))
    tbl = shape.table
    for r, row in enumerate(data):
        for c, value in enumerate(row):
            cell = tbl.cell(r, c)
            cell.text = str(value)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            if r == 0:
                cell.fill.solid()
                cell.fill.fore_color.rgb = C_ACCENT
            elif r % 2 == 0:
                cell.fill.solid()
                cell.fill.fore_color.rgb = RGBColor(0xF1, 0xF5, 0xFA)
            for p in cell.text_frame.paragraphs:
                p.font.size = Pt(font_size)
                p.alignment = PP_ALIGN.CENTER
                if r == 0:
                    p.font.bold = True
                    p.font.color.rgb = C_WHITE
                else:
                    p.font.color.rgb = C_TEXT
    return tbl


def _add_picture(slide, path, left, top, max_w, max_h):
    """Add picture maintaining aspect ratio."""
    if not path or not os.path.isfile(path):
        return None
    img = PILImage.open(path)
    aspect = img.size[0] / max(img.size[1], 1)
    width = max_w
    height = width / aspect
    if height > max_h:
        height = max_h
        width = height * aspect
    x = left + (max_w - width) / 2
    y = top + (max_h - height) / 2
    slide.shapes.add_picture(path, Inches(x), Inches(y), Inches(width), Inches(height))
    return width, height


def get_all_proportions():
    """Get cell type proportions from all4 datasets."""
    datasets = {
        'Hystera CR': HYSTERA_DIR / "Uterus_cellranger_v3_annotated.h5ad",
        'Hystera scTE': HYSTERA_DIR / "Uterus_scTE_v2_annotated.h5ad",
        'Ovaries CR': OVARIES_DIR / "ovaries_cellranger_v3_annotated.h5ad",
        'Ovaries scTE': OVARIES_DIR / "ovaries_scTE_v3_annotated.h5ad",
    }
    
    results = {}
    for name, path in datasets.items():
        if path.exists():
            adata = ad.read_h5ad(str(path))
            ct_counts = adata.obs['cell_type'].value_counts()
            results[name] = {ct: cnt/len(adata)*100 for ct, cnt in ct_counts.items()}
        else:
            results[name] = {}
    
    return results


# ── Slides ──────────────────────────────────────────────────────────────

def create_title_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    bg = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, Inches(SLIDE_W), Inches(SLIDE_H))
    bg.fill.solid()
    bg.fill.fore_color.rgb = C_NAVY
    bg.line.fill.background()
    
    _textbox(slide, 1, 1.2, 8, 1, "scRNA-seq 二次过滤与差异丰度分析", 28, True, C_WHITE, PP_ALIGN.CENTER)
    _textbox(slide, 1, 2.2, 8, 0.5, "恒河猴子宫（Hystera）与卵巢（Ovaries）", 20, False, C_ACCENT, PP_ALIGN.CENTER)
    _textbox(slide, 1, 3.0, 8, 0.5, "Cell Ranger + scTE 双平台定量", 16, False, RGBColor(0xBB, 0xBB, 0xBB), PP_ALIGN.CENTER)
    _textbox(slide, 1, 4.0, 8, 1, 
             "二次过滤流程建立 | 新细胞类型发现 | Marker基因添加 | 差异丰度分析", 
             13, False, RGBColor(0x99, 0x99, 0x99), PP_ALIGN.CENTER)
    _textbox(slide, 1, 4.8, 8, 0.3, "2026-09-07", 12, False, RGBColor(0x77, 0x77, 0x77), PP_ALIGN.CENTER)
    return slide


def create_workflow_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "分析流程")
    
    # Left: upstream
    _textbox(slide, 0.3, 0.9, 4.5, 0.3, "上游分析", 14, True, C_ACCENT)
    _bullets(slide, 0.3, 1.3, 4.5, 2, [
        "Cell Ranger: 10x Genomics标准基因表达定量",
        "scTE: 转座子（TE）表达定量",
        "质控 → 合并 → 聚类 → 注释",
        "输出：_merged.h5ad → _clustered.h5ad → _annotated.h5ad",
        "工具：scRNAseq.py（mode: qc/merge/cluster/annotate）",
    ], 10)
    
    # Right: downstream
    _textbox(slide, 5, 0.9, 4.5, 0.3, "下游二次过滤", 14, True, C_RED)
    _bullets(slide, 5, 1.3, 4.5, 2.5, [
        "1. Cluster QC分析",
        "   • 检查每个cluster的基因数、UMI、marker",
        "   • 标记低质量cluster",
        "",
        "2. 逐细胞过滤（非整cluster删除）",
        "   • 在标记cluster中过滤低质量细胞",
        "   • 参数因组织和定量方式而异",
        "",
        "3. 重新聚类 + 手动注释",
        "   • scRNAseq.py重新聚类",
        "   • annotate_all.py注释",
        "   • --manual-annotation手动指定",
    ], 9)
    
    # Bottom: parameters
    _textbox(slide, 0.3, 4.5, 9, 0.3, "过滤参数示例（因数据而异）", 11, True, C_ORANGE)
    _textbox(slide, 0.3, 4.8, 9, 0.6, 
             "Hystera: min_genes=800, min_counts=3000, max_pct_mt=20%\n"
             "Ovaries scTE: min_genes=800, min_counts=3000（不使用MT%，数据整体MT%很低）", 
             9, False, C_MUTED)
    
    return slide


def create_Uterus_proportions_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "Hystera 细胞类型比例")
    
    # CellRanger UMAP
    img_cr = HYSTERA_DIR / "plots" / "cellranger" / "cluster_annotate" / "umap_cell_type.png"
    _add_picture(slide, str(img_cr), 0.2, 0.8, 3, 2.2)
    _textbox(slide, 0.2, 3.0, 3, 0.2, "Cell Ranger", 10, True, C_ACCENT, PP_ALIGN.CENTER)
    
    # scTE UMAP
    img_sc = HYSTERA_DIR / "plots" / "scTE" / "cluster_annotate" / "umap_cell_type.png"
    _add_picture(slide, str(img_sc), 3.4, 0.8, 3, 2.2)
    _textbox(slide, 3.4, 3.0, 3, 0.2, "scTE", 10, True, C_ACCENT, PP_ALIGN.CENTER)
    
    # Proportion table
    data = [
        ["Cell Type", "CR %", "scTE %"],
        ["Stromal", "29.1", "29.8"],
        ["Endothelial", "16.4", "16.4"],
        ["Smooth_muscle", "15.6", "17.1"],
        ["Epithelial_glandular", "14.8", "4.4"],
        ["Pericyte", "12.4", "9.8"],
        ["Uterine_NK", "7.0", "6.7"],
        ["Epithelial_luminal", "1.8", "12.8"],
        ["Macrophage", "2.3", "2.2"],
        ["Schwann_cell", "0.6", "0.4"],
        ["Proliferating", "-", "0.4"],
    ]
    _table(slide, 6.6, 0.8, 3.2, 2.8, data, 8)
    
    _textbox(slide, 0.2, 3.3, 9.5, 1.5, 
             "总计：Cell Ranger 39,337细胞（9种） | scTE 41,009细胞（10种）\n"
             "注意：scTE的Epithelial_luminal比例（12.8%）远高于Cell Ranger（1.8%），可能与TE定量方式有关", 
             9, False, C_TEXT)
    
    return slide


def create_ovaries_proportions_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "Ovaries 细胞类型比例")
    
    # CellRanger UMAP
    img_cr = OVARIES_DIR / "plots" / "cellranger" / "cluster_annotate" / "umap_cell_type.png"
    _add_picture(slide, str(img_cr), 0.2, 0.8, 3, 2.2)
    _textbox(slide, 0.2, 3.0, 3, 0.2, "Cell Ranger", 10, True, C_ACCENT, PP_ALIGN.CENTER)
    
    # scTE UMAP
    img_sc = OVARIES_DIR / "plots" / "scTE" / "cluster_annotate" / "umap_cell_type.png"
    _add_picture(slide, str(img_sc), 3.4, 0.8, 3, 2.2)
    _textbox(slide, 3.4, 3.0, 3, 0.2, "scTE", 10, True, C_ACCENT, PP_ALIGN.CENTER)
    
    # Proportion table
    data = [
        ["Cell Type", "CR %", "scTE %"],
        ["Luteal", "22.8", "2.0"],
        ["Endothelial", "20.4", "17.8"],
        ["Stromal", "19.3", "21.6"],
        ["Smooth_muscle", "6.5", "10.3"],
        ["Myofibroblast", "9.2", "2.6"],
        ["Pericyte", "6.6", "5.5"],
        ["Cumulus", "5.0", "4.5"],
        ["Lymphatic_endo", "1.6", "1.5"],
        ["Mesothelial", "1.6", "1.7"],
        ["Epithelial", "2.0", "1.3"],
        ["Macrophage", "2.7", "2.4"],
        ["T_cell", "1.6", "1.0"],
        ["Proliferating", "0.7", "2.0"],
        ["ERVK_high", "-", "22.6"],
        ["Alu_high", "-", "3.2"],
    ]
    _table(slide, 6.6, 0.8, 3.2, 3.5, data, 8)
    
    _textbox(slide, 0.2, 3.3, 9.5, 1.5, 
             "总计：Cell Ranger 17,395细胞（13种） | scTE 14,943细胞（15种）\n"
             "scTE特有：ERVK_high（22.6%）、Alu_high（3.2%）\n"
             "Luteal在CR中22.8%，scTE中仅2.0%，被ERVK_high\"吞噬\"", 
             9, False, C_TEXT)
    
    return slide


def create_new_celltypes_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "新发现的细胞类型")
    
    data = [
        ["Cell Type", "Markers", "文献来源", "组织"],
        ["Myofibroblast", "ACTA2, MYH11, TAGLN, POSTN, IGFBP5, SFRP1", "PMC 2024", "H+O"],
        ["Lymphatic_endo", "MMRN1, CCL21, PROX1, LYVE1, PDPN, CAVIN2, FLT4", "Wigle 1999", "O"],
        ["Mesothelial", "MSLN, ITLN1, WT1, UPK3B", "Nature 2022", "O"],
        ["Cumulus", "FST, NR5A2, PPARG, CRHBP, GRB14, HAS2, PTX3", "卵巢卵丘", "O"],
        ["Luteal", "STAR, CYP11A1, HSD3B1, PTCH2, GPC5", "黄体细胞", "O"],
    ]
    _table(slide, 0.3, 0.9, 9.4, 2, data, 9)
    
    # Add dotplot
    img_path = HYSTERA_DIR / "plots" / "cellranger" / "cluster_annotate" / "dotplot_Uterus_cellranger_markers.png"
    if img_path.exists():
        _add_picture(slide, str(img_path), 0.3, 3.0, 9, 2.3)
    
    return slide


def create_Uterus_markers_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "Marker基因列表 - Hystera")
    
    markers = [
        ["Cell Type", "Markers"],
        ["Epithelial_luminal", "EPCAM, KRT18, KRT8, CDH1, CLDN4, PAEP, SLC2A1"],
        ["Epithelial_glandular", "EPCAM, KRT18, KRT8, PAEP, SPP1, LIF, MUC1, PAX8, ESR1, EYA2, BMPR1B, ENPP3"],
        ["Stromal", "COL1A1, COL1A2, DCN, LUM, PDGFRA, VIM, WNT4, IGFBP1"],
        ["Smooth_muscle", "ACTA2, MYH11, TAGLN, CNN1, MYLK"],
        ["Pericyte", "RGS5, PDGFRB, NOTCH3, MCAM, ABCC9, KCNJ8"],
        ["Schwann_cell", "PMP22, S100B, CDH19, PLP1, MPZ, CRYAB, MAL, GPM6B"],
        ["Endothelial", "PECAM1, VWF, CDH5, ERG, KDR, FLT1"],
        ["Lymphatic_endo", "MMRN1, CCL21, PROX1, LYVE1, PDPN, CAVIN2, FLT4"],
        ["Macrophage", "CD68, CD163, CSF1R, MRC1, MARCO, LYZ, S100A8"],
        ["Uterine_NK", "NKG7, GNLY, KLRB1, NCAM1, KLRC1, CSF1, XCL1"],
        ["T_cell", "CD3E, CD3D, CD3G, CD4, CD8A, IL7R"],
        ["Myofibroblast", "ACTA2, MYH11, TAGLN, POSTN, IGFBP5, SFRP1"],
        ["Mesothelial", "MSLN, ITLN1, WT1, UPK3B"],
        ["Proliferating", "MKI67, TOP2A, STMN1, HMGB2, NUSAP1, CENPF, TYMS, UBE2C"],
    ]
    _table(slide, 0.3, 0.9, 9.4, 4.5, markers, 8)
    
    return slide


def create_ovaries_markers_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "Marker基因列表 - Ovaries")
    
    markers = [
        ["Cell Type", "Markers"],
        ["Oocyte", "ZP1, ZP2, ZP3, ZP4, FIGLA, NOBOX, GDF9, BMP15, DAZL, DDX4"],
        ["Granulosa", "FOXL2, CYP19A1, FSHR, AMH, HSD17B1, INHA"],
        ["Cumulus", "FST, NR5A2, PPARG, CRHBP, GRB14, HAS2, PTX3"],
        ["Theca", "CYP17A1, STAR, DLK1, INSL3"],
        ["Luteal", "STAR, CYP11A1, HSD3B1, PTCH2, GPC5"],
        ["Stromal", "COL1A1, COL1A2, DCN, LUM, PDGFRA, VIM"],
        ["Endothelial", "PECAM1, VWF, CDH5, ERG, KDR, FLT1, EMCN"],
        ["Lymphatic_endo", "MMRN1, CCL21, PROX1, LYVE1, PDPN, CAVIN2, FLT4"],
        ["Smooth_muscle", "ACTA2, MYH11, TAGLN, CNN1, MYLK"],
        ["Epithelial", "EPCAM, KRT18, KRT8, CDH1, CLDN4"],
        ["Pericyte", "RGS5, PDGFRB, NOTCH3, MCAM"],
        ["Macrophage", "CD68, CD163, CSF1R, MRC1, MARCO, LYZ, S100A8, S100A9"],
        ["T_cell", "CD3E, CD3D, CD3G, CD4, CD8A, IL7R, TRAC"],
        ["Mesothelial", "MSLN, ITLN1, WT1, UPK3B"],
        ["Myofibroblast", "ACTA2, MYH11, TAGLN, POSTN, IGFBP5, SFRP1"],
        ["Proliferating", "MKI67, TOP2A, STMN1, HMGB2, NUSAP1, CENPF, TYMS, UBE2C"],
        ["ERVK_high", "MacERVK2_LTR1c, MacNERVK2-int, MacERVK2_LTR1a, MacNERVK2a-int"],
        ["Alu_high", "AluSz, AluSx, AluSx1, AluY, AluJb, AluSp, AluSg"],
    ]
    _table(slide, 0.3, 0.9, 9.4, 4.5, markers, 8)
    
    return slide


def create_da_ovaries_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "差异丰度分析 - Ovaries (21310 vs 11238)")
    
    img1 = BASE_DIR / "ovaries_volcano.png"
    _add_picture(slide, str(img1), 0.2, 0.8, 4.5, 3.5)
    
    img2 = BASE_DIR / "ovaries_bar_comparison.png"
    _add_picture(slide, str(img2), 4.9, 0.8, 4.5, 3.5)
    
    _textbox(slide, 0.2, 4.4, 9.5, 0.8, 
             "Chi-square: p = 0.00 (极显著) | ↑ ERVK_high(+3.2x), Cumulus(+2.9x), Alu_high(+1.9x)\n"
             "↓ Epithelial(-90x), Mesothelial(-9x), Smooth_muscle(-4x), Pericyte(-2.3x)", 
             10, True, C_RED)
    
    return slide


def create_da_Uterus_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "差异丰度分析 - Hystera")
    
    img = BASE_DIR / "Uterus_heatmap.png"
    _add_picture(slide, str(img), 0.2, 0.8, 9.5, 2.5)
    
    data = [
        ["比较", "增加", "减少"],
        ["11238 vs 2200671", "Stromal +6.2%, Uterine_NK +4.4%", "Endothelial -6.5%, Epithelial_glandular -5.4%"],
        ["11238 vs 1411200", "Epithelial_glandular +16.1%, Stromal +10.0%", "Endothelial -17.8%, Smooth_muscle -15.9%"],
        ["21310 vs 21224", "Stromal +24.2%, Epithelial_luminal +1.7%", "Epithelial_glandular -9.5%, Endothelial -8.7%"],
    ]
    _table(slide, 0.2, 3.4, 9.5, 1.8, data, 9)
    
    return slide


def create_conclusion_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _header(slide, "总结")
    
    _textbox(slide, 0.3, 0.9, 4.5, 0.3, "主要成果", 14, True, C_GREEN)
    _bullets(slide, 0.3, 1.3, 4.5, 3, [
        "建立二次过滤流程",
        "发现5种新细胞类型",
        "添加新marker到注释系统",
        "完成4组数据的差异丰度分析",
        "手动注释功能（--manual-annotation）",
    ], 11)
    
    _textbox(slide, 5, 0.9, 4.5, 0.3, "新发现的细胞类型", 14, True, C_ACCENT)
    _bullets(slide, 5, 1.3, 4.5, 3, [
        "Myofibroblast: SM+Stromal混合型",
        "Lymphatic_endo: 淋巴内皮细胞",
        "Mesothelial: 间皮细胞",
        "Cumulus: 卵丘细胞",
        "Luteal: 黄体细胞",
    ], 11)
    
    _textbox(slide, 0.3, 4.2, 9, 0.3, "后续工作", 14, True, C_ORANGE)
    _textbox(slide, 0.3, 4.6, 9, 0.8, 
             "• ERVK_high/Alu_high的生物学意义需进一步研究\n• 可尝试Milo方法进行更精确的差异丰度分析\n• 结合ATAC-seq验证细胞类型注释", 
             10, False, C_TEXT)
    
    return slide


def main():
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)
    
    create_title_slide(prs)
    create_workflow_slide(prs)
    create_Uterus_proportions_slide(prs)
    create_ovaries_proportions_slide(prs)
    create_new_celltypes_slide(prs)
    create_Uterus_markers_slide(prs)
    create_ovaries_markers_slide(prs)
    create_da_ovaries_slide(prs)
    create_da_Uterus_slide(prs)
    create_conclusion_slide(prs)
    
    output_path = BASE_DIR / "scRNAseq_analysis_report.pptx"
    prs.save(str(output_path))
    print(f"Report saved to: {output_path}")


if __name__ == "__main__":
    main()
