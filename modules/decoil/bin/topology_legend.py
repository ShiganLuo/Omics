#!/usr/bin/env python3
"""Generate the DECOIL topology legend (7 ecDNA topologies) with matplotlib.

All seven topologies are drawn as closed rings; graphic language:
color = fragment / source chromosome, narrow white slit = the single
excision breakpoint (simple_circle only), orange dashed arc = SV segment,
green arc with counter-clockwise chevrons = inversion, dark/light red =
duplicated copies, nested orange arcs = foldback.

Usage:
    python topology_legend.py [out.png]
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, Wedge

# Register a CJK font (the container ships only Latin fonts)
_CJK_FONT = "/tmp/fonts/NotoSansCJK-Regular.ttc"
try:
    font_manager.fontManager.addfont(_CJK_FONT)
    _CJK_NAME = font_manager.FontProperties(fname=_CJK_FONT).get_name()
except Exception:
    _CJK_NAME = "sans-serif"


# Muted palette, consistent with decoil_structure.py
BLUE = "#5B7C99"
BLUE_L = "#7A93B5"
BLUE_LL = "#9FB3CC"
ORANGE = "#D98E4A"
GREEN = "#8FB57A"
GREEN_D = "#3D6B4F"
RED = "#B5535C"
ROSE = "#B57A7F"
GRAY = "#8A94A0"
INK = "#3B4A5A"

R = 1.0          # ring outer radius
BAND = 0.22      # ring band width
SLIT = 4.0       # breakpoint slit width (degrees) for simple_circle


def add_ring(ax, arcs):
    """Draw ring arcs.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    arcs : list of (theta1, theta2, facecolor, alpha)
        Angles in degrees, counterclockwise, 0 = east.
    """
    for t1, t2, fc, alpha in arcs:
        ax.add_patch(Wedge((0, 0), R, t1, t2, width=BAND,
                           facecolor=fc, edgecolor="none", alpha=alpha))


def dashed_arc(ax, t1, t2, fc=ORANGE, dash=6.0, gap=4.5):
    """Draw a dashed arc segment inside the band (mimics a dashed stroke)."""
    t = t1
    while t < t2:
        ax.add_patch(Wedge((0, 0), R, t, min(t + dash, t2), width=BAND,
                           facecolor=fc, edgecolor="none"))
        t += dash + gap


def chevron(ax, theta_deg, fc="white", ccw=True,
            length=0.30, head_len=0.13, head_w=0.18, shaft_w=0.075):
    """Direction arrow centered on the band: rectangular shaft + triangular
    head.  Tip points along the ring direction when *ccw* else against it."""
    import numpy as np

    raise NotImplementedError("replaced by arrowed_arc")


def arrowed_arc(ax, t1, t2, fc, alpha=1.0, tip_deg=10.0, ccw=True):
    """Ring arc whose end tapers into an arrowhead (plasmid style).

    The arc spans t1 -> t2 (t1 < t2); with ``ccw=True`` the tip sits at *t2*
    (reading counter-clockwise), with ``ccw=False`` at *t1* (clockwise).
    Shape itself carries direction — no overlaid markers needed.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    t1, t2 : float
        Start / end angles in degrees (ccw, 0 = east); t1 < t2.
    fc : str
        Face color.
    alpha : float
        Face alpha.
    tip_deg : float
        Angular length of the taper.
    ccw : bool
        Tip direction: True = tip at t2 (ccw), False = tip at t1 (cw).
    """
    import numpy as np

    body_lo, body_hi = ((t1, t2 - tip_deg) if ccw else (t1 + tip_deg, t2))
    tip_at = t2 if ccw else t1
    n = max(8, int(t2 - t1))
    ts = np.linspace(body_lo, body_hi, n)
    outer = [(R * np.cos(np.radians(t)), R * np.sin(np.radians(t))) for t in ts]
    inner = [((R - BAND) * np.cos(np.radians(t)),
              (R - BAND) * np.sin(np.radians(t))) for t in ts[::-1]]
    ax.add_patch(plt.Polygon(outer + inner, closed=True, facecolor=fc,
                             edgecolor="none", alpha=alpha))
    rad_b = np.radians(body_hi if ccw else body_lo)
    rad_t = np.radians(tip_at)
    o = (R * np.cos(rad_b), R * np.sin(rad_b))
    i = ((R - BAND) * np.cos(rad_b), (R - BAND) * np.sin(rad_b))
    tip = ((R - BAND / 2.0) * np.cos(rad_t), (R - BAND / 2.0) * np.sin(rad_t))
    ax.add_patch(plt.Polygon([o, tip, i], closed=True, facecolor=fc,
                             edgecolor="none", alpha=alpha))


def pol(theta_deg, r):
    """Polar to xy helper."""
    import numpy as np

    t = np.radians(theta_deg)
    return r * np.cos(t), r * np.sin(t)


def style_axes(ax, title, desc):
    ax.set_xlim(-1.85, 1.85)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor("white")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0, 1.36, title, ha="center", va="center",
            fontsize=11.5, color=INK, fontweight="bold")
    ax.text(0, -1.30, desc, ha="center", va="center",
            fontsize=8.6, color=GRAY, linespacing=1.6)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "topology_legend.png"
    plt.rcParams["font.family"] = [_CJK_NAME, "DejaVu Sans"]

    fig, axes = plt.subplots(4, 2, figsize=(11.0, 14.2), dpi=300)
    fig.patch.set_facecolor("#F5F6F8")
    fig.suptitle("DECOIL 7 种 ecDNA 拓扑 —— 全部为闭合环",
                 fontsize=15, color=INK, fontweight="bold", y=0.985)
    for ax in axes.flat:
        ax.set_facecolor("white")
    plt.subplots_adjust(left=0.03, right=0.97, top=0.955, bottom=0.015,
                        wspace=0.06, hspace=0.10)

    # ---- 1. simple_circle: uniform ring + one narrow white slit ---------
    a = axes[0, 0]
    # hairline 1.5-deg slit centered at 30 deg: arc from 30.75 to 389.25
    add_ring(a, [(30.75, 389.25, BLUE, 1.0)])
    sx, sy = pol(30, 1.52)
    a.text(sx, sy, "断点", ha="left", va="center", fontsize=8.6, color=GRAY)
    style_axes(a, "1  simple_circle  简单切离环",
               "一段基因组整体切下、首尾相连\n白色窄缝 = 唯一的 excision 断点")

    # ---- 2. simple_sv: one dashed orange SV segment ----------------------
    a = axes[0, 1]
    add_ring(a, [(0, 360, BLUE, 1.0)])
    dashed_arc(a, 240, 300)
    style_axes(a, "2  simple_sv  单一简单 SV",
               "环上某段带一个 SV(小缺失或倒位)\n橙虚线弧 = SV 所在段")

    # ---- 3. mixed_sv: deletion + inversion ------------------------------
    a = axes[1, 0]
    # Blue ring = reference orientation (plain, no marker). The green
    # inversion arc itself tapers into a single arrowhead pointing
    # counter-ringwise (SnapGene feature convention): ONE tip, shape-
    # native, no overlaid markers, nothing to collide with.
    add_ring(a, [(0, 360, BLUE, 1.0)])
    arrowed_arc(a, 120, 210, GREEN, ccw=True, tip_deg=12)  # inversion arc
    dashed_arc(a, 240, 300)                                # deletion arc
    style_axes(a, "3  mixed_sv  缺失 + 倒位",
               "橙虚线 = 缺失段;绿弧 = 倒位段\n绿弧末端箭头示其方向与环相反(倒位)")

    # ---- 4. multi_region_intra_chr: shades of blue, no slits ------------
    a = axes[1, 1]
    add_ring(a, [(30, 150, BLUE, 1.0), (150, 270, BLUE_L, 1.0),
                 (270, 390, BLUE_LL, 1.0)])
    style_axes(a, "4  multi_region_intra_chr  同染色体多区域",
               "同一染色体 ≥2 个不相邻区域拼接\n深浅蓝区分区域;颜色交界 = 拼接点")

    # ---- 5. multi_region_inter_chr: one color per chromosome ------------
    a = axes[2, 0]
    add_ring(a, [(0, 120, BLUE, 1.0), (120, 240, ROSE, 1.0),
                 (240, 360, GREEN, 1.0)])
    # chromosome labels inside the ring cavity (one per arc, at its mid-angle)
    for ang, name, col in ((60, "chrA", BLUE), (180, "chrB", ROSE),
                           (300, "chrC", GREEN)):
        lx, ly = pol(ang, 0.55)
        a.text(lx, ly, name, ha="center", va="center",
               fontsize=9.0, color=col)
    style_axes(a, "5  multi_region_inter_chr  跨染色体",
               "不同染色体的片段拼接成环\n每种颜色 = 一种来源染色体")

    # ---- 6. simple_duplications: A / A' ---------------------------------
    a = axes[2, 1]
    add_ring(a, [(0, 90, BLUE, 1.0), (90, 170, RED, 1.0),
                 (170, 250, RED, 0.5), (250, 360, BLUE_L, 1.0)])
    ax_, ay_ = pol(130, 1.42)
    a.text(ax_, ay_, "A", ha="center", va="center", fontsize=9.5,
           color=RED, fontweight="bold")
    bx, by = pol(210, 1.42)
    a.text(bx, by, "A′", ha="center", va="center", fontsize=9.5,
           color=RED, fontweight="bold")
    style_axes(a, "6  simple_duplications  串联重复",
               "同一片段在环内出现 ≥2 次\n深红 / 浅红 = 两次拷贝(A 与 A′ 同源)")

    # ---- 7. foldbacks: hairpin of two same-color arcs at a fold point ----
    a = axes[3, 0]
    # Foldback = the circle re-enters the SAME segment in reverse. Genome-
    # coordinate arrows of the two copies CONVERGE at the fold point:
    #   copy 1 (B)  runs forward  -> its 3' arrow leads INTO the fold
    #   copy 2 (BR) is reversed   -> its coordinate arrow points BACK
    #                                into the fold as well
    # Blue base covers ONLY the rest of the circle (138 -> 42 via 180/0):
    # the foldback sector belongs to the orange copies — no phantom blue
    # may show between the converging tips. Solid orange dot at the apex
    # = the hinge; deep/light orange = forward vs reversed copy.
    add_ring(a, [(138, 402, BLUE, 1.0)])                                # rest
    # Block-arrow style: long full-width body + SHORT head (8 deg) — the
    # inverted-repeat glyph (▶◀). Voids shrink to a crisp notch at the apex.
    arrowed_arc(a, 42, 90, ORANGE, alpha=1.0, ccw=True, tip_deg=8)     # copy 1
    arrowed_arc(a, 90, 138, ORANGE, alpha=0.55, ccw=False, tip_deg=8)  # copy 2, tip at fold
    a.add_patch(Circle(pol(90, R - BAND / 2.0), 0.042, facecolor=ORANGE,
                       edgecolor="none", zorder=6))
    fx, fy = pol(90, 1.18)
    a.text(fx, fy, "foldback", ha="center", va="center",
           fontsize=8.6, color=ORANGE)
    style_axes(a, "7  foldbacks  折回连接",
               "两段橙弧 = 同一片段的两份拷贝(深=正向,浅=反向)\n尖端相向汇聚于折回点(橙点):第二拷贝方向相反")

    # ---- 8. graphic-language summary -------------------------------------
    a = axes[3, 1]
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.axis("off")
    a.set_facecolor("white")
    lines = [
        ("全部为闭合环(ecDNA 本身即环状分子)", INK),
        ("白色窄缝 = excision 断点(仅 1 号需要)", GRAY),
        ("橙虚线弧 = SV 所在段(缺失)", ORANGE),
        ("绿弧箭头 = 倒位段(弧端箭头示方向翻转)", GREEN_D),
        ("深浅红 A/A′ = 同一片段的重复拷贝", RED),
        ("每种颜色 = 一种来源染色体(4/5 号)", BLUE),
        ("两橙弧尖端相向 = foldback 折回重复", ORANGE),
    ]
    a.text(0.5, 0.93, "图形语言", ha="center", va="center",
           fontsize=11.5, color=INK, fontweight="bold")
    for i, (txt, col) in enumerate(lines):
        a.text(0.5, 0.80 - i * 0.105, txt, ha="center", va="center",
               fontsize=9.2, color=col)

    fig.savefig(out, facecolor=fig.get_facecolor())
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
