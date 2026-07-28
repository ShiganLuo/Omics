import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Optional, Dict
from pathlib import Path
from scipy.signal import savgol_filter
from itertools import cycle
def run_gsva(
        count_matrix:str,
        gmt_pathway:str,
        outdir:str
):
    es = gp.gsva(
    data = count_matrix,
    gene_sets = gmt_pathway,
    outdir = outdir
    )
    return es

def plot_gsea_from_csv(
    csv_path: str,
    ranked_genes: List[str],
    out_png: str,
    top_n: Optional[int] = None,
    fdr_cutoff: float = 0.05,
    label_col: str = "Term",
    color_map: Optional[Dict[str, str]] = None,
    label_font: int = 9,
    title_font: int = 14,
    fig_size: tuple = (12, 7),
    legend_bottom: float = -0.18
):
    """
    绘制 GSEA 显著通路累积 ES 曲线（多彩、平滑、底部图例），图例列数根据通路数量和图宽自动调整。
    """

    # ===== 读取 CSV =====
    df = pd.read_csv(csv_path)

    # 定位 FDR 列
    fdr_col_candidates = [c for c in df.columns if "fdr" in c.lower()]
    if not fdr_col_candidates:
        raise ValueError("CSV 中未找到 FDR 列(例如 'FDR q-val')")
    fdr_col = fdr_col_candidates[0]

    # 只保留显著通路
    df_sig = df[df[fdr_col] < fdr_cutoff].copy()
    if df_sig.empty:
        print(f"没有显著通路 (FDR < {fdr_cutoff})，不绘图。")
        return

    df_sig = df_sig.sort_values(fdr_col)
    if top_n is not None:
        df_sig = df_sig.head(top_n)

    # ===== 绘图设置 =====
    plt.figure(figsize=fig_size)
    ax = plt.gca()
    N = len(ranked_genes)
    ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.6)

    # ===== 颜色循环 =====
    if color_map is None:
        color_map = {}
        colors_cycle = cycle(plt.get_cmap("tab20").colors)
    else:
        colors_cycle = None

    # ===== 绘制每条通路 =====
    for _, row in df_sig.iterrows():
        term = str(row[label_col])
        nes = row["NES"]
        lead_genes = str(row["Lead_genes"]).split(";")

        hits = np.array([1 if g in lead_genes else 0 for g in ranked_genes])
        nh = hits.sum()
        if nh == 0:
            continue
        no = N - nh
        running_es = np.cumsum(hits / nh - (1 - hits) / no)

        # 平滑
        win = min(len(running_es) - (1 - len(running_es) % 2), 101)
        if win >= 11:
            smooth_es = savgol_filter(running_es, window_length=win, polyorder=3)
        else:
            smooth_es = running_es

        # 获取颜色
        if term not in color_map:
            color = next(colors_cycle)
            color_map[term] = color
        else:
            color = color_map[term]

        ax.plot(smooth_es, linewidth=2, color=color, label=f"{term} (NES={nes:.2f})")

    # ===== 美化 =====
    ax.set_xlabel("Ranked Genes", fontsize=label_font)
    ax.set_ylabel("Enrichment Score (ES)", fontsize=label_font)
    # ax.set_title(f"GSEA Significant Pathways (FDR < {fdr_cutoff})", fontsize=title_font)
    ax.grid(True, alpha=0.3)

    # ===== 自动计算图例列数 =====
    n_terms = len(df_sig)
    # 根据图形宽度、字体大小、通路数量自动估算列数
    approx_char_per_col = 25  # 每列大约能放多少字符
    fig_width_inch = fig_size[0]
    legend_ncol = max(1, min(n_terms, int((fig_width_inch * 5) / approx_char_per_col)))

    # 图例底部
    ax.legend(
        fontsize=label_font,
        loc="upper center",
        bbox_to_anchor=(0.5, legend_bottom),
        ncol=legend_ncol,
        frameon=False
    )

    # 自动调整 subplots bottom
    plt.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=max(0.2, -legend_bottom + 0.05))
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"图保存到: {out_png}")
    return color_map


def plot_biax_line(
        local_alignment_file: str,
        global_alignment_file: str,
        out_png: str = "alignment_comparison_line.png"
    ):
    """
    Plot alignment comparison using dual y-axis line chart.
    """

    # 读取数据
    local_df = pd.read_csv(local_alignment_file, sep="\t")
    global_df = pd.read_csv(global_alignment_file, sep="\t")

    # 合并数据
    df = pd.merge(local_df, global_df, on="AminoAcid", suffixes=("_Local", "_Global"))

    # 排序（可选）
    df = df.sort_values("Similarity_Global", ascending=False)

    # 创建画布
    fig, ax1 = plt.subplots(figsize=(10,5))

    # 创建第二个y轴
    ax2 = ax1.twinx()

    # Local 折线
    line1 = ax1.plot(
        df["AminoAcid"],
        df["Similarity_Local"],
        marker="o",
        color="steelblue",
        label="Local Alignment"
    )

    # Global 折线
    line2 = ax2.plot(
        df["AminoAcid"],
        df["Similarity_Global"],
        marker="s",
        color="orange",
        label="Global Alignment"
    )

    # 轴标签
    ax1.set_ylabel("Max Score", color="steelblue")
    ax1.tick_params(axis='y', colors="steelblue")
    ax1.spines['left'].set_color("steelblue")

    ax2.set_ylabel("Global Similarity", color="orange")
    ax2.tick_params(axis='y', colors="orange")
    ax2.spines['right'].set_color("orange")
    ax2.spines['left'].set_color("none")  # 隐藏左侧边框

    # x轴旋转
    plt.xticks(rotation=45)

    # 合并图例
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper right")

    # 保存
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
def cor_plot(
	infile: str,
	cor_col: List[Tuple[str,str]],
	outprefix:str,
	font_size: float = 14.0
):
	"""Plot correlation scatter charts for one or multiple column pairs.

	Correlation and trendline method
	-------------------------------
	For each ``(x_col, y_col)`` pair:
	1. Both columns are converted to numeric with ``errors='coerce'``.
	2. Rows with NaN in either column are removed (pairwise complete cases).
	3. Pearson correlation coefficient is computed by
	   ``plot_df[x_col].corr(plot_df[y_col])``.
	4. A first-order linear trend line is fitted via least squares using
	   ``numpy.polyfit(x, y, 1)`` and overlaid on the scatter plot.

	Parameters
	----------
	infile : str
		Input TSV/CSV file path.
	cor_col : list[tuple[str, str]]
		Column pairs to plot, e.g. [("x1", "y1"), ("x2", "y2")].
	outprefix : str
		Output file prefix; final image is "{outprefix}_cor_scatter.png".
	font_size : float, default=14.0
		Base font size for title, axes labels, ticks, and annotations.

	Returns
	-------
	None
		Save one combined PNG figure as ``{outprefix}_cor_scatter.png``.
	"""
	if not cor_col:
		raise ValueError("cor_col is empty. Please provide at least one column pair.")
	if font_size <= 0:
		raise ValueError(f"font_size must be > 0, got {font_size}")

	# Ensure Chinese labels/titles render correctly in this plotting entrypoint.
	configure_chinese_font()

	in_path = Path(infile)
	if not in_path.exists():
		raise FileNotFoundError(f"Input file not found: {in_path}")

	sep = "\t" if in_path.suffix.lower() in {".tsv", ".txt"} else ","
	df = pd.read_csv(in_path, sep=sep)

	missing_cols = [
		col_name
		for x_col, y_col in cor_col
		for col_name in (x_col, y_col)
		if col_name not in df.columns
	]
	if missing_cols:
		missing_unique = list(dict.fromkeys(missing_cols))
		raise ValueError(f"Missing required columns: {missing_unique}")

	n_plots = len(cor_col)
	ncols = 2 if n_plots > 1 else 1
	nrows = int(np.ceil(n_plots / ncols))
	fig, axes = plt.subplots(
		nrows=nrows,
		ncols=ncols,
		figsize=(9 * ncols, 7 * nrows),
		squeeze=False,
	)
	flat_axes = axes.flatten()

	for idx, (x_col, y_col) in enumerate(cor_col):
		ax = flat_axes[idx]
		plot_df = df[[x_col, y_col]].copy()
		plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
		plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors="coerce")
		plot_df = plot_df.dropna(subset=[x_col, y_col])

		if plot_df.empty:
			ax.text(0.5, 0.5, "No numeric data", ha="center", va="center", fontsize=font_size)
			ax.set_title(f"{x_col} vs {y_col}", fontsize=font_size + 2)
			ax.set_xlabel(x_col, fontsize=font_size)
			ax.set_ylabel(y_col, fontsize=font_size)
			ax.tick_params(axis="both", labelsize=max(font_size - 1, 8))
			continue

		x_vals = plot_df[x_col].to_numpy(dtype=float)
		y_vals = plot_df[y_col].to_numpy(dtype=float)

		ax.scatter(x_vals, y_vals, alpha=0.7, s=24, edgecolors="none")
		ax.set_title(f"{x_col} vs {y_col}", fontsize=font_size + 2)
		ax.set_xlabel(x_col, fontsize=font_size)
		ax.set_ylabel(y_col, fontsize=font_size)
		ax.tick_params(axis="both", labelsize=max(font_size - 1, 8))
		ax.grid(alpha=0.25, linestyle="--")

		# Add a linear trend line when x has variability.
		if len(plot_df) >= 2 and np.unique(x_vals).size >= 2:
			slope, intercept = np.polyfit(x_vals, y_vals, 1)
			x_line = np.linspace(np.min(x_vals), np.max(x_vals), 100)
			y_line = slope * x_line + intercept
			ax.plot(x_line, y_line, color="#D62728", linewidth=1.6, linestyle="-")

		if len(plot_df) >= 2:
			corr = plot_df[x_col].corr(plot_df[y_col])
			ax.text(
				0.03,
				0.97,
				f"r = {corr:.3f}",
				transform=ax.transAxes,
				va="top",
				ha="left",
				fontsize=max(font_size - 1, 8),
				bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75},
			)

	for idx in range(n_plots, len(flat_axes)):
		flat_axes[idx].axis("off")

	fig.tight_layout()
	out_path = f"{outprefix}_cor_scatter.png"
	fig.savefig(out_path, dpi=300)
	plt.close(fig)
	print(f"Saved correlation plot: {out_path}")

if __name__ == "__main__":
    # run_gsva(
    #     count_matrix= "/home/luosg/Data/genomeStability/output/result/matrix/log2tpm229.tsv",
    #     gmt_pathway= "/home/luosg/Data/genomeStability/data/final.gmt",
    #     outdir= "/home/luosg/Data/genomeStability/output/result/gsva230"
    # )
    ################ gsea preprank
    infiles = {
        "ci8CLC": "/disk5/luosg/Totipotent20251031/output/result/ci8CLC/DESeq2/TEcount_Gene_name.tsv",
        "hTBLC": "/disk5/luosg/Totipotent20251031/output/result/hTBLC/DESeq2/TEcount_Gene_name.tsv",
        "TLSC": "/disk5/luosg/Totipotent20251031/output/result/TLSC/DESeq2/TEcount_Gene_name.tsv",
        "ciTotiSC": "/disk5/luosg/Totipotent20251031/output/result/ciTotiSC/DESeq2/TEcount_Gene_name.tsv"
    }
    human_cell = ["ci8CLC","hTBLC"]
    mouse_cell = ["TLSC","ciTotiSC"]
    for cell,infile in infiles.items():
        df = pd.read_csv(infile,sep="\t",index_col=0)
        geneRank = df["log2FoldChange"]
        outdir = Path(infile).parent.parent / "gsea"
        # if cell in human_cell:
        #     gsea_results = run_gsea(geneRank,
        #                             "/disk5/luosg/Totipotent20251031/data/geneset/GSI_human.gmt",
        #                             str(outdir))
        # elif cell in mouse_cell:
        #     gsea_results = run_gsea(geneRank,
        #                 "/disk5/luosg/Totipotent20251031/data/geneset/GSI_mouse.gmt",
        #                 str(outdir))
        # else:
        #     raise ValueError("not support cell")
        report = outdir / "gseapy.gene_set.prerank.report.csv"
        rnk = outdir / "prerank_data.rnk"
        df_gsea = pd.read_csv(str(report))
        df_rnk = pd.read_csv(str(rnk),header=None,sep="\t")
        geneRnk = df_rnk[0].to_list()
        outfile = outdir / f"{Path(infile).parent.parent.name}.png"
        i = 0
        if i == 0:
            color_map = plot_gsea_from_csv(report,
                            geneRnk,outfile,
                            fig_size=(12,8),
                            title_font=15,
                            label_font=12,
                            legend_bottom= -0.22)
        else:
            plot_gsea_from_csv(report,
                            geneRnk,outfile,
                            fig_size=(12,8),
                            color_map=color_map,
                            title_font=15,
                            label_font=12,
                            legend_bottom= -0.12)
        i += 1
   
    ########## multiple curve plot
    # df_gsea = pd.read_csv("/home/luosg/Data/genomeStability/output/result/gsea/gseapy.gene_set.prerank.report.csv")
    # df_rnk = pd.read_csv("/home/luosg/Data/genomeStability/output/result/gsea/prerank_data.rnk",header=None,sep="\t")
    # geneRnk = df_rnk[0].to_list()
    # plot_gsea_from_csv("/home/luosg/Data/genomeStability/output/result/gsea/gseapy.gene_set.prerank.report.csv",
    #                    geneRnk,"/home/luosg/Data/genomeStability/output/result/gsea/all.png",5)

    
