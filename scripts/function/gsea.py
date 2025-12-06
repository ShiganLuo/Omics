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

    
