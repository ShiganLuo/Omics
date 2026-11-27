import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Optional, Dict
from pathlib import Path
from scipy.signal import savgol_filter
from itertools import cycle
def run_gsva(
        count_matrix:pd.DataFrame,
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

    df = pd.read_csv(csv_path)

    fdr_col_candidates = [c for c in df.columns if "fdr" in c.lower()]
    if not fdr_col_candidates:
        raise ValueError("CSV 中未找到 FDR 列(例如 'FDR q-val')")
    fdr_col = fdr_col_candidates[0]

    df_sig = df[df[fdr_col] < fdr_cutoff].copy()
    if df_sig.empty:
        print(f"没有显著通路 (FDR < {fdr_cutoff})，不绘图。")
        return

    df_sig = df_sig.sort_values(fdr_col)
    if top_n is not None:
        df_sig = df_sig.head(top_n)

    plt.figure(figsize=fig_size)
    ax = plt.gca()
    N = len(ranked_genes)
    ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.6)

    if color_map is None:
        color_map = {}
        colors_cycle = cycle(plt.get_cmap("tab20").colors)
    else:
        colors_cycle = None

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

        win = min(len(running_es) - (1 - len(running_es) % 2), 101)
        if win >= 11:
            smooth_es = savgol_filter(running_es, window_length=win, polyorder=3)
        else:
            smooth_es = running_es

        if term not in color_map:
            color = next(colors_cycle)
            color_map[term] = color
        else:
            color = color_map[term]

        ax.plot(smooth_es, linewidth=2, color=color, label=f"{term} (NES={nes:.2f})")

    ax.set_xlabel("Ranked Genes", fontsize=label_font)
    ax.set_ylabel("Enrichment Score (ES)", fontsize=label_font)
    ax.grid(True, alpha=0.3)

    n_terms = len(df_sig)
    approx_char_per_col = 25
    fig_width_inch = fig_size[0]
    legend_ncol = max(1, min(n_terms, int((fig_width_inch * 5) / approx_char_per_col)))

    ax.legend(
        fontsize=label_font,
        loc="upper center",
        bbox_to_anchor=(0.5, legend_bottom),
        ncol=legend_ncol,
        frameon=False
    )

    plt.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=max(0.2, -legend_bottom + 0.05))
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"图保存到: {out_png}")
    return color_map

if __name__ == "__main__":
    df_count = pd.read_csv("/disk5/luosg/Totipotent20251031/PRJNA663159/count/GSE157874_All_TPM.csv.gz",index_col=0)
    run_gsva(
        count_matrix= df_count,
        gmt_pathway= "/disk5/luosg/Totipotent20251031/data/geneset/GSI_mouse.gmt",
        outdir= "/disk5/luosg/Totipotent20251031/PRJNA663159/count/gsva"
    )
    infiles = {
        "ci8CLC": "/disk5/luosg/Totipotent20251031/Totipotent/result/ci8CLC/DESeq2/TEcount_Gene_name.tsv",
        "hTBLC": "/disk5/luosg/Totipotent20251031/Totipotent/result/hTBLC/DESeq2/TEcount_Gene_name.tsv",
        "TLSC": "/disk5/luosg/Totipotent20251031/Totipotent/result/TLSC/DESeq2/TEcount_Gene_name.tsv",
        "ciTotiSC": "/disk5/luosg/Totipotent20251031/Totipotent/result/ciTotiSC/DESeq2/TEcount_Gene_name.tsv"
    }
    human_cell = ["ci8CLC","hTBLC"]
    mouse_cell = ["TLSC","ciTotiSC"]
    for cell,infile in infiles.items():
        df = pd.read_csv(infile,sep="\t",index_col=0)
        geneRank = df["log2FoldChange"]
        outdir = Path(infile).parent.parent / "gsea"
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
                            fig_size=(12.1,8),
                            title_font=16,
                            label_font=13,
                            legend_bottom= -0.12)
        else:
            plot_gsea_from_csv(report,
                            geneRnk,outfile,
                            fig_size=(12.1,8),
                            color_map=color_map,
                            title_font=16,
                            label_font=13,
                            legend_bottom= -0.12)
        i += 1
   

    
