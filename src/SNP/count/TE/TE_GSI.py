import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_gi_te_heatmap(file_path, gi_markers, out_file,
                        padj_thresh=0.05, log2fc_thresh=0.0,
                        sample_prefix="GSM", figsize=(10,8), cmap="RdBu_r"):

    df = pd.read_csv(file_path, sep="\t", header=0)
    df.rename(columns={df.columns[0]: 'TE'}, inplace=True)
    df['subfamily'] = df['TE'].apply(lambda x: str(x).split(':')[0])

    # 匹配 GI marker（不区分大小写）
    mask = df['subfamily'].str.upper().isin([m.upper() for m in gi_markers])
    df_gi = df.loc[mask].copy()

    # 样本列
    sample_cols = [c for c in df_gi.columns if c.startswith(sample_prefix)]
    df_counts = df_gi[sample_cols]
    df_log = np.log2(df_counts + 1)

    # 差异显著
    df_gi['sig'] = (df_gi['padj'] < padj_thresh) & (df_gi['log2FoldChange'] > log2fc_thresh)

    # 绘制热图
    plt.figure(figsize=figsize)
    ax = sns.heatmap(df_log, cmap=cmap, center=df_log.values.mean(),
                     cbar_kws={'label':'log2(TMM counts)'},
                     yticklabels=df_gi['TE'])

    ax.set_ylabel("GI-TE Marker")
    ax.set_xlabel("Samples")

    # 获取每行 y 坐标
    yticks = ax.get_yticks()

    # 标注显著 marker 星号
    for i, sig in enumerate(df_gi['sig']):
        if sig:
            ax.text(len(sample_cols)+0.1, yticks[i], "*", color='red', fontsize=12,
                    va='center', ha='left')

    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


if __name__ == "__main__":
    mouse_gi_markers = [
        "B1_Mus1","B1_Mur1","B1_Mm","B1_Mus2",
        "B2_Mm1a","B2_Mm1t",
        "L1Md_T","L1Md_A","L1Md_Gf","L1Md_F2",
        "MERVL","MERVL-int","RLTR10","RLTR9D"
    ]

    plot_gi_te_heatmap(
        file_path="/disk5/luosg/Totipotent20251031/output/result/TLSC/DESeq2/TEcount_TE.tsv",
        gi_markers=mouse_gi_markers,
        out_file="GI_TE_heatmap.png"
    )
