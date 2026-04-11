import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
def parse_te_attributes_fast(gtf_path: str) -> pd.DataFrame:
    """
    快速解析 GTF 第9列 attributes，矢量化提取 subfamily/family/class 并去重。
    适合内存足够、文件不是极大场景。
    返回 DataFrame(columns=['subfamily','family','class'])
    """
    # 只读取第9列（index 8）
    attrs = pd.read_csv(
        gtf_path,
        sep="\t",
        header=None,
        usecols=[8],
        names=["attr"],
        comment="#",
        quoting=3,           # 避免引号解析问题
        dtype={"attr": "string"},
        compression="infer"
    )["attr"].fillna("")

    # 矢量化提取
    sub = attrs.str.extract(r'gene_id\s+"([^"]+)"')[0]
    sub.fillna(attrs.str.extract(r'transcript_id\s+"([^"]+)"')[0], inplace=True)
    sub = sub.fillna("NA")

    fam = attrs.str.extract(r'family_id\s+"([^"]+)"')[0].fillna("NA")
    # cls = attrs.str.extract(r'class_id\s+"([^"]+)"')[0].fillna("NA")

    out = pd.DataFrame({
        "subfamily": sub,
        "family": fam,
    })

    # 去重并重置索引
    out = out.drop_duplicates().reset_index(drop=True)
    return out

def parse_te_attributes_chunks(gtf_path: str, chunksize: int = 200_000) -> pd.DataFrame:
    """
    分块解析属性字段，适合内存受限或非常大的 GTF。
    使用集合来累积唯一 (sub,family,class) 三元组，最后转 DataFrame。
    """
    unique = set()
    reader = pd.read_csv(
        gtf_path,
        sep="\t",
        header=None,
        usecols=[8],
        names=["attr"],
        comment="#",
        quoting=3,
        dtype={"attr": "string"},
        compression="infer",
        chunksize=chunksize
    )

    for chunk in reader:
        attrs = chunk["attr"].fillna("")
        sub = attrs.str.extract(r'gene_id\s+"([^"]+)"')[0]
        sub.fillna(attrs.str.extract(r'transcript_id\s+"([^"]+)"')[0], inplace=True)
        sub = sub.fillna("NA")

        fam = attrs.str.extract(r'family_id\s+"([^"]+)"')[0].fillna("NA")
        cls = attrs.str.extract(r'class_id\s+"([^"]+)"')[0].fillna("NA")

        # 将每一行加入集合以去重
        for s, f, c in zip(sub, fam, cls):
            unique.add((s if pd.notna(s) else "NA",
                        f if pd.notna(f) else "NA",
                        c if pd.notna(c) else "NA"))

    out = pd.DataFrame(list(unique), columns=["subfamily", "family", "class"])
    out = out.sort_values(["family", "subfamily"]).reset_index(drop=True)
    return out

def TEfamily(
    DEG_file: str,
    outPrefix: str,
    logfc_col: str = "log2FoldChange",
    p_col: str = "pvalue",
    figsize:tuple = (10,6)
):
    df_TE = pd.read_csv(DEG_file, sep="\t")

    # 分拆 subfamily, family, class
    parts = df_TE.iloc[:, 0].str.extract(r'^([^:]+):([^:]+):([^:]+)$')
    parts.columns = ["subfamily", "family", "class"]
    new_df = pd.concat([df_TE, parts], axis=1)

    # 筛选显著
    significant_df = new_df[
        (abs(new_df[logfc_col]) > 0.58) &
        (-np.log10(new_df[p_col]) > 1.3)
    ].copy()

    if significant_df.empty:
        raise ValueError("The filtered DataFrame of significant TEs is empty.")

    # Up / Down
    significant_df["regulation"] = np.where(significant_df[logfc_col] > 0, "Up", "Down")

    result = (
        significant_df.groupby(["family", "regulation"])
        .size()
        .reset_index(name="count")
    )

    # 保存
    result_pivot = result.pivot(index="family", columns="regulation", values="count").fillna(0).astype(int)
    result_pivot.to_csv(f"{outPrefix}_TEfamily.tsv", sep="\t", header=True)

    # ---------------- 绘图 ----------------
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figsize)

    up_df   = result[result["regulation"] == "Up"].sort_values("count", ascending=False)
    down_df = result[result["regulation"] == "Down"].sort_values("count", ascending=False)

    up_families   = up_df["family"].tolist()
    down_families = down_df["family"].tolist()

    up_colors   = sns.color_palette("dark", len(up_families))
    down_colors = sns.color_palette("pastel", len(down_families))

    # 绘制 Up（左侧）
    for i, (fam, cnt) in enumerate(up_df[["family", "count"]].values):
        ax.bar(i, cnt, color=up_colors[i])
        ax.text(i, cnt + 0.1, str(cnt), ha="center", va="bottom", fontsize=9)

    # 分隔
    offset = len(up_families) + 2

    # 绘制 Down（右侧）
    for i, (fam, cnt) in enumerate(down_df[["family", "count"]].values):
        ax.bar(offset + i, cnt, color=down_colors[i])
        ax.text(offset + i, cnt + 0.1, str(cnt), ha="center", va="bottom", fontsize=9)

    xticks = list(range(len(up_families))) + list(range(offset, offset + len(down_families)))
    xlabels = up_families + down_families
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha="right",fontsize=13)

    ax.set_ylabel("Count", fontsize=15)

    # 最高值
    y_max = max(result["count"])

    # 小标题：Up 和 Down
    ax.hlines(y=y_max*1.1, xmin=0, xmax=len(up_families)-1, color="black", linewidth=1.5)
    ax.text((len(up_families)-1)/2, y_max*1.13, "Upregulated", ha="center", fontsize=15)

    ax.hlines(y=y_max*1.1, xmin=offset, xmax=offset+len(down_families)-1, color="black", linewidth=1.5)
    ax.text(offset + (len(down_families)-1)/2, y_max*1.13, "Downregulated", ha="center", fontsize=15)

    ax.set_ylim(0, y_max * 1.25)

    plt.tight_layout()
    plt.savefig(f"{outPrefix}_TEfamily.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    # TEfamily("/disk5/luosg/Totipotent20251031/Totipotent/result/ci8CLC/DESeq2/TEcount_TE.tsv",
    #          "/disk5/luosg/Totipotent20251031/Totipotent/result/ci8CLC/DESeq2/plot/TEfamily/ci8CLC")
    # TEfamily("/disk5/luosg/Totipotent20251031/Totipotent/result/ciTotiSC/DESeq2/TEcount_TE.tsv",
    #         "/disk5/luosg/Totipotent20251031/Totipotent/result/ciTotiSC/DESeq2/plot/TEfamily/ciTotiSC")
    TEfamily("/disk5/luosg/Totipotent20251031/Totipotent/result/hTBLC/DESeq2/TEcount_TE.tsv",
            "/disk5/luosg/Totipotent20251031/Totipotent/result/hTBLC/DESeq2/plot/TEfamily/hTBLC",
            figsize=(12,7.2))    
    # TEfamily("/disk5/luosg/Totipotent20251031/Totipotent/result/TLSC/DESeq2/TEcount_TE.tsv",
    #         "/disk5/luosg/Totipotent20251031/Totipotent/result/TLSC/DESeq2/plot/TEfamily/TLSC")   