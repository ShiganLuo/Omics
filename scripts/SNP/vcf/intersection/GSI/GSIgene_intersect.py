from typing import Dict, Set
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
from matplotlib.patches import Patch
# ================================
# 1. Functional mutation 定义
# ================================
#     "synonymous SNV"
functional_terms = {
    "nonsynonymous SNV",
    "stopgain",
    "stoploss",
    "frameshift insertion",
    "frameshift deletion",
    "frameshift block substitution"
    "nonframeshift insertion",
    "nonframeshift deletion",
    "nonframeshift block substitution"
}

# ================================
# 2. 判断突变是否 functional
# ================================
def classify_mutation(row,func_region_col:str="Func.refGene",func_type_col:str="ExonicFunc.refGene"):
    func = row[func_region_col]
    exonic = row[func_type_col]

    # 非编码区 → 一律非功能
    if func not in ["exonic", "splicing"]:
        return "nonfunctional"

    # 编码区但未命中 functional category
    if exonic not in functional_terms:
        return "nonfunctional"

    return "functional"

# ================================
# 3. 主函数
# ================================
def analyze_gsi_mutation(
        variant_file:str,
        GSI_genes:list,
        gene_col:str = "GeneName.symbol",
        func_type_col:str = "ExonicFunc.refGene",
        func_region_col:str = "Func.refGene",
        outdir:str="gsi_mutation",
        sep:str = ","):
    df = pd.read_csv(variant_file,sep=sep)
    df["is_GSI"] = df[gene_col].isin(GSI_genes)

    df["impact"] = df.apply(classify_mutation, axis=1,func_type_col=func_type_col,func_region_col=func_region_col)

    # ================================
    # 图 1：GSI vs background 的 functional proportion
    # ================================
    prop = df.groupby("is_GSI")["impact"].apply(
        lambda x: (x == "functional").mean()
    ).reset_index()

    prop["group"] = prop["is_GSI"].map({True: "GSI genes", False: "Background genes"})

    plt.figure(figsize=(6, 5))
    sns.barplot(data=prop, x="group", y="impact")
    plt.ylabel("Functional mutation proportion")
    plt.title("Functional mutation rate in GSI vs Background genes")
    plt.tight_layout()
    plt.savefig(f"{outdir}/functional_proportion.png", dpi=300)
    plt.close()

    # ================================
    # 图 2：GSI 基因内部 mutation 注释分布
    # ================================
    gsi_df = df[df["is_GSI"] & (df[func_type_col] != '.')]
    gsi_df.to_csv(f"{outdir}/gsi_exonic.csv")

    plt.figure(figsize=(8, 5))
    sns.countplot(data=gsi_df, y=f"{func_type_col}", order=gsi_df[func_type_col].value_counts().index)
    plt.title("Mutation category distribution in GSI genes")
    plt.xlabel("Count")
    plt.tight_layout()
    plt.savefig(f"{outdir}/GSI_mutation_category.png", dpi=300)
    plt.close()

    # ================================
    # 图 3：每个 GSI gene 的 functional mutation 数
    # ================================
    gsi_func = gsi_df[gsi_df["impact"] == "functional"]

    count_per_gene = gsi_func[gene_col].value_counts().reset_index()
    count_per_gene = count_per_gene.head(10)
    count_per_gene.columns = ["Gene", "Functional_mutations"]

    plt.figure(figsize=(8, max(4, len(count_per_gene)/3)))
    sns.barplot(data=count_per_gene, y="Gene", x="Functional_mutations")
    plt.title("Functional mutations per GSI gene")
    plt.tight_layout()
    plt.savefig(f"{outdir}/GSI_functional_per_gene.png", dpi=300)
    plt.close()

    return prop, count_per_gene, gsi_df

def read_gmt_file(gmt_file_path: str) -> Dict[str, Set[str]]:
    """
    读取 GMT 文件并将其内容解析为基因集字典。

    Args:
        gmt_file_path: GMT 文件的完整路径。

    Returns:
        一个字典，键为基因集名称 (Gene Set Name)，值为该基因集包含的基因集合 (Set of Genes)。
    """
    gene_sets: Dict[str, Set[str]] = {}
    
    try:
        with open(gmt_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    gene_set_name = parts[0]
                    genes = set(parts[2:])
                    
                    genes.discard('') 
                    
                    if genes:
                        gene_sets[gene_set_name] = genes
                        
        print(f"成功读取 {len(gene_sets)} 个基因集。")
        return gene_sets
        
    except FileNotFoundError:
        print(f"错误：文件未找到，请检查路径: {gmt_file_path}")
        return {}
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return {}


def GSI_gene(
        gmt_file:str
):
    all_gene_sets = read_gmt_file(gmt_file)
    total_gmt_genes = set()
    for gene_set in all_gene_sets.values():
        total_gmt_genes.update(gene_set)
    return total_gmt_genes

def plot_mutation_heatmap(files: Dict[str, str],
                          group_map: Dict[str, str],
                          out_file: str = "mutation_heatmap.png",
                          func_type_col:str = "ExonicFunc.refGene",
                          percentage: bool = False):
    """
    使用 sns.heatmap 绘制多个数据集外显子突变功能类型热图
    支持 dataset 按 group 排序、group row_colors、y 轴标签完全显示、列标签斜放、colorbar 居中

    Args:
        files: dict, {dataset_name: file_path}
        group_map: dict, {dataset_name: group_name}
        out_file: str, 输出图片路径
        percentage: bool, 是否绘制百分比热图

    Returns:
        heatmap_data: pd.DataFrame, 行=dataset，列=ExonicFunc.refGene，值=数量或百分比
    """
    # -------------------------
    # 读取数据并合并
    # -------------------------
    df_list = []
    for name, f in files.items():
        df = pd.read_csv(f)
        df["dataset"] = name
        df["group"] = group_map.get(name, "Unknown")
        df_list.append(df)
    all_df = pd.concat(df_list, ignore_index=True)

    # -------------------------
    # 统计每个 dataset 每种突变类型数量
    # -------------------------
    heatmap_data = all_df.groupby(["dataset", func_type_col]).size().unstack(fill_value=0)

    if percentage:
        plot_data = heatmap_data.div(heatmap_data.sum(axis=1), axis=0) * 100
        fmt = ".1f"
        cbar_label = "% of mutations"
    else:
        plot_data = heatmap_data
        fmt = "d"
        cbar_label = "Mutation count"

    # -------------------------
    # dataset 按 group 排序
    # -------------------------
    dataset_groups = pd.Series([group_map.get(ds, "Unknown") for ds in plot_data.index], index=plot_data.index)
    sorted_index = dataset_groups.sort_values().index
    plot_data = plot_data.loc[sorted_index]
    dataset_groups = dataset_groups.loc[sorted_index]

    # -------------------------
    # group 对应颜色
    # -------------------------
    unique_groups = dataset_groups.unique()
    palette = sns.color_palette("Set2", len(unique_groups))
    group_colors_map = dict(zip(unique_groups, palette))
    row_colors = dataset_groups.map(group_colors_map)

    # -------------------------
    # 绘制热图
    # -------------------------
    sns.set(style="white")
    fig, ax = plt.subplots(figsize=(max(12, len(plot_data.columns)*0.8), max(6, len(plot_data)*0.5)))

    sns.heatmap(plot_data, ax=ax, cmap="Reds", annot=True, fmt=fmt,
                cbar_kws={'label': cbar_label, 'orientation':'vertical'},
                linewidths=0.5, linecolor='gray')

    # -------------------------
    # 添加 row_colors 左侧的 group 颜色条
    # -------------------------
    # for y, color in enumerate(row_colors):
    #     ax.add_patch(plt.Rectangle((-0.5, y), width=-0.2, height=1, color=color, transform=ax.transData, clip_on=False))

    # -------------------------
    # 列标签斜放
    # -------------------------
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # y 轴标签显示
    ax.set_yticklabels(plot_data.index, rotation=0, ha='right')

    # -------------------------
    # 添加 group 图例
    # -------------------------
    
    # handles = [Patch(facecolor=group_colors_map[g], label=g) for g in unique_groups]
    # ax.legend(handles=handles, title="Group", bbox_to_anchor=(-0.2, 0), loc='upper left',frameon=False)

    plt.title("")
    plt.xlabel("")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.show()

    return heatmap_data

if __name__ == "__main__":
    # human_files = ["/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/table/toti_only.vcf.GRCh38_multianno_with_names.csv",
    #                "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/hTBLC/table/toti_only.vcf.GRCh38_multianno_with_names.csv"]
    # mouse_files = ["/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/TLSC/table/toti_only.vcf.GRCm39_multianno_with_names.csv",
    #                "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ciTotiSC/table/toti_only.vcf.GRCm39_multianno_with_names.csv"]
    # GSI_human = GSI_gene("/disk5/luosg/Totipotent20251031/data/geneset/GSI_human.gmt")
    # GSI_mouse = GSI_gene("/disk5/luosg/Totipotent20251031/data/geneset/GSI_mouse.gmt")
    # outdir = Path("/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate")
    # for infile in human_files:
    #     cell = Path(infile).parent.parent.name
    #     cell_outdir = outdir / f"{cell}/functional"
    #     cell_outdir.mkdir(parents=True,exist_ok=True)
    #     prop, count_per_gene, gsi_mut = analyze_gsi_mutation(
    #         variant_file=infile,
    #         GSI_genes=GSI_human,
    #         outdir=str(cell_outdir)
    #     )
    # for infile in mouse_files:
    #     cell = Path(infile).parent.parent.name
    #     cell_outdir = outdir / f"{cell}/functional"
    #     cell_outdir.mkdir(parents=True,exist_ok=True)
    #     prop, count_per_gene, gsi_mut = analyze_gsi_mutation(
    #         variant_file=infile,
    #         GSI_genes=GSI_mouse,
    #         outdir=str(cell_outdir)
    #     )

    ### heapmap
    # files = {
    #     "ci8CLC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/functional/gsi_exonic.csv",
    #     "hTBLC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/hTBLC/functional/gsi_exonic.csv",
    #     "TLSC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/TLSC/functional/gsi_exonic.csv",
    #     "ciTotiSC": "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ciTotiSC/functional/gsi_exonic.csv"
    # }
    # group_map = {
    #     "ci8CLC": "Human",
    #     "hTBLC": "Human",
    #     "TLSC": "Mouse",
    #     "ciTotiSC": "Mouse"
    # }
    # heatmap_data = plot_mutation_heatmap(
    #     files=files,
    #     group_map=group_map,
    #     out_file="/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/GSI_heatmap.png",
    #     percentage=False  # True 绘制百分比热图
    # )
    GSI_human = GSI_gene("/disk5/luosg/Totipotent20251031/data/geneset/GSI_human.gmt")
    GSI_mouse = GSI_gene("/disk5/luosg/Totipotent20251031/data/geneset/GSI_mouse.gmt")
    # mouse = ["TLSC","ciTotiSC"]
    # human = ["ci8CLC","hTBLC"]
    # indir = "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect"
    # infiles = glob.glob(f"{indir}/**/*_multianno_with_names.csv",recursive=True)
    # target_names = ["pluripotency","totipotency"]
    # for infile in infiles:
    #     infile = Path(infile)
    #     outdir = infile.parent.parent / "functional"
    #     outdir.mkdir(parents=True,exist_ok=True)
    #     print(infile)
    #     print(outdir)
    #     if infile.parent.parent.name in target_names:
    #         if infile.parent.parent.parent.name in human:
    #             prop, count_per_gene, gsi_mut = analyze_gsi_mutation(
    #                 variant_file=infile,
    #                 GSI_genes=GSI_human,
    #                 outdir=str(outdir)
    #             )
    #             print("human\n")
    #         elif infile.parent.parent.parent.name in mouse:
    #             print("mouse\n")
    #             prop, count_per_gene, gsi_mut = analyze_gsi_mutation(
    #                 variant_file=infile,
    #                 GSI_genes=GSI_mouse,
    #                 outdir=str(outdir)
    #             )
    #     else:
    #         print("skip:\t" + str(infile))
    # vcf_files = {
    #     "SL_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777769.vcf.gz",
    #     "SL_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777770.vcf.gz",
    #     "2iL_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777771.vcf.gz",
    #     "2iL_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777772.vcf.gz",
    #     "A2iL_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777773.vcf.gz",
    #     "A2iL_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777774.vcf.gz",
    #     "LCDM_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777775.vcf.gz",
    #     "LCDM_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777776.vcf.gz",
    #     "2iL-F_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777777.vcf.gz",
    #     "2iL-F_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777778.vcf.gz"
    # }
    # sample_groups = {
    #     "SL_3": "SL",
    #     "SL_5": "SL",
    #     "2iL_3": "2iL",
    #     "2iL_5": "2iL",
    #     "A2iL_3": "A2iL",
    #     "A2iL_5": "A2iL",
    #     "LCDM_3": "LCDM",
    #     "LCDM_5": "LCDM",
    #     "2iL-F_3": "2iL-F",
    #     "2iL-F_5": "2iL-F"
    # }
    # basedir = Path("/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/gsi")
    # for sample,infile in vcf_files.items():
    #     infile = Path(infile)
    #     group = sample_groups[sample]
    #     outdir = basedir / group
    #     outdir.mkdir(parents=True,exist_ok=True)
    #     prop, count_per_gene, gsi_mut = analyze_gsi_mutation(
    #         variant_file=infile,
    #         GSI_genes=GSI_mouse,
    #         outdir=str(outdir),
    #         sep = "\t",
    #         gene_col = "GeneName",
    #         func_type_col = "ExonicFunc",
    #         func_region_col = "Func"
    #     )
    infiles = {
        "SL": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/gsi/2iL/gsi_exonic.csv",
        "LCDM": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/gsi/LCDM/gsi_exonic.csv",
        "2iL": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/gsi/2iL/gsi_exonic.csv",
        "2iL-F": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/gsi/2iL-F/gsi_exonic.csv",
        "A2iL": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/gsi/A2iL/gsi_exonic.csv"
    }
    group_map = {
        "SL": "Human",
        "LCDM": "Human",
        "2iL": "Mouse",
        "2iL-F": "Mouse",
        "A2iL": "Mouse",
    }
    plot_mutation_heatmap(infiles,group_map,func_type_col="ExonicFunc",out_file = "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/gsi/heatmap.png")
