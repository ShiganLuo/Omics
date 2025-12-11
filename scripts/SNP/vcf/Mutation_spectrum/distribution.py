import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def somatic_mutation_analysis_by_group(
    vcf_files,           # 样本名: VCF 文件字典
    sample_groups,       # 样本名:组字典，例如 {"SL_3":"SL","SL_5":"SL"}
    output_prefix,       # 输出前缀
    qual_threshold=30,   # 质量阈值
    filter_pass=True,    # 是否只保留 FILTER==PASS
    save_csv=True
):
    """
    读取多个体细胞突变文件，按分组统计功能分布。
    分析：
      1. 基因区域分布 (ExonicFunc, GeneName)
      2. 重复序列 vs 非重复序列
      3. 疾病/功能通路富集 (OMIM, GO_BP, KEGG, Reactome)
    """
    all_df = []
    for sample, vcf_file in vcf_files.items():
        df = pd.read_csv(vcf_file, sep="\t")
        # 过滤 SNV
        df = df[(df['REF'].str.len()==1) & (df['ALT'].str.len()==1)]
        # if filter_pass:
        #     df = df[df['FILTER']=='PASS']
        # df = df[df['QUAL']>qual_threshold]
        df['Sample'] = sample
        df['Group'] = sample_groups[sample]
        all_df.append(df)

    df_grouped = pd.concat(all_df, ignore_index=True)

    # -----------------------------
    # 1. 基因区域分布
    # -----------------------------
    gene_func_count = df_grouped.groupby(['Group','ExonicFunc','GeneName']).size().reset_index(name='Count')
    exonic_func_count = df_grouped.groupby(['Group','ExonicFunc']).size().reset_index(name='Count')

    plt.figure(figsize=(10,6))
    sns.barplot(data=exonic_func_count, x='ExonicFunc', y='Count', hue='Group', palette="Set2")
    plt.ylabel("Mutation count")
    plt.xlabel("Exonic Function")
    plt.title("Distribution of Mutations by Exonic Function (Grouped)")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_mutation_by_exonicfunc_grouped.png", dpi=300)
    plt.close()

    # -----------------------------
    # 2. 重复序列分布分析
    # -----------------------------
    df_grouped['RepeatStatus'] = df_grouped['Repeat'].apply(lambda x: 'Repeat' if pd.notna(x) else 'Non-Repeat')
    repeat_count = df_grouped.groupby(['Group','RepeatStatus']).size().reset_index(name='Count')

    plt.figure(figsize=(8,5))
    sns.barplot(data=repeat_count, x='RepeatStatus', y='Count', hue='Group', palette="Set2")
    plt.ylabel("Mutation count")
    plt.xlabel("Repeat status")
    plt.title("Mutations in Repeat vs Non-Repeat Regions (Grouped)")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_mutation_repeat_distribution_grouped.png", dpi=300)
    plt.close()

    # -----------------------------
    # 3. 疾病/功能通路富集分析
    # -----------------------------
    def plot_top_terms_grouped(df, col, title, filename, top_n=20):
        df_exp = df[['Group', col]].dropna()
        df_exp = df_exp.assign(Term=df_exp[col].str.split(';')).explode('Term')
        top_df = df_exp.groupby(['Group','Term']).size().reset_index(name='Count')
        top_terms = top_df.groupby('Term')['Count'].sum().sort_values(ascending=False).head(top_n).index
        plot_df = top_df[top_df['Term'].isin(top_terms)]
        plt.figure(figsize=(10,6))
        sns.barplot(data=plot_df, y='Term', x='Count', hue='Group', palette='Set2')
        plt.xlabel("Number of mutations")
        plt.ylabel(title)
        plt.title(f"Top {top_n} {title} (Grouped)")
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        plt.close()
        return plot_df

    omim_count = plot_top_terms_grouped(df_grouped, 'OMIM', 'OMIM disease', f"{output_prefix}_omim_top20_grouped.png")
    go_bp_count = plot_top_terms_grouped(df_grouped, 'GO_BP', 'GO Biological Process', f"{output_prefix}_go_bp_top20_grouped.png")
    kegg_count = plot_top_terms_grouped(df_grouped, 'KEGG_PATHWAY', 'KEGG Pathway', f"{output_prefix}_kegg_top20_grouped.png")
    reactome_count = plot_top_terms_grouped(df_grouped, 'REACTOME_PATHWAY', 'Reactome Pathway', f"{output_prefix}_reactome_top20_grouped.png")

    # -----------------------------
    # 保存统计表
    # -----------------------------
    if save_csv:
        gene_func_count.to_csv(f"{output_prefix}_gene_func_count_grouped.csv", index=False)
        repeat_count.to_csv(f"{output_prefix}_repeat_count_grouped.csv", index=False)
        omim_count.to_csv(f"{output_prefix}_omim_count_grouped.csv", index=False)
        go_bp_count.to_csv(f"{output_prefix}_go_bp_count_grouped.csv", index=False)
        kegg_count.to_csv(f"{output_prefix}_kegg_count_grouped.csv", index=False)
        reactome_count.to_csv(f"{output_prefix}_reactome_count_grouped.csv", index=False)

    return {
        "gene_func_count": gene_func_count,
        "repeat_count": repeat_count,
        "omim_count": omim_count,
        "go_bp_count": go_bp_count,
        "kegg_count": kegg_count,
        "reactome_count": reactome_count
    }

if __name__ == "__main__":
    vcf_files = {
        "SL_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777769.vcf.gz",
        "SL_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777770.vcf.gz",
        "2iL_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777771.vcf.gz",
        "2iL_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777772.vcf.gz",
        "A2iL_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777773.vcf.gz",
        "A2iL_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777774.vcf.gz",
        "LCDM_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777775.vcf.gz",
        "LCDM_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777776.vcf.gz",
        "2iL-F_3": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777777.vcf.gz",
        "2iL-F_5": "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/GSM4777778.vcf.gz"
    }

    sample_groups = {
        "SL_3": "SL",
        "SL_5": "SL",
        "2iL_3": "2iL",
        "2iL_5": "2iL",
        "A2iL_3": "A2iL",
        "A2iL_5": "A2iL",
        "LCDM_3": "LCDM",
        "LCDM_5": "LCDM",
        "2iL-F_3": "2iL-F",
        "2iL-F_5": "2iL-F"
    }

    somatic_mutation_analysis_by_group(vcf_files, sample_groups, 
                                       "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/distribution/ESC")
