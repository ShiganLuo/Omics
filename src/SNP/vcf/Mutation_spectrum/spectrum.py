import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from Bio import SeqIO

# -------- 基本突变顺序 --------
BASE_MUT_ORDER = ["C>A","C>G","C>T","T>A","T>C","T>G"]

def load_reference_genome(fasta_file):
    """加载参考基因组序列"""
    print("Loading reference genome...")
    ref_seq = {record.id: str(record.seq).upper() for record in SeqIO.parse(fasta_file, "fasta")}
    print("Reference genome loaded.")
    return ref_seq

def get_trinucleotide(ref_seq, chrom, pos, ref, alt):
    """获取突变上下文 trinucleotide 并标准化为 pyrimidine 方向"""
    seq = ref_seq.get(chrom)
    if not seq:
        return None
    pos0 = pos - 1
    if pos0 < 1 or pos0+1 >= len(seq):
        return None
    tri = seq[pos0-1:pos0+2]
    ref_base = ref.upper()
    alt_base = alt.upper()
    # pyrimidine 方向统一
    if ref_base in ["A","G"]:
        complement = str.maketrans("ACGT","TGCA")
        tri = tri.translate(complement)[::-1]
        ref_base = ref_base.translate(complement)
        alt_base = alt_base.translate(complement)
    return f"{tri[0]}[{ref_base}>{alt_base}]{tri[2]}"

def process_vcf(vcf_file, ref_seq):
    """读取 VCF 文件前6列，统计 trinucleotide 突变"""
    cols = ["CHROM","POS","ID","REF","ALT","QUAL"]
    df = pd.read_csv(vcf_file, sep="\t", comment="#", header=None, usecols=range(6), names=cols)
    df = df[df["ALT"].str.len() == 1]
    df["trinuc"] = df.apply(lambda x: get_trinucleotide(ref_seq, x["CHROM"], int(x["POS"]), x["REF"], x["ALT"]), axis=1)
    df = df.dropna(subset=["trinuc"])
    return Counter(df["trinuc"])

def build_group_dataframe(all_counts, sample_groups):
    """构建按组累加的 trinucleotide DataFrame"""
    group_counts = {}
    for sample, counts in all_counts.items():
        group = sample_groups[sample]
        if group not in group_counts:
            group_counts[group] = Counter()
        group_counts[group] += counts
    
    # 按 6 类基本突变顺序整理 trinucleotide
    tri_list = []
    for base_mut in BASE_MUT_ORDER:
        for tri in sorted(set().union(*[c.keys() for c in group_counts.values()])):
            if f"[{base_mut}]" in tri:
                tri_list.append(tri)
    
    tri_df = pd.DataFrame(index=tri_list, columns=group_counts.keys()).fillna(0)
    for group, counts in group_counts.items():
        for tri, num in counts.items():
            if tri in tri_df.index:
                tri_df.loc[tri, group] = num
    return tri_df

def plot_stacked_bar(tri_df, output_file):
    """绘制堆积柱状图"""
    tri_df.plot(kind="bar", stacked=True, figsize=(20,6))
    plt.ylabel("Mutation count")
    plt.xlabel("Trinucleotide context")
    # plt.title("Mutation Spectrum (Stacked by Group)")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()

def plot_heatmap(tri_df, output_file):
    """绘制热图"""
    plt.figure(figsize=(12,8))
    sns.heatmap(tri_df, cmap="Reds", annot=False)
    # plt.title("Mutation Spectrum Heatmap")
    plt.ylabel("Trinucleotide context")
    plt.xlabel("Group")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()

def mutation_spectrum_analysis(vcf_files, sample_groups, ref_fasta, output_prefix):
    ref_seq = load_reference_genome(ref_fasta)
    all_counts = {}
    for sample, vcf in vcf_files.items():
        print(f"Processing {sample} ...")
        all_counts[sample] = process_vcf(vcf, ref_seq)
    
    tri_df = build_group_dataframe(all_counts, sample_groups)
    tri_df.to_csv(f"{output_prefix}.csv")
    plot_stacked_bar(tri_df, f"{output_prefix}_stacked_bar.png")
    plot_heatmap(tri_df, f"{output_prefix}_heatmap.png")
    
    return tri_df

if __name__ == "__main__":
    # 样本名: VCF 文件
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
    
    # 样本分组
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
    
    ref_fasta = "/disk5/luosg/Reference/GENCODE/mouse/GRCm38/GRCm38.primary_assembly.genome.fa"
    output_prefix = "/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/spectrum/mutation_spectrum_group"

    tri_df = mutation_spectrum_analysis(vcf_files, sample_groups, ref_fasta, output_prefix)
