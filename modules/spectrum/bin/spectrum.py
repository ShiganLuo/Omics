import argparse
import csv
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

def plot_stacked_bar(tri_df, output_file, show_plot=True):
    """绘制堆积柱状图"""
    tri_df.plot(kind="bar", stacked=True, figsize=(20,6))
    plt.ylabel("Mutation count")
    plt.xlabel("Trinucleotide context")
    # plt.title("Mutation Spectrum (Stacked by Group)")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    if show_plot:
        plt.show()

def plot_heatmap(tri_df, output_file, show_plot=True):
    """绘制热图"""
    plt.figure(figsize=(12,8))
    sns.heatmap(tri_df, cmap="Reds", annot=False)
    # plt.title("Mutation Spectrum Heatmap")
    plt.ylabel("Trinucleotide context")
    plt.xlabel("Group")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    if show_plot:
        plt.show()

def mutation_spectrum_analysis(vcf_files, sample_groups, ref_fasta, output_prefix, show_plot=True):
    ref_seq = load_reference_genome(ref_fasta)
    all_counts = {}
    for sample, vcf in vcf_files.items():
        print(f"Processing {sample} ...")
        all_counts[sample] = process_vcf(vcf, ref_seq)
    
    tri_df = build_group_dataframe(all_counts, sample_groups)
    tri_df.to_csv(f"{output_prefix}.csv")
    plot_stacked_bar(tri_df, f"{output_prefix}_stacked_bar.png", show_plot=show_plot)
    plot_heatmap(tri_df, f"{output_prefix}_heatmap.png", show_plot=show_plot)
    
    return tri_df

def read_mapping_tsv(file_path, key_col, value_col):
    mapping = {}
    with open(file_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if key_col not in reader.fieldnames or value_col not in reader.fieldnames:
            raise ValueError(f"Missing required columns: {key_col}, {value_col}")
        for row in reader:
            key = row[key_col].strip()
            value = row[value_col].strip()
            if key:
                mapping[key] = value
    return mapping

def parse_args():
    parser = argparse.ArgumentParser(description="Mutation spectrum analysis")
    parser.add_argument("--vcf-map", help="TSV with columns: experiment_sample_id, vcf")
    parser.add_argument("--group-map", help="TSV with columns: experiment_sample_id, group")
    parser.add_argument("--ref-fasta", required=True, help="Reference FASTA")
    parser.add_argument("--output-prefix", required=True, help="Output prefix")
    parser.add_argument("--no-show", action="store_true", help="Do not show plots")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if not args.vcf_map or not args.group_map:
        raise SystemExit("--vcf-map and --group-map are required")

    vcf_files = read_mapping_tsv(args.vcf_map, "experiment_sample_id", "vcf")
    sample_groups = read_mapping_tsv(args.group_map, "experiment_sample_id", "group")
    mutation_spectrum_analysis(
        vcf_files,
        sample_groups,
        args.ref_fasta,
        args.output_prefix,
        show_plot=not args.no_show,
    )
