import pandas as pd
import json
import os
import sys
from typing import Dict
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from src.annotation.gene_id2name import load_gtf_gene_map
def serum_ksr_toti(
    toti_ksr:str,
    toti_serum:str,
    ksr_serum:str,
    output_json:str,
    gtf_path:str,
    lfc_cut:float = 0.58,
    padj_cut:float = 0.05
):
    """
    Function: get seed gene list
    Parameters:
        - toti_ksr: The DESeq2 results of Totipotent sample vs KSR sample
        - toti_serum: The DESeq2 results of Totipotent sample vs serum sample
        - ksr_serum: The DESeq2 results of KSR sample sample vs serum sample
    """
    df_tr = pd.read_csv(toti_ksr,sep="\t", index_col=0)
    df_ts = pd.read_csv(toti_serum,sep="\t", index_col=0)
    df_km = pd.read_csv(ksr_serum,sep="\t", index_col=0)
    gene_name_map = load_gtf_gene_map(gtf_path)
    for df in [df_tr, df_ts, df_km]:
        df.reset_index(names="gene_id", inplace=True)
        df["gene_name"] = df["gene_id"].map(gene_name_map)
    toti_vs_ksr_up = df_tr[
        (df_tr["log2FoldChange"] >= lfc_cut) & (df_tr["padj"] < padj_cut)
    ]["gene_name"].tolist()
    toti_vs_serum_up = df_ts[
        (df_ts["log2FoldChange"] >= lfc_cut) & (df_ts["padj"] < padj_cut)
    ]["gene_name"].tolist()
    toti_specific = list(set(toti_vs_ksr_up) & set(toti_vs_serum_up))

    ksr_vs_serum_up = df_km[
        (df_km["log2FoldChange"] >= lfc_cut) & (df_km["padj"] < padj_cut)
    ]["gene_name"].tolist()
    toti_vs_ksr_not_down = df_tr[df_tr["log2FoldChange"] >= -lfc_cut][
        "gene_name"
    ].tolist()
    ksr_maintained = list(set(ksr_vs_serum_up) & set(toti_vs_ksr_not_down))

    serum_vs_ksr_up = df_km[
        (df_km["log2FoldChange"] <= -lfc_cut) & (df_km["padj"] < padj_cut)
    ]["gene_name"].tolist()
    serum_vs_toti_up = df_ts[
        (df_ts["log2FoldChange"] <= -lfc_cut) & (df_ts["padj"] < padj_cut)
    ]["gene_name"].tolist()
    serum_specific = list(set(serum_vs_ksr_up) & set(serum_vs_toti_up))

    seed_genes = list(set(toti_specific + ksr_maintained + serum_specific))
    outdir = os.path.dirname(output_json)
    os.makedirs(outdir, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(
            {
                "seed_genes": seed_genes,
                "toti_specific": toti_specific,
                "ksr_maintained": ksr_maintained,
                "serum_specific": serum_specific,
            },
            f,
            indent=4
        )
    return {
        "seed_genes": seed_genes,
        "toti_specific": toti_specific,
        "ksr_maintained": ksr_maintained,
        "serum_specific": serum_specific,
    }

def filter_TF(
    seed_genes:Dict,
    tf_list:str,
    output_json:str
):
    df_tf = pd.read_csv(tf_list, sep="\t")
    tf_genes = df_tf["Symbol"].tolist()
    for key,gene_list in seed_genes.items():
        seed_genes[key] = list(set(gene_list) & set(tf_genes))
        with open(output_json.replace(".json", f"_{key}_gene_list.txt"), "w") as f:
            f.write("\n".join(seed_genes[key]))
    outdir = os.path.dirname(output_json)
    os.makedirs(outdir, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(seed_genes, f, indent=4)
    return seed_genes
if __name__ == "__main__":
    seed_genes = serum_ksr_toti(
        toti_ksr="/data/pub/zhousha/20260207_Exome/output/RNAseq/results/DEG/TLSCS_vs_KSR/DESeq2/TEcount_Gene.tsv",
        toti_serum="/data/pub/zhousha/20260207_Exome/output/RNAseq/results/DEG/TLSCS_vs_E14/DESeq2/TEcount_Gene.tsv",
        ksr_serum="/data/pub/zhousha/20260207_Exome/output/RNAseq/results/DEG/KSR_vs_E14/DESeq2/TEcount_Gene.tsv",
        gtf_path="/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf",
        output_json="/data/pub/zhousha/20260207_Exome/output/RNAseq/results/network/seed_genes_raw.json"
    )
    seed_tfs = filter_TF(
        seed_genes=seed_genes,
        tf_list="/data/pub/zhousha/Reference/mouse/TF/Mus_musculus_TF.txt",
        output_json="/data/pub/zhousha/20260207_Exome/output/RNAseq/results/network/seed_genes_TF.json"
    )