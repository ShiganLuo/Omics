#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Somatic mutation positive selection analysis for multiple samples
Tailored for somatic mutations (not QTL)
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def somatic_selection_pipeline(
    mutation_files,           # list of CSV files, each for one sample
    output_dir="somatic_analysis",
    vaf_threshold=0.3,
    cadd_threshold=20,
    perform_dn_ds=True
):
    """
    Multi-sample somatic mutation analysis for positive selection signals.
    
    Parameters
    ----------
    mutation_files : list of str
        CSV files containing somatic mutations. Required columns:
        Sample, Gene, Chr, Pos, Ref, Alt, VAF, MutationType, CADD
    output_dir : str
        Directory to save results and plots.
    vaf_threshold : float
        Minimum VAF to consider for candidate mutations.
    cadd_threshold : float
        Minimum CADD score to consider for candidate mutations.
    perform_dn_ds : bool
        Whether to perform dN/dS analysis using dndscv.
    """
    
    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------
    # 1. Load and merge all mutation files
    # -----------------------------
    all_mutations = []
    for f in mutation_files:
        df = pd.read_csv(f)
        all_mutations.append(df)
    mutations = pd.concat(all_mutations, ignore_index=True)
    print(f"Total mutations across {len(mutation_files)} files: {mutations.shape[0]}")

    # -----------------------------
    # 2. VAF distribution
    # -----------------------------
    plt.figure(figsize=(8,5))
    sns.histplot(mutations['VAF'], bins=50, color='skyblue')
    plt.xlabel('Variant Allele Frequency')
    plt.ylabel('Count')
    plt.title('VAF Distribution (All Samples)')
    plt.savefig(os.path.join(output_dir, "VAF_distribution_all_samples.png"))
    plt.close()

    # -----------------------------
    # 3. Hotspot gene analysis
    # -----------------------------
    gene_counts = mutations['Gene'].value_counts()
    hotspot_genes = gene_counts[gene_counts > 1]

    plt.figure(figsize=(12,5))
    sns.barplot(x=hotspot_genes.index, y=hotspot_genes.values, color='salmon')
    plt.xticks(rotation=90)
    plt.xlabel('Gene')
    plt.ylabel('Mutation count')
    plt.title('Mutation Hotspot Genes (All Samples)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "hotspot_genes_all_samples.png"))
    plt.close()

    # -----------------------------
    # 4. High-impact mutation filtering
    # -----------------------------
    high_impact = mutations[mutations['CADD'] > cadd_threshold]
    print(f"High-impact mutations (CADD>{cadd_threshold}): {high_impact.shape[0]}")

    # -----------------------------
    # 5. dN/dS analysis (optional)
    # -----------------------------
    selected_genes_dn_ds = mutations['Gene'].unique().tolist()
    if perform_dn_ds:
        try:
            from dndscv import dndscv
            dds = dndscv(mutations[['Sample','Gene','Chr','Pos','Ref','Alt']])
            results_dn_ds = dds.results_genes
            results_dn_ds.to_csv(os.path.join(output_dir, "dnds_results_all_samples.csv"), index=False)
            selected_genes_dn_ds = results_dn_ds[results_dn_ds['dNdS>1']]['gene_name'].tolist()
            print("dN/dS analysis completed")
        except ModuleNotFoundError:
            print("dndscv not installed. Skipping dN/dS analysis.")

    # -----------------------------
    # 6. Candidate mutations (综合筛选)
    # -----------------------------
    candidate_mutations = mutations[
        (mutations['VAF'] > vaf_threshold) &
        (mutations['CADD'] > cadd_threshold) &
        (mutations['Gene'].isin(selected_genes_dn_ds))
    ]
    candidate_mutations.to_csv(os.path.join(output_dir, "candidate_mutations_all_samples.csv"), index=False)
    print(f"Candidate mutations after filtering: {candidate_mutations.shape[0]}")

    # -----------------------------
    # 7. Candidate VAF distribution
    # -----------------------------
    plt.figure(figsize=(8,5))
    sns.histplot(candidate_mutations['VAF'], bins=30, color='green')
    plt.xlabel('VAF')
    plt.ylabel('Count')
    plt.title('Candidate Mutations VAF Distribution')
    plt.savefig(os.path.join(output_dir, "candidate_VAF_distribution_all_samples.png"))
    plt.close()

    # -----------------------------
    # 8. Candidate hotspot genes
    # -----------------------------
    candidate_gene_counts = candidate_mutations['Gene'].value_counts()
    plt.figure(figsize=(12,5))
    sns.barplot(x=candidate_gene_counts.index, y=candidate_gene_counts.values, color='purple')
    plt.xticks(rotation=90)
    plt.xlabel('Gene')
    plt.ylabel('Count')
    plt.title('Candidate Mutation Hotspot Genes')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "candidate_hotspot_genes_all_samples.png"))
    plt.close()

    print("Analysis completed. All results saved to:", output_dir)
    
    return candidate_mutations

if __name__ == "__main__":
    mutation_files = [
        "sample1_somatic.csv",
        "sample2_somatic.csv",
        "sample3_somatic.csv"
    ]

    candidates = somatic_selection_pipeline(
        mutation_files,
        output_dir="somatic_analysis_results",
        vaf_threshold=0.25,
        cadd_threshold=15,
        perform_dn_ds=True
    )
