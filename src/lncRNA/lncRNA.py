import pandas as pd
from typing import List,Dict
from dataclasses import dataclass
import matplotlib.pyplot as plt
import numpy as np
@dataclass
class TSamples:
    control: List[str]
    experiment: List[str]

def normalize_by_gapdh(df: pd.DataFrame, gene_col: str, all_samples: List[str], gapdh_name:str="Gapdh"):
    """
    Function: Normalize gene expression by GAPDH expression in each sample.
    Parameters:
    - df: DataFrame containing gene expression data with gene identifiers and sample columns.
    - gene_col: Column name in df that contains gene identifiers.
    - all_samples: List of column names in df corresponding to all samples to be normalized.
    - gapdh_name: The gene name of GAPDH in the gene_col column. Default is "Gapdh".
    Returns:
    - DataFrame with normalized gene expression values.
    """
    gapdh_row = df[df[gene_col] == gapdh_name]

    if gapdh_row.empty:
        raise ValueError(f"{gapdh_name} not found in {gene_col}")
    gapdh_values = gapdh_row.loc[:, all_samples].iloc[0]
    df.loc[:, all_samples] = df.loc[:, all_samples].div(gapdh_values)

    return df
def estimate_transport(
    df: pd.DataFrame,
    df_gene_col: str,
    cytoplasm_samples:TSamples,
    nuclueus_samples:TSamples,
    cell_samples:TSamples,
    df_anno: pd.DataFrame,
    anno_gene_col: str = "gene_name",
    gene_type_col: str = "gene_type",
):
    """
    Function: Estimate RNA transport status based on gene expression in cytoplasm, nucleus and whole cell samples.
    Parameters:
    - df: DataFrame containing gene expression data with gene identifiers and sample columns.
    - df_gene_col: Column name in df that contains gene identifiers.
    - cytoplasm_samples: Dictionary where keys are condition names (e.g., "control", "experiment") and values are lists of column names in df corresponding to cytoplasm samples for each condition.
    - nuclueus_samples: Dictionary where keys are condition names (e.g., "control", "experiment") and values are lists of column names in df corresponding to nucleus samples for each condition.
    - cell_samples: Dictionary where keys are condition names (e.g., "control", "experiment") and values are lists of column names in df corresponding to whole cell samples for each condition.
    - df_anno: DataFrame containing gene annotations with gene identifiers and gene types.
    - anno_gene_col: Column name in df_anno that contains gene identifiers. Default is "gene_name".
    - gene_type_col: Column name in df_anno that contains gene types. Default is "gene_type".
    Returns:
    - DataFrame with an additional column for gene type annotation.
    """
    cytoplasm_samples_set = set(cytoplasm_samples.experiment + cytoplasm_samples.control)
    nuclueus_samples_set = set(nuclueus_samples.experiment + nuclueus_samples.control)
    cell_samples_set = set(cell_samples.experiment + cell_samples.control)
    if any([cytoplasm_samples_set.intersection(nuclueus_samples_set), cytoplasm_samples_set.intersection(cell_samples_set), nuclueus_samples_set.intersection(cell_samples_set)]):
        raise ValueError("Cytoplasm 、nucleus and cell samples should not overlap.")
    if not cytoplasm_samples_set.issubset(df.columns) or not nuclueus_samples_set.issubset(df.columns) or not cell_samples_set.issubset(df.columns):
        raise ValueError("Some samples in cytoplasm_samples, nuclueus_samples or cell_samples are not present in the dataframe columns.")
    df_anno = df_anno.drop_duplicates(subset=anno_gene_col)
    df[gene_type_col] = df[df_gene_col].map(df_anno.set_index(anno_gene_col)[gene_type_col])
    df = df[(df[gene_type_col] == "lncRNA") | (df[gene_type_col] == "protein_coding")]
    df[gene_type_col] = df[gene_type_col].map({"lncRNA": "lncRNA", "protein_coding": "mRNA"})
    # select unchanged genes in cell samples
    df = df[(df[cell_samples.experiment].mean(axis=1) > 0) & (df[cell_samples.control].mean(axis=1) > 0)]
    df["fold_change"] = df[cell_samples.experiment].mean(axis=1) / (df[cell_samples.control].mean(axis=1))
    df = df[(df["fold_change"] > 0.5) & (df["fold_change"] < 2.0)]
    all_samples = cytoplasm_samples.experiment + cytoplasm_samples.control + nuclueus_samples.experiment + nuclueus_samples.control
    df = normalize_by_gapdh(df, gene_col=df_gene_col, all_samples=all_samples, gapdh_name="Gapdh")
    df["nuclear_cytoplasmic_ratio_experiment"] = df[nuclueus_samples.experiment].mean(axis=1) / (df[cytoplasm_samples.experiment].mean(axis=1))
    df["nuclear_cytoplasmic_ratio_control"] = df[nuclueus_samples.control].mean(axis=1) / (df[cytoplasm_samples.control].mean(axis=1))
    df["ratio_kd_wt"] = df["nuclear_cytoplasmic_ratio_experiment"] / (df["nuclear_cytoplasmic_ratio_control"])

    df["export_state"] = np.select(
        [
            df["ratio_kd_wt"] < 0.5,
            (df["ratio_kd_wt"] >= 0.5) & (df["ratio_kd_wt"] < 2),
            df["ratio_kd_wt"] >= 2
        ],
        [
            "export (LC < 0.5)",
            "unchanged (0.5 ≤ LC < 2)",
            "nuclear retention (LC ≥ 2)"
        ],
        default=np.nan
    )

    return df[[df_gene_col, gene_type_col, "ratio_kd_wt", "export_state"]]


def transport_plot(
    df,
    output_path: str,
    gene_type_col: str = "gene_type",
    ratio_col: str = "ratio_kd_wt",
    state_col: str = "export_state",
):
    """
    Plot transport ratio by gene type with export state coloring.

    - Legend shows ONLY export states
    - State proportions are annotated next to corresponding scatter clusters
    - NA values removed
    - Log-scaled ratio axis
    """

    df = df.dropna(subset=[gene_type_col, ratio_col, state_col])
    df = df[df[ratio_col] > 0]

    fig, ax = plt.subplots(figsize=(7, 5))

    gene_types = df[gene_type_col].unique()
    states = df[state_col].unique()

    cmap = dict(zip(states, plt.cm.Set2.colors[:len(states)]))

    for i, gene_type in enumerate(gene_types):
        group = df[df[gene_type_col] == gene_type]
        total = len(group)

        for state in states:
            sub = group[group[state_col] == state]
            if len(sub) == 0:
                continue

            x = sub[ratio_col].to_numpy()
            x = x[np.isfinite(x) & (x > 0)]
            if len(x) == 0:
                continue

            y = np.random.normal(i, 0.08, len(x))

            ax.scatter(
                x,
                y,
                color=cmap[state],
                alpha=0.6,
                s=12,
                zorder=2
            )

            pct = len(sub) / total * 100 if total > 0 else 0
            x_pos = np.median(x)

            ax.text(
                x_pos,
                i + 0.40,
                f"{pct:.1f}%",
                color=cmap[state],
                ha="left",
                va="center",
                fontsize=8,
                zorder=10
            )

    ax.set_yticks(np.arange(len(gene_types)))
    ax.set_yticklabels(gene_types)

    ax.set_xscale("log")
    ax.set_xlabel("Normalized nuclear/cytoplasmic ratio (KD / WT)")
    ax.set_ylabel("")

    handles = [
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=cmap[s], markersize=6)
        for s in states
    ]

    ax.legend(
        handles,
        states,
        title="Export state",
        loc="center left",
        frameon=False
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
if __name__ == "__main__":
    cytoplasm_samples = TSamples(control=["GSM4260783"], experiment=["GSM4260780"])
    nuclueus_samples = TSamples(control=["GSM4260784"], experiment=["GSM4260781"])
    cell_samples = TSamples(control=["GSM4260785"], experiment=["GSM4260782"])
    df_anno = pd.read_csv("/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/geneIDAnnotation.csv",sep="\t")
    df = pd.read_csv("/data/pub/zhousha/20260417_RNAseq/output/RNAseq/featureCounts/all_fpkm.tsv", sep="\t")
    df = estimate_transport(
        df=df,
        df_gene_col="gene_name",
        cytoplasm_samples=cytoplasm_samples,
        nuclueus_samples=nuclueus_samples,
        cell_samples=cell_samples,
        df_anno=df_anno,
        anno_gene_col="gene_name",
        gene_type_col="gene_type"
    )
    output_path = "/data/pub/zhousha/20260417_RNAseq/output/RNAseq/featureCounts/lncRNA_transport.png"
    transport_plot(df, output_path=output_path, gene_type_col="gene_type", ratio_col="ratio_kd_wt", state_col="export_state")