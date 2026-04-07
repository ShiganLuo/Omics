from utils.VEP_SV import read_vep_tab
from utils.LogUtil import setup_logger
import pandas as pd
from typing import List, Dict, Optional
import logging
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import numpy as np

def tab_parser(
        table_file: str,
        cosmic_file:str,
        tab_gene_col: str = "SYMBOL",
        cosmic_gene_col:str = "GENE_SYMBOL",
        sv_type_pattern:str = r"(DEL|DUP|INV|INS|TRA)",
        **kwargs
):
    """
    Function: Parse a VEP annotation table and filter for cancer genes.
    Parameters:
        - table_file: path to the VEP annotation table (tab-delimited format)
        - cosmic_file: path to the COSMIC cancer gene list (tab-delimited format)
        - tab_gene_col: column name in the VEP table that contains gene symbols (default "SYMBOL")
        - cosmic_gene_col: column name in the COSMIC file that contains gene symbols (default "GENE_SYMBOL")
    Returns:    
        - DataFrame containing only the rows from the VEP table where the gene is in the COSMIC cancer gene list.
    """
    df_tab = read_vep_tab(table_file, **kwargs)
    
    df_coismic = pd.read_csv(cosmic_file, sep="\t")

    cancer_genes = set(df_coismic[cosmic_gene_col].unique())
    cancer_genes = set(g.upper() for g in cancer_genes if isinstance(g, str))

    df_tab["is_cancer_gene"] = df_tab[tab_gene_col].apply(lambda x: x.upper() in cancer_genes)
    df_tab = df_tab[df_tab["is_cancer_gene"] == True]

    col_line_prefix = kwargs.get("col_line_prefix", "#Uploaded_variation")
    df_tab["svtype"] = df_tab[col_line_prefix].str.extract(sv_type_pattern).fillna("OTHER")

    return df_tab[[tab_gene_col, "svtype"]]

def prepare_oncoprint_data(
    tab_files:Dict[str, str],
    cosmic_file:str,
    gene_col:str = "SYMBOL",
    sample_order:Optional[List[str]] = None,
    **kwargs
):
    """
    Function: Prepare data for OncoPrint visualization by parsing multiple VEP annotation tables and creating a binary matrix of gene alterations across samples.
    """
    df_list = []
    for sample, tab_file in tab_files.items():
        df_sample = tab_parser(tab_file, cosmic_file=cosmic_file, tab_gene_col=gene_col, **kwargs)
        df_sample["sample"] = sample
        df_list.append(df_sample)
    
    df = pd.concat(df_list, ignore_index=True)

    df = (
        df.groupby([gene_col, "sample"])["svtype"]
        .apply(lambda x: ",".join(sorted(set(x))))
        .reset_index()
    )

    # pivot 成 OncoPrint 矩阵
    matrix = df.pivot(index=gene_col, columns="sample", values="svtype")
    if sample_order is not None:
        matrix = matrix.reindex(columns=sample_order)
    return matrix
        
def Deseq2_oncoprint_data(
    deseq2_files: List[str],
    condition_names: List[str],
    oncoprint_file: str,
    oncoprint_gene_col: str = "SYMBOL",
    deseq2_gene_col: str = "name",
    padj_col: str = "padj",
    padj_threshold: Optional[float] = 0.05,
    log2fc_col: str = "log2FoldChange",
    sep: str = "\t",
    logger: logging.Logger = setup_logger(__name__),
) -> pd.DataFrame:
    """
    Summarize log2 fold change values of COSMIC genes across multiple conditions.

    This function reads multiple DESeq2 result files and extracts log2 fold change
    values for genes present in an oncoprint gene list. It supports an arbitrary
    number of conditions and merges results into a single wide-format DataFrame.
    Missing values are filled with NA.

    Parameters
    ----------
    deseq2_files : List[str]
        List of DESeq2 result file paths.
    condition_names : List[str]
        Names of conditions corresponding to each DESeq2 file.
    oncoprint_file : str
        Path to oncoprint gene file.
    oncoprint_gene_col : str, default="SYMBOL"
        Column name for gene symbols in oncoprint file.
    deseq2_gene_col : str, default="name"
        Column name for gene symbols in DESeq2 files.
    padj_col : str, default="padj"
        Column name for adjusted p-value.
    padj_threshold : float, default=0.05
        Threshold for filtering significant genes.
    log2fc_col : str, default="log2FoldChange"
        Column name for log2 fold change.
    sep : str, default="\\t"
        Separator for DESeq2 files.
    logger : logging.Logger, default=setup_logger(__name__)
        Logger instance for debug and trace.

    Returns
    -------
    pd.DataFrame
        A DataFrame with genes as rows and conditions as columns,
        containing fold change values. Missing values are NA.

    Raises
    ------
    ValueError
        If input file list and condition names length mismatch.
    """

    if len(deseq2_files) != len(condition_names):
        raise ValueError("Length of deseq2_files must match condition_names")

    logger.info("Loading oncoprint gene list")
    df_oncoprint = pd.read_csv(oncoprint_file)
    gene_list = (
        df_oncoprint[oncoprint_gene_col]
        .dropna()
        .astype(str)
        .str.capitalize()
        .unique()
    )
    gene_set = set(gene_list)
    logger.info(f"Total genes in oncoprint: {len(gene_set)}")

    result_df = pd.DataFrame({oncoprint_gene_col: sorted(gene_set)})

    for file, condition in zip(deseq2_files, condition_names):
        logger.info(f"Processing condition: {condition}, file: {file}")

        df = pd.read_csv(file, sep=sep)

        # Standardize gene name
        df[deseq2_gene_col] = df[deseq2_gene_col].astype(str).str.capitalize()

        # Filter genes and padj
        df_filtered = df[
            df[deseq2_gene_col].isin(gene_set) &
            (df[padj_col] < padj_threshold if padj_threshold is not None else True)
        ]

        logger.info(
            f"{condition}: {df_filtered.shape[0]} genes passed padj < {padj_threshold}"
        )

        # Keep only necessary columns
        df_filtered = df_filtered[[deseq2_gene_col, log2fc_col]].drop_duplicates(
            subset=deseq2_gene_col
        )

        # Rename column to condition name
        df_filtered = df_filtered.rename(
            columns={
                deseq2_gene_col: oncoprint_gene_col,
                log2fc_col: condition,
            }
        )
        df_filtered[condition] = 2 ** df_filtered[condition]
        # Merge
        result_df = result_df.merge(
            df_filtered,
            on=oncoprint_gene_col,
            how="left"
        )
    logger.info("Merging completed")

    return result_df

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Tuple, Optional, Callable
from scipy.stats import chi2_contingency, fisher_exact, ttest_ind, mannwhitneyu


# =========================
# statistical test utilities
# =========================
def compute_pvalue(
    table: np.ndarray,
    method: str = "chi2"
) -> float:
    """
    Compute p-value for a 2x2 contingency table or two-sample comparison.

    Parameters
    ----------
    table : np.ndarray
        Input data:
        - For "chi2" / "fisher": shape (2, 2)
        - For "t-test" / "mannwhitney": shape (2, n)
    method : str, default="chi2"
        Statistical test method. Supported:
        - "chi2" : Chi-square test
        - "fisher" : Fisher's exact test
        - "t-test" : Independent t-test
        - "mannwhitney" : Mann-Whitney U test

    Returns
    -------
    float
        P-value.
    """
    if method == "chi2":
        _, p, _, _ = chi2_contingency(table)
    elif method == "fisher":
        _, p = fisher_exact(table)
    elif method == "t-test":
        p = ttest_ind(table[0], table[1], equal_var=False).pvalue
    elif method == "mannwhitney":
        p = mannwhitneyu(table[0], table[1], alternative="two-sided").pvalue
    else:
        raise ValueError(f"Unsupported method: {method}")
    return p


def p_to_star(p: float) -> str:
    """
    Convert p-value to significance star label.
    """
    if p < 1e-4:
        return "****"
    elif p < 1e-3:
        return "***"
    elif p < 1e-2:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def compute_significance(
    pivot: pd.DataFrame,
    group_order: Tuple[str, str],
    method: str = "chi2"
) -> Dict[str, str]:
    """
    Compute significance stars per SV type.

    Parameters
    ----------
    pivot : pandas.DataFrame
        Pivot table with SV types as index and groups as columns.
    group_order : tuple of str
        Two groups to compare.
    method : str, default="chi2"
        Statistical test method.

    Returns
    -------
    dict
        Mapping: svtype -> significance stars.
    """
    stars = {}
    g1, g2 = group_order

    for sv in pivot.index:
        if method in ("chi2", "fisher"):
            table = np.array([
                [pivot.loc[sv, g1], pivot[g1].sum() - pivot.loc[sv, g1]],
                [pivot.loc[sv, g2], pivot[g2].sum() - pivot.loc[sv, g2]],
            ])
        else:
            # fallback: treat as single values (not ideal, but generic)
            table = np.array([
                [pivot.loc[sv, g1]],
                [pivot.loc[sv, g2]],
            ])

        p = compute_pvalue(table, method)
        stars[sv] = p_to_star(p)

    return stars


# =========================
# plotting core
# =========================
def plot_comparison_broken_bar(
    df: pd.DataFrame,
    out_png: str,
    group_col: str = "group",
    svtype_col: str = "svtype",
    count_col: str = "count",
    group_order: Tuple[str, str] = ("Control", "Experiment"),
    svtype_order: Tuple[str, ...] = ("BND", "DEL", "DUP", "INS", "INV"),
    legend_map: Optional[Dict[str, str]] = None,
    figsize: Tuple[int, int] = (9, 5),
    ylabel: str = "SV count",
    dpi: int = 300,
    test_method: str = "chi2",
    do_test: bool = True,
    use_broken_axis: bool = True,
) -> None:

    if legend_map is None:
        legend_map = {g: g for g in group_order}

    pivot = (
        df.pivot(index=svtype_col, columns=group_col, values=count_col)
        .reindex(svtype_order)
        .fillna(0)
    )

    # ---------- significance ----------
    stars = {}
    sig_sv = []
    if do_test:
        # 限制：当前函数仅支持计数型检验
        if test_method not in ("chi2", "fisher"):
            raise ValueError("Only 'chi2' and 'fisher' are supported for count data")

        stars = compute_significance(pivot, group_order, test_method)
        sig_sv = [sv for sv, s in stars.items() if s != "ns"]

    x = np.arange(len(pivot.index))
    width = 0.36

    colors = {
        group_order[0]: "#4C72B0",
        group_order[1]: "#DD8452",
    }

    # =========================
    # axis setup
    # =========================
    low_max = None  # ✅ FIX: 统一初始化

    if use_broken_axis:
        if len(sig_sv) > 0:
            low_max = pivot.loc[sig_sv].values.max() * 1.15
        else:
            low_max = np.median(pivot.values) * 1.5

        global_max = pivot.values.max()
        high_min = low_max * 1.1
        high_max = global_max * 1.15

        fig, (ax_top, ax_bottom) = plt.subplots(
            2, 1, sharex=True,
            figsize=figsize,
            gridspec_kw={"height_ratios": [1, 3]},
        )

        axes = (ax_top, ax_bottom)

        for ax in axes:
            ax.bar(x - width / 2, pivot[group_order[0]], width,
                   color=colors[group_order[0]], label=legend_map[group_order[0]])
            ax.bar(x + width / 2, pivot[group_order[1]], width,
                   color=colors[group_order[1]], label=legend_map[group_order[1]])

        ax_bottom.set_ylim(0, low_max)
        ax_top.set_ylim(high_min, high_max)

        # broken marks
        d = 0.008
        ax_top.plot((-d, +d), (-d, +d), transform=ax_top.transAxes, color="black", clip_on=False)
        ax_bottom.plot((-d, +d), (1 - d, 1 + d), transform=ax_bottom.transAxes, color="black", clip_on=False)

    else:
        fig, ax = plt.subplots(figsize=figsize)
        ax.bar(x - width / 2, pivot[group_order[0]], width,
               color=colors[group_order[0]], label=legend_map[group_order[0]])
        ax.bar(x + width / 2, pivot[group_order[1]], width,
               color=colors[group_order[1]], label=legend_map[group_order[1]])

        axes = (ax,)
        ax_top = ax_bottom = ax  # 统一接口

    # =========================
    # significance annotation
    # =========================
    if do_test:
        LEG_PT = 8
        TEXT_PT = 3

        fig.canvas.draw()

        for i, sv in enumerate(pivot.index):
            if stars.get(sv, "ns") == "ns":
                continue

            y1 = pivot.loc[sv, group_order[0]]
            y2 = pivot.loc[sv, group_order[1]]
            y_base = max(y1, y2)

            if use_broken_axis:
                ax = ax_bottom if y_base <= low_max else ax_top
            else:
                ax = ax_top

            trans = ax.transData
            inv = ax.transData.inverted()

            _, y_disp = trans.transform((0, y_base))
            _, y_hat = inv.transform((0, y_disp + LEG_PT))
            _, y_text = inv.transform((0, y_disp + LEG_PT + TEXT_PT))

            x1 = x[i] - width / 2
            x2 = x[i] + width / 2

            ax.plot([x1, x1], [y1, y_hat], lw=1.2, c="black")
            ax.plot([x2, x2], [y2, y_hat], lw=1.2, c="black")
            ax.plot([x1, x2], [y_hat, y_hat], lw=1.2, c="black")

            ax.text(x[i], y_text, stars[sv],
                    ha="center", va="bottom",
                    fontsize=12, fontweight="bold")

    # =========================
    # formatting
    # =========================
    ax_bottom.set_xticks(x)
    ax_bottom.set_xticklabels(pivot.index)
    ax_bottom.set_ylabel(ylabel)

    if use_broken_axis:
        ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
        ax_top.spines["bottom"].set_visible(False)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].legend(frameon=False)

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.close()


def main():
    vcfs = {
        "PlaB06_vs_DMSO06": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB_only.vcf",
        "PlaB20_vs_DMSO20": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB_only.vcf",
    }
    # vep_sv = VEP_SV(
    #     vep_cache_dir="~/.vep",
    #     species="mus_musculus",
    #     assembly="GRCm39"
    # )
    # for sample, vcf in vcfs.items():
    #     out_tab = f"/data/pub/zhousha/Totipotent20251031/PacBio/SV/{sample}/{sample}_annotated.tab"
    #     vep_sv.annotate_sv_vep(vcf, out_tab,"tab")
    # tab_files = {
    #     "PlaB06_vs_DMSO06": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB06_vs_DMSO06_annotated.tab",
    #     "PlaB20_vs_DMSO20": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB20_vs_DMSO20_annotated.tab"
    # }
    # cosmic_file = "/data/pub/zhousha/Totipotent20251031/data/Cancer/Cosmic_CancerGeneCensusHallmarksOfCancer_v103_GRCh38.tsv"
    # oncoprint_matrix = prepare_oncoprint_data(tab_files, cosmic_file=cosmic_file,sample_order=["PlaB06_vs_DMSO06", "PlaB20_vs_DMSO20"])
    # oncoprint_matrix.to_csv("/data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint/oncoprint_matrix.csv")
    Deseq2_files = [
        "/data/pub/zhousha/Totipotent20251031/data/Pacbio/RNAseq/P6.tsv",
        "/data/pub/zhousha/Totipotent20251031/data/Pacbio/RNAseq/P20.tsv"
    ]
    condition_names = ("P6", "P20")
    result_df = Deseq2_oncoprint_data(
        deseq2_files=Deseq2_files,
        condition_names=condition_names,
        oncoprint_file="/data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint/oncoprint_matrix.csv",
        oncoprint_gene_col="SYMBOL",
        deseq2_gene_col="name",
        padj_col="padj",
        log2fc_col="log2FoldChange",
        sep="\t",
        logger=setup_logger("Deseq2_oncoprint"),
        padj_threshold = None
    )
    result_df.to_csv("/data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint/RNAseq/deseq2_oncoprint_fold_change_all.csv", index=False)
    result_df = result_df.dropna()
    plot_comparison_broken_bar(
        df=result_df.melt(id_vars="SYMBOL", var_name="group", value_name="fold_change"),
        out_png="/data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint/RNAseq/deseq2_oncoprint_fold_change_all.png",
        group_col="group",
        svtype_col="SYMBOL",
        count_col="fold_change",
        group_order=condition_names,
        svtype_order=result_df["SYMBOL"].tolist(),
        legend_map={"P6": "PlaB P6", "P20": "PlaB P20"},
        figsize=(15, max(5, int(len(result_df) * 0.3))),
        ylabel="Fold Change (2^log2FC)",
        do_test=False,
        use_broken_axis=False
    )
if __name__ == "__main__":
    main()