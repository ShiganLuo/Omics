import argparse
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.LogUtil import setup_logger
from utils.VEP_SV import read_vep_tab
from utils.SV_TYPE_plot import plot_group_type_comparison
import pandas as pd
from typing import List, Dict, Optional, Tuple, Callable, Literal
PlotFormat = Literal["png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg", "pgf", "raw", "rgba"]
import logging
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact, ttest_ind, mannwhitneyu

logger = setup_logger("OncoPrint", level=logging.INFO)

def tab_parser(
        table_file: str,
        cosmic_file:str,
        tab_gene_col: str = "SYMBOL",
        cosmic_gene_col:str = "GENE_SYMBOL",
        sv_type_pattern:str = r"(DEL|DUP|INV|INS|TRA)",
        **kwargs
):
    """
    Parse a VEP annotation table and filter for COSMIC cancer genes.

    Parameters
    ----------
    table_file : str
        Path to the VEP annotation table (tab-delimited format).
    cosmic_file : str
        Path to the COSMIC cancer gene list (tab-delimited format).
    tab_gene_col : str, default="SYMBOL"
        Column name in the VEP table that contains gene symbols.
    cosmic_gene_col : str, default="GENE_SYMBOL"
        Column name in the COSMIC file that contains gene symbols.
    sv_type_pattern : str, default=r"(DEL|DUP|INV|INS|TRA)"
        Regex pattern to extract SV type from the variation column.
    **kwargs
        Additional keyword arguments passed to ``read_vep_tab``.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``[tab_gene_col, 'svtype']`` containing only
        rows where the gene symbol is present in the COSMIC cancer gene list.
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
    Prepare data for OncoPrint visualization from multiple VEP tables.

    Parses each sample's VEP annotation table, filters for COSMIC cancer
    genes, and builds a gene-by-sample matrix of SV types suitable for
    OncoPrint plotting.

    Parameters
    ----------
    tab_files : dict of str to str
        Mapping of sample names to VEP annotation file paths.
    cosmic_file : str
        Path to the COSMIC cancer gene list file.
    gene_col : str, default="SYMBOL"
        Column name for gene symbols in the VEP tables.
    sample_order : list of str, optional
        Explicit ordering of sample columns in the output matrix.
        If ``None``, columns follow insertion order of ``tab_files``.
    **kwargs
        Additional keyword arguments forwarded to ``tab_parser``.

    Returns
    -------
    pd.DataFrame
        Gene-by-sample matrix where values are comma-separated SV type
        strings (e.g. ``"DEL,DUP"``) and ``NaN`` indicates no alteration.
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
    condition_files: Dict[str, str],
    oncoprint_file: str,
    oncoprint_gene_col: str = "SYMBOL",
    deseq2_gene_col: str = "name",
    padj_col: str = "padj",
    padj_threshold: Optional[float] = 0.05,
    log2fc_col: str = "log2FoldChange",
    sep: str = "\t"
) -> pd.DataFrame:
    """
    Summarize log2 fold change values of COSMIC genes across multiple conditions.

    This function reads multiple DESeq2 result files and extracts log2 fold change
    values for genes present in an oncoprint gene list. It supports an arbitrary
    number of conditions and merges results into a single wide-format DataFrame.
    Missing values are filled with NA.

    Parameters
    ----------
    condition_files : Dict[str, str]
        Mapping of condition names to DESeq2 result file paths,
        e.g. ``{"P6": "p6.tsv", "P20": "p20.tsv"}``.
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

    Returns
    -------
    pd.DataFrame
        A DataFrame with genes as rows and conditions as columns,
        containing fold change values. Missing values are NA.
    """

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

    for condition, file in condition_files.items():
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
    Convert a p-value to a significance star label.

    Parameters
    ----------
    p : float
        The p-value to convert.

    Returns
    -------
    str
        Significance label: ``"****"`` for p < 1e-4, ``"***"`` for p < 1e-3,
        ``"**"`` for p < 1e-2, ``"*"`` for p < 0.05, and ``"ns"`` otherwise.
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
    group_order: tuple,
    method: str = "chi2"
) -> Dict[str, str]:
    """
    Compute significance stars per SV type via pairwise comparisons.

    Parameters
    ----------
    pivot : pandas.DataFrame
        Pivot table with SV types as index and groups as columns.
    group_order : tuple of str
        Group names to compare.
    method : str, default="chi2"
        Statistical test method.

    Returns
    -------
    dict
        Mapping: svtype -> best significance star across all pairs.
    """
    stars = {}
    group_list = list(group_order)

    for sv in pivot.index:
        best_p = 1.0
        for i, g1 in enumerate(group_list):
            for g2 in group_list[i + 1:]:
                if method in ("chi2", "fisher"):
                    table = np.array([
                        [pivot.loc[sv, g1], pivot[g1].sum() - pivot.loc[sv, g1]],
                        [pivot.loc[sv, g2], pivot[g2].sum() - pivot.loc[sv, g2]],
                    ])
                else:
                    table = np.array([
                        [pivot.loc[sv, g1]],
                        [pivot.loc[sv, g2]],
                    ])

                p = compute_pvalue(table, method)
                best_p = min(best_p, p)

        stars[sv] = p_to_star(best_p)

    return stars



# =========================
# OncoPrint visualization
# =========================
def plot_oncoprint(
    matrix: pd.DataFrame,
    out_prefix: str,
    sv_colors: Optional[Dict[str, str]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: int = 300,
    image_formats: Optional[List[PlotFormat]] = None,
    title: str = "",
) -> None:
    """Plot an OncoPrint heatmap from a gene-by-sample SV type matrix.

    Draws a grid where rows are genes (top → bottom), columns are samples,
    and cells are colored by SV type.  Margin bar charts show per-gene
    mutation frequency (right) and per-sample mutation counts (top).

    Parameters
    ----------
    matrix : pd.DataFrame
        Gene-by-sample matrix. Index = gene names, columns = sample names.
        Values are comma-separated SV type strings (e.g. ``"DEL,DUP"``)
        or ``NaN`` for no alteration.
    out_prefix : str
        Output file path prefix. Each format is saved as
        ``{out_prefix}.{fmt}``.
    sv_colors : dict, optional
        Mapping from SV type to color. Defaults to a curated palette
        for common SV types (DEL, DUP, INS, INV, TRA, OTHER).
    figsize : tuple of float, optional
        Figure size ``(width, height)``. Auto-calculated if ``None``.
    dpi : int, optional
        Resolution in dots per inch. Default is ``300``.
    image_formats : list of str, optional
        Output image formats. Defaults to ``["png"]``.
    title : str, optional
        Plot title. Default is ``"OncoPrint"``.

    Returns
    -------
    None
        Saves the figure to ``{out_prefix}.{fmt}`` for each format.
    """
    if image_formats is None:
        image_formats = ["png"]

    if sv_colors is None:
        sv_colors = {
            "DEL": "#E74C3C",
            "DUP": "#3498DB",
            "INS": "#2ECC71",
            "INV": "#9B59B6",
            "TRA": "#F39C12",
            "OTHER": "#95A5A6",
        }

    # Collect all SV types present
    all_types = set()
    for val in matrix.values.flat:
        if pd.notna(val):
            all_types.update(str(val).split(","))

    genes = matrix.index.tolist()
    samples = matrix.columns.tolist()
    n_genes = len(genes)
    n_samples = len(samples)

    if figsize is None:
        figsize = (max(6, n_samples * 0.6 + 2), max(4, n_genes * 0.35 + 1.5))

    fig = plt.figure(figsize=figsize)

    # Use fixed ratios: main grid gets most space, margins are small
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[4, 1],
        height_ratios=[1, 4],
        wspace=0.04, hspace=0.04,
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[1, 0])
    ax_right = fig.add_subplot(gs[1, 1])

    # ========== Main grid ==========
    # Use imshow-friendly coordinate system:
    #   row 0 = top gene, row n-1 = bottom gene
    #   col 0 = left sample, col n-1 = right sample
    # We place rectangles at integer (col, row) centers.
    default_color = "#ECF0F1"
    ax_main.set_facecolor("#FAFAFA")

    for i in range(n_genes):          # row index (top=0)
        for j in range(n_samples):    # col index
            ax_main.add_patch(plt.Rectangle(
                (j - 0.5, i - 0.5), 1, 1,
                facecolor=default_color, edgecolor="white", linewidth=0.5,
            ))

    for i, gene in enumerate(genes):
        for j, sample in enumerate(samples):
            val = matrix.loc[gene, sample]
            if pd.isna(val):
                continue
            types = str(val).split(",")
            n_types = len(types)
            for k, sv_type in enumerate(types):
                sv_type = sv_type.strip()
                color = sv_colors.get(sv_type, sv_colors.get("OTHER", "#95A5A6"))
                # Stack: k=0 at bottom, k=n-1 at top of the cell
                y_lo = i - 0.5 + k / n_types
                height = 1.0 / n_types
                ax_main.add_patch(plt.Rectangle(
                    (j - 0.5, y_lo), 1, height,
                    facecolor=color, edgecolor="none",
                ))

    ax_main.set_xlim(-0.5, n_samples - 0.5)
    ax_main.set_ylim(n_genes - 0.5, -0.5)    # row 0 at top
    ax_main.set_xticks(range(n_samples))
    ax_main.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax_main.set_yticks(range(n_genes))
    ax_main.set_yticklabels(genes, fontsize=8)
    ax_main.set_xlabel("")
    ax_main.set_ylabel("")
    ax_main.tick_params(axis="both", length=0)

    # ========== Right margin: mutation frequency per gene ==========
    mut_counts = matrix.notna().sum(axis=1)
    mut_freq = mut_counts / n_samples * 100

    ax_right.set_ylim(n_genes - 0.5, -0.5)    # same orientation as main
    ax_right.barh(range(n_genes), mut_freq, color="#2C76FF", height=0.6)
    ax_right.set_xlim(0, 100)
    ax_right.set_xlabel("%", fontsize=8)
    ax_right.set_yticks([])
    ax_right.tick_params(axis="y", length=0)
    ax_right.tick_params(axis="x", labelsize=7)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.invert_xaxis()

    # ========== Top margin: mutation count per sample ==========
    sample_counts = matrix.notna().sum(axis=0)

    ax_top.set_xlim(-0.5, n_samples - 0.5)    # same orientation as main
    ax_top.bar(range(n_samples), sample_counts, color="#2C76FF", width=0.6)
    ax_top.set_xticks([])
    ax_top.tick_params(axis="x", length=0)
    ax_top.tick_params(axis="y", labelsize=7)
    ax_top.set_ylabel("Mutations", fontsize=8)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    # ========== Legend ==========
    legend_types = [t for t in sv_colors if t in all_types]
    if not legend_types:
        legend_types = list(all_types)
    patches = [mpatches.Patch(color=sv_colors.get(t, "#95A5A6"), label=t)
               for t in legend_types]
    ax_top.legend(
        handles=patches, loc="upper right",
        fontsize=7, frameon=False, ncol=len(legend_types),
        bbox_to_anchor=(1.0, 1.3),
    )

    ax_main.set_title(title, fontsize=11, pad=30)

    plt.tight_layout()
    outdir = os.path.dirname(out_prefix)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    for fmt in image_formats:
        plt.savefig(f"{out_prefix}.{fmt}", dpi=dpi, bbox_inches="tight")
    plt.close()
    logger.info(f"OncoPrint saved: {out_prefix}.[{','.join(image_formats)}]")


def parser_args(
        
) -> argparse.Namespace:
    """
    Parse command-line arguments for OncoPrint analysis.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes:

        - ``tab_files`` : list of str -- VEP annotation tab files (``sample=path``).
        - ``cosmic_file`` : str -- path to COSMIC cancer gene list file.
        - ``oncoprint_file`` : str -- path to output oncoprint matrix file.
        - ``deseq2_files`` : list of str -- DESeq2 result files (``condition=path``).
        - ``formats`` : list of str -- image output formats (default ``["png"]``).
    """
    parser = argparse.ArgumentParser(description="Run OncoPrint analysis for SV data")
    parser.add_argument(
        "-g", "--group",
        action="append",
        required=True,
        metavar="NAME:VCF",
        help="Group in format 'name:vep_annotation_file(tab)'. Can be specified multiple times. "
             "Example: -g Control:ctrl.vep.txt -g Experiment:exp.vep.txt",
    )
    parser.add_argument("--cosmic_file", required=True, help="Path to COSMIC cancer gene list file")
    parser.add_argument("-o","--outdir", required=True, help="Path to output directory")
    parser.add_argument(
        "-d", "--deseq2_files",
        action="append",
        required=True,
        metavar="CONDITION:PATH",
        help="DESeq2 result files (format: condition:path)",
    )
    parser.add_argument("-f", "--image_format", action="append", dest="image_format", metavar="FMT", help="Image output format (png, pdf, svg, ...). Can be specified multiple times. Default: png.")
    return parser.parse_args()

def main():
    """Run the OncoPrint analysis pipeline from command-line arguments.

    Parses CLI arguments, builds the OncoPrint matrix from VEP annotation
    tables and COSMIC gene list, summarizes DESeq2 fold changes, and
    generates a comparison bar plot.
    """
    args = parser_args()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    # Parse group_files from "sample:path" format
    tab_files = {}
    for item in args.group:
        if ":" not in item:
            raise ValueError(f"Invalid group files format '{item}', expected 'sample:path'")
        sample, path = item.split(":", 1)
        tab_files[sample] = path
    sample_order = list(tab_files.keys())

    # Build oncoprint matrix
    oncoprint_file = os.path.join(outdir,"oncoprint_matrix.csv")

    oncoprint_matrix = prepare_oncoprint_data(
        tab_files,
        cosmic_file=args.cosmic_file,
        sample_order=sample_order,
    )

    oncoprint_matrix.to_csv(oncoprint_file)
    logger.info(f"Oncoprint matrix saved to: {oncoprint_file}")

    # Plot OncoPrint
    plot_oncoprint(
        matrix=oncoprint_matrix,
        out_prefix=os.path.join(outdir, "oncoprint"),
        image_formats=args.image_format,
    )

    # Parse deseq2_files from "condition=path" format
    condition_files = {}
    for item in args.deseq2_files:
        if ":" not in item:
            raise ValueError(f"Invalid deseq2_files format '{item}', expected 'condition:path'")
        cond, path = item.split(":", 1)
        condition_files[cond] = path

    # DESeq2 fold change summary
    result_df = Deseq2_oncoprint_data(
        condition_files=condition_files,
        oncoprint_file=oncoprint_file,
    )

    out_csv = os.path.join(outdir, "deseq2_oncoprint_fold_change_all.csv")
    result_df.to_csv(out_csv, index=False)
    logger.info(f"DESeq2 fold change table saved to: {out_csv}")

    # Plot
    result_df = result_df.dropna()
    for fmt in args.image_format:
        out_plot = os.path.join(outdir, f"deseq2_oncoprint_fold_change_all.{fmt}")
        logger.info(f"Plotting fold change comparison to: {out_plot}")
        plot_group_type_comparison(
            df=result_df.melt(id_vars="SYMBOL", var_name="group", value_name="fold_change"),
            out_png=out_plot,
            group_col="group",
            type_col="SYMBOL",
            count_col="fold_change",
            group_order=list(condition_files.keys()),
            type_order=result_df["SYMBOL"].tolist(),
            figsize=(15, max(5, int(len(result_df) * 0.3))),
            ylabel="Fold Change (2^log2FC)",
            do_test=False,
            use_broken_axis=False,
        )
if __name__ == "__main__":
    main()