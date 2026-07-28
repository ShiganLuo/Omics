import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import chi2_contingency, fisher_exact, ttest_ind, mannwhitneyu
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
