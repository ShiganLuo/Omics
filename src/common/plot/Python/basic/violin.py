from itertools import combinations
from statannotations.Annotator import Annotator
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd
from scipy.stats import mannwhitneyu
import itertools
import numpy as np
from typing import Optional, List, Tuple, Dict, Any
import logging
logger = logging.getLogger(__name__)

def single_violin_plot(
    df: pd.DataFrame,
    value: str,
    key: Optional[str] = None,
    outfile: str = "violin.png",
    xlabel: str = "",
    ylabel: str = "Value",
    order: Optional[List[str]] = None,
    sort_keys: bool = True,
    palette: str = "Set2",
    figsize: Tuple[int, int] = (6, 6),
    violin_alpha: float = 0.85,
    point_size: float = 2,
    show_median: bool = False,
    threshold: Optional[float] = None,
):
    """Draw a single (or per-category) violin + strip plot.

    When ``key`` is ``None`` a single violin is drawn for the entire
    ``value`` column.  When ``key`` is provided, one violin is drawn
    per unique category in that column.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data.
    value : str
        Column name for the numeric values (y-axis).
    key : str or None, optional
        Column name for x-axis categories.  ``None`` draws one
        overall violin.
    outfile : str, optional
        Output image path (PNG/PDF).
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    order : list of str or None, optional
        Explicit category ordering on the x-axis.  Ignored when
        ``key`` is ``None``.
    sort_keys : bool, optional
        Sort categories alphabetically when ``order`` is ``None``.
    palette : str, optional
        Seaborn / matplotlib colormap name.
    figsize : Tuple[int, int], optional
        Figure size in inches.
    violin_alpha : float, optional
        Fill transparency of violin bodies.
    point_size : float, optional
        Strip-plot dot size.
    show_median : bool, optional
        Draw a horizontal median line inside each violin.
    threshold : float or None, optional
        When set, draws a red dashed horizontal line at this value
        and displays in the legend what percentage of data points
        fall below it (overall when ``key`` is ``None``, per
        category otherwise).

    Raises
    ------
    ValueError
        If ``value`` (or ``key``, when provided) is not in
        ``df.columns``.
    """
    if value not in df.columns:
        raise ValueError(f"Value column not found: {value}")
    if key is not None and key not in df.columns:
        raise ValueError(f"Category column not found: {key}")
    logger.info(f"Drawing violin plot for column: {value}")
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    sns.set_style("whitegrid")

    if key is not None:
        if key not in df.columns:
            raise ValueError(f"Category column not found: {key}")
        plot_df = df[[key, value]].dropna().copy()
        if order is None:
            observed = plot_df[key].unique().tolist()
            order = sorted(observed) if sort_keys else observed
    else:
        plot_df = df[[value]].dropna().copy()
        plot_df["_cat"] = ""
        key = "_cat"
        order = [""]

    n_cats = len(order)
    base_colors = sns.color_palette(palette, n_colors=max(n_cats, 1))

    fig, ax = plt.subplots(figsize=figsize)

    # ── violin ──────────────────────────────────────────────
    vp = sns.violinplot(
        data=plot_df,
        x=key,
        y=value,
        cut=0,
        order=order,
        palette=base_colors,
        inner="box",
        linewidth=1.2,
        width=0.7,
        bw_adjust=1.2,
        density_norm="area",
        ax=ax,
    )
    for poly in ax.collections:
        try:
            poly.set_alpha(violin_alpha)
            poly.set_edgecolor("black")
            poly.set_linewidth(1.0)
        except Exception:
            pass

    # ── strip ───────────────────────────────────────────────
    sns.stripplot(
        data=plot_df,
        x=key,
        y=value,
        order=order,
        color="black",
        size=point_size,
        alpha=0.9,
        jitter=0.12,
        linewidth=0,
        ax=ax,
    )

    # ── median line ─────────────────────────────────────────
    if show_median:
        for i, cat in enumerate(order):
            vals = plot_df.loc[plot_df[key] == cat, value]
            if len(vals) > 0:
                med = float(np.median(vals))
                ax.hlines(
                    y=med,
                    xmin=i - 0.18,
                    xmax=i + 0.18,
                    color="black",
                    linewidth=2.5,
                    zorder=5,
                )

    # ── threshold line ──────────────────────────────────────
    if threshold is not None:
        ax.axhline(
            y=threshold,
            color="red",
            linestyle="--",
            linewidth=1.5,
            alpha=0.8,
            zorder=4,
        )
        # Compute proportion below threshold
        below = plot_df[value] < threshold
        pct_below = 100.0 * below.sum() / len(plot_df)
        if key == "_cat":
            # Single violin: one overall percentage
            label = f"< {threshold}: {pct_below:.1f}%"
        else:
            # Per-category percentages
            parts = []
            for cat in order:
                cat_vals = plot_df.loc[plot_df[key] == cat, value]
                cat_pct = 100.0 * (cat_vals < threshold).sum() / len(cat_vals) if len(cat_vals) > 0 else 0.0
                parts.append(f"{cat} {cat_pct:.1f}%")
            label = f"< {threshold}: " + ", ".join(parts)
        ax.legend(
            [plt.Line2D([0], [0], color="red", linestyle="--", linewidth=1.5)],
            [label],
            loc="upper right",
            fontsize=10,
            frameon=True,
            facecolor="white",
            edgecolor="none",
            framealpha=0.9,
        )

    # ── axes ────────────────────────────────────────────────
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.25)

    if key == "_cat":
        ax.set_xticks([])

    fig.tight_layout()
    fig.savefig(outfile, dpi=600, bbox_inches="tight")
    plt.close(fig)

def paired_violin_plot(
    df: pd.DataFrame,
    key: str,
    values: List[str],
    outfile: str,
    xlabel: str = "Group",
    ylabel: str = "Value",
    group_names: Optional[List[str]] = None,
    palette: Tuple[str, str] = ("#2166AC", "#B2182B"),
    figsize: Tuple[int, int] = (7, 8),
    violin_alpha: float = 0.25,
    violin_width: float = 0.7,
    point_size: float = 28,
    line_alpha: float = 0.35,
    line_width: float = 0.8,
    jitter: float = 0.04,
):
    if key not in df.columns:
        raise ValueError(
            f"sample id column not found: {key}"
        )

    if len(values) != 2:
        raise ValueError(
            "`values` must contain exactly two columns"
        )

    missing_cols = [
        v for v in values
        if v not in df.columns
    ]

    if missing_cols:
        raise ValueError(
            f"columns not found: {missing_cols}"
        )
    if group_names is None:
        group_names = values

    if len(group_names) != 2:
        raise ValueError(
            "`group_names` must contain two names"
        )
    df_plot = (
        df[[key] + values]
        .dropna(subset=values)
        .copy()
    )
    long_parts = []

    for value_col, group_name in zip(
        values,
        group_names,
    ):

        sub_df = pd.DataFrame(
            {
                "sample_id": df_plot[key],
                "group": group_name,
                "value": df_plot[value_col],
            }
        )

        long_parts.append(sub_df)

    long_df = pd.concat(
        long_parts,
        ignore_index=True,
    )
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    sns.set_style("whitegrid")

    fig, ax = plt.subplots(
        figsize=figsize
    )
    violin = sns.violinplot(
        data=long_df,
        x="group",
        y="value",
        order=group_names,
        palette=list(palette),
        cut=1,
        inner=None,
        linewidth=1.5,
        width=violin_width,
        bw_adjust=1.2,
        density_norm="area",
        ax=ax,
    )
    for poly in violin.collections:

        try:
            poly.set_alpha(violin_alpha)
            poly.set_edgecolor("black")
            poly.set_linewidth(1.2)
        except Exception:
            pass
    x_positions = {
            group_names[0]: 0,
            group_names[1]: 1,
        }
    rng = np.random.default_rng(123)

    for _, row in df_plot.iterrows():


        x1 = (
            x_positions[group_names[0]]
            + rng.uniform(-jitter, jitter)
        )

        x2 = (
            x_positions[group_names[1]]
            + rng.uniform(-jitter, jitter)
        )

        y1 = row[values[0]]
        y2 = row[values[1]]

        # paired line
        ax.plot(
            [x1, x2],
            [y1, y2],
            color="gray",
            linewidth=line_width,
            alpha=line_alpha,
            zorder=1,
        )

        # left point
        ax.scatter(
            x1,
            y1,
            s=point_size,
            color="black",
            alpha=0.9,
            zorder=3,
        )

        # right point
        ax.scatter(
            x2,
            y2,
            s=point_size,
            color="black",
            alpha=0.9,
            zorder=3,
        )
    medians = [
        np.median(df_plot[v])
        for v in values
    ]

    for i, median in enumerate(medians):

        ax.hlines(
            y=median,
            xmin=i - 0.18,
            xmax=i + 0.18,
            color="black",
            linewidth=2.5,
            zorder=5,
        )
    ax.set_xlabel(
        xlabel,
        fontsize=14,
    )

    ax.set_ylabel(
        ylabel,
        fontsize=14,
    )
    ax.tick_params(
        axis="x",
        labelsize=12,
    )

    ax.tick_params(
        axis="y",
        labelsize=11,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.25,
    )
    fig.tight_layout()
    fig.savefig(
        outfile,
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)

def adjust_group_color(
    color,
    group_idx,
    n_groups,
    lightness_strength: float = 0.25,
    saturation_strength: float = 0.35,
    max_lightness: float = 0.75,
):
    """
    Adjust group color while keeping
    mutation hue stable.

    Parameters
    ----------
    color :
        Base RGB tuple.

    group_idx :
        Group index.

    n_groups :
        Total group count.

    lightness_strength :
        How much brighter later groups become.

        Larger value:
            stronger lightness contrast.

    saturation_strength :
        How much saturation decreases.

        Larger value:
            stronger pale effect.

    max_lightness :
        Upper limit of brightness.

        Prevents colors becoming nearly white.
    """

    r, g, b = color

    h, l, s = colorsys.rgb_to_hls(
        r,
        g,
        b,
    )

    if n_groups == 1:
        return color

    ratio = group_idx / (n_groups - 1)

    # increase lightness
    l = min(
        max_lightness,
        l + ratio * lightness_strength,
    )

    # reduce saturation
    s = max(
        0.05,
        s - ratio * saturation_strength,
    )

    return colorsys.hls_to_rgb(
        h,
        l,
        s,
    )


def violin_plot_advanced(
    df: pd.DataFrame,
    key: str,
    values: List[str],
    outfile: str,
    xlabel: str = "SV Type",
    ylabel: str = "Delta Frequency",
    order: Optional[List[str]] = None,
    sort_keys: bool = True,
    group_names: Optional[List[str]] = None,
    group_order: Optional[List[str]] = None,
    base_palette: str = "Set2",
    figsize: Tuple[int, int] = (16, 8),
):
    """
    Draw violin + strip plot.

    Features
    --------
    - Different mutations:
        different base colors

    - Different groups:
        different saturation/lightness
        within same mutation

    - Black solid strip points

    - Publication-style aesthetics
    """

    # ==========================================================
    # matplotlib settings
    # ==========================================================

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    sns.set_style("whitegrid")

    # ==========================================================
    # validation
    # ==========================================================

    if key not in df.columns:
        raise ValueError(
            f"x-axis column not found: {key}"
        )

    if not values:
        raise ValueError(
            "`values` cannot be empty"
        )

    missing_cols = [
        v for v in values
        if v not in df.columns
    ]

    if missing_cols:
        raise ValueError(
            f"y-axis columns not found: {missing_cols}"
        )

    if (
        group_names is not None
        and len(group_names) != len(values)
    ):
        raise ValueError(
            "`group_names` length must equal "
            "`values` length"
        )

    # ==========================================================
    # build long dataframe
    # ==========================================================

    long_parts = []

    for i, v_col in enumerate(values):

        group_name = (
            group_names[i]
            if group_names is not None
            else v_col
        )

        sub_df = (
            df[[key, v_col]]
            .dropna(subset=[key, v_col])
            .copy()
        )

        sub_df = sub_df.rename(
            columns={v_col: "_y"}
        )

        sub_df["group"] = group_name

        long_parts.append(sub_df)

    if not long_parts:
        raise ValueError(
            "No valid plotting data"
        )

    df_plot = pd.concat(
        long_parts,
        ignore_index=True,
    )

    # ==========================================================
    # x-axis order
    # ==========================================================

    if order is None:

        observed_order = (
            df_plot[key]
            .dropna()
            .unique()
            .tolist()
        )

        if sort_keys:
            order = sorted(observed_order)
        else:
            order = observed_order

    # ==========================================================
    # group order
    # ==========================================================

    if group_order is None:

        group_order = (
            df_plot["group"]
            .dropna()
            .unique()
            .tolist()
        )

    # ==========================================================
    # base colors
    # ==========================================================

    base_colors = sns.color_palette(
        base_palette,
        n_colors=len(order),
    )

    mutation_color_map = {
        mutation: color
        for mutation, color in zip(
            order,
            base_colors,
        )
    }

    # ==========================================================
    # color adjustment
    # ==========================================================

    

    # ==========================================================
    # plotting
    # ==========================================================

    fig, ax = plt.subplots(
        figsize=figsize
    )

    # ----------------------------------------------------------
    # violin plot
    # ----------------------------------------------------------

    sns.violinplot(
        data=df_plot,
        x=key,
        y="_y",
        hue="group",
        order=order,
        hue_order=group_order,
        dodge=True,
        inner="quartile",
        linewidth=1.2,
        saturation=1,
        bw_adjust=1.2,
        density_norm="area",
        width=0.95,
        ax=ax,
    )

    # ==========================================================
    # recolor violins
    # ==========================================================

    violin_bodies = [
        c for c in ax.collections
        if isinstance(
            c,
            mcollections.PolyCollection
        )
    ]

    violin_idx = 0

    for mutation in order:

        base_color = mutation_color_map[
            mutation
        ]

        for group_idx, group in enumerate(
            group_order
        ):

            if violin_idx >= len(
                violin_bodies
            ):
                break

            color = adjust_group_color(
                base_color,
                group_idx,
                len(group_order),
            )

            violin = violin_bodies[
                violin_idx
            ]

            violin.set_facecolor(color)

            violin.set_edgecolor(
                "black"
            )

            violin.set_alpha(0.95)

            violin.set_linewidth(1)

            violin_idx += 1

    # ==========================================================
    # black solid strip points
    # ==========================================================

    sns.stripplot(
        data=df_plot,
        x=key,
        y="_y",
        hue="group",
        order=order,
        hue_order=group_order,
        dodge=True,
        jitter=0.12,
        color="black",
        size=2.5,
        alpha=0.9,
        linewidth=0,
        ax=ax,
    )

    # ==========================================================
    # remove duplicated legends
    # ==========================================================

    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # ==========================================================
    # explanatory legend
    # ==========================================================

    example_mutation = order[0]

    example_base_color = mutation_color_map[
        example_mutation
    ]

    legend_handles = []

    n_groups = len(group_order)

    if n_groups == 2:

        shade_words = [
            "Dark",
            "Light",
        ]

    else:

        shade_words = [
            f"Shade {i + 1}"
            for i in range(n_groups)
        ]

    for group_idx, group in enumerate(
        group_order
    ):

        legend_color = adjust_group_color(
            example_base_color,
            group_idx,
            n_groups,
        )

        label = (
            f"{shade_words[group_idx]} "
            f"{group}"
        )

        legend_handles.append(
            mpatches.Patch(
                facecolor=legend_color,
                edgecolor="black",
                label=label,
            )
        )

    ax.legend(
        handles=legend_handles,
        title="Color meaning",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=11,
        title_fontsize=12,
    )

    # ==========================================================
    # labels
    # ==========================================================

    ax.set_xlabel(
        xlabel,
        fontsize=14,
    )

    ax.set_ylabel(
        ylabel,
        fontsize=14,
    )

    # ==========================================================
    # ticks
    # ==========================================================

    ax.tick_params(
        axis="x",
        rotation=45,
        labelsize=11,
    )

    ax.tick_params(
        axis="y",
        labelsize=11,
    )

    # ==========================================================
    # style
    # ==========================================================

    ax.spines["top"].set_visible(False)

    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
    )

    # ==========================================================
    # layout
    # ==========================================================

    fig.tight_layout()

    # ==========================================================
    # save
    # ==========================================================

    fig.savefig(
        outfile,
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)

def violin_plot(
    df: pd.DataFrame,
    key: str,
    value: str,
    outfile:str,
    xlabel: str = "SV Type",
    ylabel: str = "Delta Frequency (Freq - ddPCR_AF)",
):
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.violinplot(
        x=key,
        y=value,
        data=df,
        inner='box',
        linewidth=1.5,
        hue=key,
        palette='Set2',
        ax=ax
    )
    sns.stripplot(
        x=key,
        y=value,
        data=df,
        color='black',
        size=3,
        alpha=0.5,
        ax=ax
    )
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(axis='x', labelrotation=45)
    fig.savefig(outfile, dpi=300, bbox_inches='tight')

def violin_test_df(
        df:pd.DataFrame,
        outplot:str,
        key:str,
        value:str
):
    unique_cell_types = df[key].unique().tolist()
    box_pairs = list(combinations(unique_cell_types, 2))
    # combinations(unique_cell_types, 2) 会生成所有不重复的 (组A, 组B) 组合

    plt.figure(figsize=(10, 8))
    ax = sns.violinplot(
        x=key,
        y=value,
        data=df,
        palette='Set2',
        inner='box',
        linewidth=1.5
    )
    annotator = Annotator(
        ax,
        box_pairs,
        data=df,
        x=key,
        y=value,
        order=unique_cell_types # 确保 order 与图中的顺序一致
    )
    # 注意：进行多重比较时，通常需要进行多重检验校正 (Multiple Testing Correction)。
    annotator.configure(
        test='Mann-Whitney',
        text_format='star',
        loc='inside',
        # 多重检验校正：例如使用 Benjamini/Hochberg FDR
        comparisons_correction='fdr_bh', # <--- 关键参数！
        pvalue_thresholds=[
            [0.05, "*"],
            [0.01, "**"],
            [0.001, "***"],
            [0.0001, "****"]
        ]
    )

    annotator.apply_and_annotate()
    ax.set_xlabel('') 
    plt.ylabel('GSI')
    plt.tight_layout()
    
    plt.savefig(outplot,dpi=300)


def plot_comparison_with_violin(
    result: Dict[str, Any],
    scatter_size: int = 10,
    jitter_width:float = 0.3,
    title: str = "Mapping Rate Comparison",
    xlabel: str = "Sample Type",
    ylabel: str = "Mapping Rate (%)",
    figsize: Tuple[int, int] = (10, 7),
    palette: Tuple[str, str] = ("#3498db", "#e74c3c"),
    show_stats: bool = True,
    save_path: Optional[str] = None,
    dpi: int = 300,
    verbose: bool = True
) -> Figure:
    """
    Create a violin plot with embedded box plot for comparing two groups.

    Parameters
    ----------
    result : Dict[str, Any]
        Result dictionary from prepare_data().
    title : str, optional
        Plot title.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    figsize : Tuple[int, int], optional
        Figure size.
    palette : Tuple[str, str], optional
        Colors for the two groups.
    show_stats : bool, optional
        Whether to show statistics annotation.
    save_path : str, optional
        Path to save the figure.
    dpi : int, optional
        DPI for saved figure.
    verbose : bool, optional
        Whether to logger.info information.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure object.
    """
    # Configure Chinese font
    selected_font = configure_chinese_font()
    if verbose and selected_font:
        logger.info(f"Using Chinese font: {selected_font}")

    # Extract data
    group1_name = result['group1_name']
    group2_name = result['group2_name']
    group1_data = result['group1_data']
    group2_data = result['group2_data']
    comparison = result['comparison']

    # Set style
    sns.set_style("whitegrid")
    apply_font_after_style(selected_font)
    plt.rcParams['font.size'] = 12

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data for seaborn
    df_plot = pd.DataFrame({
        'value': group1_data + group2_data,
        'group': [group1_name] * len(group1_data) + [group2_name] * len(group2_data)
    })

    # Create violin plot using seaborn (better KDE control)
    sns.violinplot(
        data=df_plot,
        x='group',
        y='value',
        palette=palette,
        bw_adjust=0.8,   # Adjust bandwidth for smoother KDE
        inner='box',     # Show box plot inside violin
        alpha=0.7,
        linewidth=1.5,
        ax=ax
    )

    # Add individual points with jitter
    for i, (data, pos, color) in enumerate(zip([group1_data, group2_data], [0, 1], palette)):
        jitter = np.random.uniform(-jitter_width, jitter_width, size=len(data))
        ax.scatter(
            np.full_like(data, pos) + jitter,
            data,
            color=color,
            alpha=0.6,
            s=scatter_size,
            edgecolors='black',
            linewidths=0.5,
            zorder=5
        )

    # Set labels with explicit font properties
    from matplotlib.font_manager import FontProperties
    if selected_font:
        font_prop = FontProperties(family=selected_font, size=12)
        font_prop_bold = FontProperties(family=selected_font, size=14, weight='bold')
        font_prop_title = FontProperties(family=selected_font, size=16, weight='bold')
    else:
        font_prop = FontProperties(size=12)
        font_prop_bold = FontProperties(size=14, weight='bold')
        font_prop_title = FontProperties(size=16, weight='bold')

    # Seaborn uses 0-based positions
    ax.set_xticklabels([f"{group1_name}\n(n={len(group1_data)})", f"{group2_name}\n(n={len(group2_data)})"],
                       fontproperties=font_prop)
    ax.set_xlabel(xlabel, fontproperties=font_prop_bold)
    ax.set_ylabel(ylabel, fontproperties=font_prop_bold)
    ax.set_title(title, fontproperties=font_prop_title, pad=20)

    # Add statistics if requested
    if show_stats:
        p_value = comparison['p_value']
        significant = comparison['significant']
        effect_size = comparison['effect_size']
        effect_interp = comparison['effect_size_interpretation']

        if p_value < 0.001:
            sig_stars = "***"
        elif p_value < 0.01:
            sig_stars = "**"
        elif p_value < 0.05:
            sig_stars = "*"
        else:
            sig_stars = "ns"

        y_max = max(max(group1_data), max(group2_data))
        y_min = min(min(group1_data), min(group2_data))
        y_range = y_max - y_min

        bracket_y = y_max + y_range * 0.15
        star_y = bracket_y + y_range * 0.05

        # Use 0-based positions for seaborn
        ax.plot([0, 0, 1, 1], [bracket_y, bracket_y + y_range * 0.02, bracket_y + y_range * 0.02, bracket_y],
                color='black', linewidth=1.5)

        ax.text(0.5, star_y, sig_stars, ha='center', va='bottom', fontsize=16, fontweight='bold')

        p_text = f"p = {p_value:.4f}" if p_value >= 0.001 else f"p < 0.001"
        ax.text(0.5, star_y - y_range * 0.08, p_text, ha='center', va='top', fontsize=11, style='italic')

        # Add effect size info - use axes transform for better positioning
        effect_text = f"{comparison['effect_size_name']}: {effect_size:.2f} ({effect_interp})"
        ax.text(0.98, 0.98, effect_text, transform=ax.transAxes, fontsize=10, 
                color='gray', ha='right', va='top', style='italic')

        # Add test name - use axes transform for better positioning
        test_text = f"Test: {comparison['test_name']}"
        ax.text(0.02, 0.02, test_text, transform=ax.transAxes, fontsize=9, color='gray', style='italic')

        # Remove ylim restriction to avoid truncating violin plot
        # ax.set_ylim(y_min - y_range * 0.1, y_max + y_range * 0.35)

    # Add mean annotations with explicit font properties
    from matplotlib.font_manager import FontProperties
    if selected_font:
        font_prop = FontProperties(family=selected_font)
    else:
        font_prop = FontProperties()

    for i, (data, pos, name) in enumerate(zip([group1_data, group2_data], [0, 1], [group1_name, group2_name])):
        mean_val = np.mean(data)
        ax.text(pos, mean_val, f"  μ={mean_val:.1f}%", va='center', fontsize=10, fontweight='bold', fontproperties=font_prop)

    # Clean up
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        if verbose:
            logger.info(f"Figure saved to: {save_path}")

    return fig



def violin_by_method_celltype_step(df_dict: Dict, outplot: str):
    all_dfs = []
    for method, df in df_dict.items():
        tmp = df.copy()
        tmp['Method'] = method
        all_dfs.append(tmp)
    combined_df = pd.concat(all_dfs, axis=0)
    logger.info(f"Combined DataFrame shape: {combined_df.shape}")
    methods = combined_df['Method'].unique().tolist()
    cell_types = combined_df['Cell_type'].unique().tolist()
    palette = sns.color_palette("Set2", n_colors=len(cell_types))
    cell_color_map = dict(zip(cell_types, palette))

    plt.figure(figsize=(12, 8))
    ax = sns.violinplot(
        x='Method',
        y='LV1',
        hue='Cell_type',
        data=combined_df,
        palette=cell_color_map,
        inner='box',
        linewidth=1.5,
        dodge=True
    )

    # 计算每个方法+cell_type的小提琴中心点
    violin_positions = {}
    for method_idx, method in enumerate(methods):
        cells_in_method = combined_df[combined_df['Method']==method]['Cell_type'].unique()
        n_cell_types = len(cells_in_method)
        offsets = np.linspace(-0.2, 0.2, n_cell_types)
        for i, ct in enumerate(cells_in_method):
            violin_positions[(method, ct)] = method_idx + offsets[i]

    # 遍历每个方法做两两检验
    for method in methods:
        subset = combined_df[combined_df['Method'] == method]
        cells_in_method = subset['Cell_type'].unique().tolist()
        if len(cells_in_method) < 2:
            continue
        pairs = list(itertools.combinations(cells_in_method, 2))
        step = subset['LV1'].max() * 0.05

        for i, (c1, c2) in enumerate(pairs):
            data1 = subset[subset['Cell_type']==c1]['LV1']
            data2 = subset[subset['Cell_type']==c2]['LV1']
            stat, p = mannwhitneyu(data1, data2, alternative='two-sided')

            x1 = violin_positions[(method, c1)]
            x2 = violin_positions[(method, c2)]
            y_base = max(data1.max(), data2.max())
            y = y_base + (i+1)*step
            y_top = y + step*0.3  # 梯形顶部高度

            # 画梯形：竖直-水平-竖直
            ax.plot([x1, x1], [y, y_top], color='black', lw=2)
            ax.plot([x1, x2], [y_top, y_top], color='black', lw=2)
            ax.plot([x2, x2], [y_top, y], color='black', lw=2)

            # 显示星号
            ax.text((x1+x2)/2, y_top + step*0.05,
                    "****" if p<0.0001 else
                    "***" if p<0.001 else
                    "**" if p<0.01 else
                    "*" if p<0.05 else "ns",
                    ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Method')
    ax.set_ylabel('LV1 score')
    plt.tight_layout()
    plt.savefig(outplot, dpi=300)
    plt.close()
    print(f"图已保存到 {outplot}")

def combine_run():
    Latent = {
        "pca": "/home/luosg/Data/genomeStability/output/result/latent/fa_factors.csv",
        "nmf": "/home/luosg/Data/genomeStability/output/result/latent/nmf_factors.csv",
        "factor": "/home/luosg/Data/genomeStability/output/result/latent/pca_factors.csv"
    }
    # Latent = {
    #     "pca": "/home/luosg/Data/genomeStability/output/result1/latent/pca_factors.csv",
    #     "nmf": "/home/luosg/Data/genomeStability/output/result1/latent/nmf_factors.csv",
    #     "factor": "/home/luosg/Data/genomeStability/output/result1/latent/fa_factors.csv" 
    # }
    outdir="/home/luosg/Data/genomeStability/output/result/latent"
    # outdir = "/home/luosg/Data/genomeStability/output/result1/latent"
    df_dict = {}
    for method, infile in Latent.items():
        df = pd.read_csv(infile, index_col=0)
        df_annotation = pd.read_csv("/home/luosg/Data/genomeStability/data/target_fq.tsv", sep="\t")
        df_annotation.drop_duplicates(subset=["Sample_id"], inplace=True)
        Sample_to_Status_map = df_annotation.set_index('Sample_id')['Status']
        new_index = df.index.map(Sample_to_Status_map)
        df.index = new_index
        df_scaled = (df - df.min()) / (df.max() - df.min())
        df_scaled = df_scaled.reset_index(names="Cell_type")
        df_dict[method] = df_scaled

    outplot = f"{outdir}/all_methods_violin.png"
    violin_by_method_celltype_step(df_dict, outplot)
if __name__ == "__main__":
    # single_run()
    combine_run()
