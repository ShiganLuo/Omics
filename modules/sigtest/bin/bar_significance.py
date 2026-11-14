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


# =========================
# bracket drawing (display-point layout)
# =========================
# Screen-unit geometry, identical on linear and log axes: the first
# bracket sits LEG_PT above the cluster's tallest bar top; each further
# stacked bracket adds BRACKET_STEP_PT; the star label sits TEXT_PT
# above its bracket line; FONT_ASCENT_PT approximates the label font's
# ascent for headroom reservation.
LEG_PT = 10.0
TEXT_PT = 4.0
FONT_ASCENT_PT = 8.0
BRACKET_STEP_PT = 14.0


def y_to_px(ax, y_data: float) -> float:
    """
    Map a data y value on ``ax`` to display coordinates (pixels).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes whose transform is used.
    y_data : float
        Data-coordinate y value.

    Returns
    -------
    float
        Display-coordinate y (pixels).
    """
    return float(ax.transData.transform((0.0, y_data))[1])


def px_to_y(ax, y_disp: float) -> float:
    """
    Map a display-coordinate y (pixels) to a data value on ``ax``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes whose transform is used.
    y_disp : float
        Display-coordinate y (pixels).

    Returns
    -------
    float
        Data-coordinate y value.
    """
    return float(ax.transData.inverted().transform((0.0, y_disp))[1])


def fit_bracket_headroom(ax, bases, n_brackets: int) -> None:
    """
    Raise the y upper limit so the tallest bracket stack fits.

    Finalize the axes layout (ticks, legend, tight_layout, canvas.draw)
    BEFORE calling this: display transforms must be stable afterwards.
    The new top is bisected so that the tallest base lands exactly the
    stack height below the axes-top pixel. Only ``ylim[1]`` moves, so
    this works identically on linear and log axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes already containing the bars and error bars.
    bases : array-like of float
        Per-cluster y of the tallest bar top including its error bar.
    n_brackets : int
        Maximum number of stacked brackets in any cluster.

    Returns
    -------
    None
    """
    if not n_brackets:
        return
    fig = ax.get_figure()
    px_per_pt = fig.dpi / 72.0
    stack_px = (LEG_PT + (n_brackets - 1) * BRACKET_STEP_PT + TEXT_PT + FONT_ASCENT_PT + 2.0) * px_per_pt
    positives = [float(b) for b in bases if b > 0]
    if not positives:
        return
    base = max(positives)
    axes_top_px = y_to_px(ax, ax.get_ylim()[1])  # top of the axes box
    target_px = axes_top_px - stack_px
    if y_to_px(ax, base) <= target_px:
        return  # already fits
    bottom = ax.get_ylim()[0]
    is_log = ax.get_yscale() == "log"
    lo = ax.get_ylim()[1]
    hi = lo * 10.0 if is_log else lo + (lo - bottom)
    while y_to_px(ax, base) > target_px:  # expand until it fits
        lo = hi
        hi = hi * 10.0 if is_log else hi + (hi - bottom)
        ax.set_ylim(bottom=bottom, top=hi)
        if hi > 1e12 * lo:
            break
    for _ in range(60):  # bisect the top value
        mid = np.sqrt(lo * hi) if is_log else (lo + hi) / 2.0
        ax.set_ylim(bottom=bottom, top=mid)
        if y_to_px(ax, base) > target_px:
            lo = mid
        else:
            hi = mid
    ax.set_ylim(bottom=bottom, top=hi)


def draw_significance_bracket(
    ax,
    x1: float,
    x2: float,
    y_base: float,
    label: str,
    level: int = 0,
    lw: float = 0.7,
    fontsize: float = 8.0,
    color: str = "#333333",
) -> float:
    """
    Draw one 宝盖头 significance bracket above a bar cluster.

    The bracket is placed ``LEG_PT + level * BRACKET_STEP_PT`` display
    points above ``y_base`` (screen units, so spacing is identical on
    linear and log axes). The horizontal line spans ``x1`` to ``x2``
    (bar centers) with small outward diagonal ticks at both ends, drawn
    as a single polyline to avoid line-cap overlap at junctions.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes (layout already finalized, see fit_bracket_headroom).
    x1, x2 : float
        Data x of the left/right bar centers.
    y_base : float
        Data y of the cluster's tallest bar top including its error bar.
    label : str
        Star label (from p_to_star or custom text).
    level : int, default=0
        Stack level (0 = closest to the bars).
    lw : float, default=0.7
        Bracket line width.
    fontsize : float, default=8.0
        Star label font size.
    color : str, default="#333333"
        Bracket and label color.

    Returns
    -------
    float
        Display y (pixels) of the star label top, for further stacking.
    """
    fig = ax.get_figure()
    px_per_pt = fig.dpi / 72.0
    y_hat_disp = y_to_px(ax, y_base) + (LEG_PT + level * BRACKET_STEP_PT) * px_per_pt
    y_hat = px_to_y(ax, y_hat_disp)
    y_text = px_to_y(ax, y_hat_disp + TEXT_PT * px_per_pt)
    y_tick = px_to_y(ax, y_hat_disp - 3.0 * px_per_pt)
    ax.plot([x1, x1, x2, x2], [y_tick, y_hat, y_hat, y_tick],
            lw=lw, color=color, clip_on=False)
    ax.text((x1 + x2) / 2, y_text, label,
            ha="center", va="bottom", fontsize=fontsize, color=color)
    return y_hat_disp + (TEXT_PT + FONT_ASCENT_PT) * px_per_pt


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
    """
        Plot grouped bar charts with a broken y-axis to compare SV type counts across conditions.

        This function visualizes structural variant (SV) counts for multiple groups using
        side-by-side bar plots. A broken y-axis is applied to simultaneously display both
        low and high count ranges. Statistical significance between groups for each SV type
        is evaluated using a chi-square test and annotated as star labels.

        Parameters
        ----------
        df : pandas.DataFrame
            Input DataFrame containing SV counts.
        out_png : str
            Output file path for the saved plot (PNG format).
        group_col : str, default="group"
            Column name indicating group labels (e.g., Control, Experiment).
        svtype_col : str, default="svtype"
            Column name indicating SV types.
        count_col : str, default="count"
            Column name containing count values.
        group_order : tuple of str, default=("Control", "Experiment")
            Order of groups for plotting and comparison.
        svtype_order : tuple of str, default=("BND", "DEL", "DUP", "INS", "INV")
            Order of SV types on the x-axis.
        legend_map : dict, optional
            Mapping from group names to legend labels. If None, uses group names directly.
        figsize : tuple, default=(9, 5)
            Figure size in inches.
        ylabel : str, default="SV count"
            Label for the y-axis.
        dpi : int, default=300
            Resolution of the output image.

        Notes
        -----
        - Missing combinations of (SV type, group) are filled with zero.
        - Significance is computed per SV type using chi-square contingency test.
        - Significance levels:
            * p < 1e-4 : "****"
            * p < 1e-3 : "***"
            * p < 1e-2 : "**"
            * p < 0.05 : "*"
            * otherwise: "ns"
        - The broken y-axis improves visualization when count ranges are highly skewed.

        Returns
        -------
        None
            The plot is saved to `out_png`.
    """
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

if __name__ == "__main__":
    result_df = pd.read_csv("/data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint/RNAseq/deseq2_oncoprint_fold_change_all.csv", index=False)
    result_df = result_df.dropna()
    condition_names = ("P6", "P20")
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