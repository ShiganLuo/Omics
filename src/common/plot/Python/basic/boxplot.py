from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import seaborn as sns
def plot_comparison(
    result: Dict[str, Any],
    title: str = "Mapping Rate Comparison",
    xlabel: str = "Sample Type",
    ylabel: str = "Mapping Rate (%)",
    figsize: Tuple[int, int] = (10, 7),
    palette: Tuple[str, str] = ("#3498db", "#e74c3c"),
    show_points: bool = True,
    show_stats: bool = True,
    save_path: Optional[str] = None,
    dpi: int = 300,
    verbose: bool = True
) -> Figure:
    """
    Create a publication-quality comparison plot for two groups.

    Parameters
    ----------
    result : dict
        Result dictionary from prepare_data() containing:
        - group1_name, group2_name: Names of the groups
        - group1_data, group2_data: Data values for each group
        - comparison: Statistical comparison results
    title : str, optional
        Plot title (default: "Mapping Rate Comparison").
    xlabel : str, optional
        X-axis label (default: "Sample Type").
    ylabel : str, optional
        Y-axis label (default: "Mapping Rate (%)").
    figsize : tuple of int, optional
        Figure size (default: (10, 7)).
    palette : tuple of str, optional
        Colors for the two groups (default: ("#3498db", "#e74c3c")).
    show_points : bool, optional
        Whether to show individual data points (default: True).
    show_stats : bool, optional
        Whether to show statistics annotation (default: True).
    save_path : str, optional
        Path to save the figure (default: None, does not save).
    dpi : int, optional
        DPI for saved figure (default: 300).
    verbose : bool, optional
        Whether to logger.info information about the plot (default: True).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure object.

    Examples
    --------
    >>> from prepare_data import prepare_data
    >>> from plot_comparison import plot_comparison
    >>> result = prepare_data("RNAseq.tsv")
    >>> fig = plot_comparison(result, save_path="comparison.png")
    """
    # Configure Chinese font
    selected_font = configure_chinese_font()
    if verbose and selected_font:
        logger.info(f"Using Chinese font: {selected_font}")

    # Extract data from result
    group1_name = result['group1_name']
    group2_name = result['group2_name']
    group1_data = result['group1_data']
    group2_data = result['group2_data']
    comparison = result['comparison']

    # Set style
    sns.set_style("whitegrid")
    apply_font_after_style(selected_font)
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data for plotting
    data_to_plot = [group1_data, group2_data]
    positions = [1, 2]

    # Create box plot
    bp = ax.boxplot(
        data_to_plot,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showmeans=True,
        meanprops={"marker": "D", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": 8},
        medianprops={"color": "black", "linewidth": 2},
        whiskerprops={"color": "black", "linewidth": 1.5},
        capprops={"color": "black", "linewidth": 1.5},
        flierprops={"marker": "o", "markerfacecolor": "gray", "markersize": 6, "alpha": 0.5}
    )

    # Color the boxes
    for patch, color in zip(bp['boxes'], palette):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor("black")
        patch.set_linewidth(1.5)

    # Add individual data points if requested
    if show_points:
        for i, (data, pos, color) in enumerate(zip(data_to_plot, positions, palette)):
            # Add jitter to x-coordinates
            jitter = np.random.uniform(-0.15, 0.15, size=len(data))
            ax.scatter(
                np.full_like(data, pos) + jitter,
                data,
                color=color,
                alpha=0.6,
                s=50,
                edgecolors='black',
                linewidths=0.5,
                zorder=5
            )

    # Set x-axis labels with explicit font properties
    from matplotlib.font_manager import FontProperties
    if selected_font:
        font_prop = FontProperties(family=selected_font, size=12)
        font_prop_bold = FontProperties(family=selected_font, size=14, weight='bold')
        font_prop_title = FontProperties(family=selected_font, size=16, weight='bold')
    else:
        font_prop = FontProperties(size=12)
        font_prop_bold = FontProperties(size=14, weight='bold')
        font_prop_title = FontProperties(size=16, weight='bold')

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{group1_name}\n(n={len(group1_data)})", f"{group2_name}\n(n={len(group2_data)})"],
                       fontproperties=font_prop)

    # Set axis labels
    ax.set_xlabel(xlabel, fontproperties=font_prop_bold)
    ax.set_ylabel(ylabel, fontproperties=font_prop_bold)
    ax.set_title(title, fontproperties=font_prop_title, pad=20)

    # Add statistics annotation if requested
    if show_stats:
        # Get statistical results
        p_value = comparison['p_value']
        significant = comparison['significant']
        effect_size = comparison['effect_size']
        effect_interp = comparison['effect_size_interpretation']

        # Determine significance stars
        if p_value < 0.001:
            sig_stars = "***"
        elif p_value < 0.01:
            sig_stars = "**"
        elif p_value < 0.05:
            sig_stars = "*"
        else:
            sig_stars = "ns"

        # Calculate position for annotation
        y_max = max(max(group1_data), max(group2_data))
        y_min = min(min(group1_data), min(group2_data))
        y_range = y_max - y_min

        # Position the bracket
        bracket_y = y_max + y_range * 0.15
        star_y = bracket_y + y_range * 0.05

        # Draw bracket
        ax.plot([1, 1, 2, 2], [bracket_y, bracket_y + y_range * 0.02, bracket_y + y_range * 0.02, bracket_y],
                color='black', linewidth=1.5)

        # Add significance stars
        ax.text(1.5, star_y, sig_stars, ha='center', va='bottom', fontsize=16, fontweight='bold')

        # Add p-value text
        p_text = f"p = {p_value:.4f}" if p_value >= 0.001 else f"p < 0.001"
        ax.text(1.5, star_y - y_range * 0.08, p_text, ha='center', va='top', fontsize=11, style='italic')

        # Add effect size info - use axes transform for better positioning
        effect_text = f"{comparison['effect_size_name']}: {effect_size:.2f} ({effect_interp})"
        ax.text(0.98, 0.98, effect_text, transform=ax.transAxes, fontsize=10, 
                color='gray', ha='right', va='top', style='italic')

        # Add test name - use axes transform for better positioning
        test_text = f"Test: {comparison['test_name']}"
        ax.text(0.02, 0.02, test_text, transform=ax.transAxes, fontsize=9, color='gray', style='italic')

    # Add mean and SD annotations
    for i, (data, pos, name) in enumerate(zip(data_to_plot, positions, [group1_name, group2_name])):
        mean_val = np.mean(data)
        std_val = np.std(data, ddof=1)
        ax.annotate(
            f"Mean: {mean_val:.1f}%\nSD: {std_val:.1f}%",
            xy=(pos, np.median(data)),
            xytext=(pos + 0.4 if i == 0 else pos - 0.4, np.median(data)),
            fontsize=9,
            ha='left' if i == 0 else 'right',
            va='center',
            arrowprops=dict(arrowstyle='->', color='gray', lw=1),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8)
        )

    # Adjust y-axis limits to make room for annotation
    ax.set_ylim(y_min - y_range * 0.1, y_max + y_range * 0.35)

    # Add grid
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Tight layout
    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        if verbose:
            logger.info(f"Figure saved to: {save_path}")

    if verbose:
        logger.info(f"Plot created successfully")
        logger.info(f"  Groups: {group1_name} (n={len(group1_data)}) vs {group2_name} (n={len(group2_data)})")
        logger.info(f"  p-value: {p_value:.4f} ({'significant' if significant else 'not significant'})")

    return fig

