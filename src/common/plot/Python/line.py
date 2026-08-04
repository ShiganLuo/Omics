# Author: ShiganLuo
# Email: 25303020102@qq.com
# GitHub: https://github.com/ShiganLuo
import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from typing import List, Optional, Dict
from itertools import cycle
import os
from pathlib import Path
import colorsys
import textwrap
from typing import Iterable, Literal, Sequence, Tuple, cast, Union
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline, PchipInterpolator, UnivariateSpline
from scipy.signal import savgol_filter
import pandas as pd
import logging
logger = logging.getLogger(__name__)
SmoothMethod = Literal["none", "linear", "cubic", "pchip", "spline", "asymptotic"]
FloatArray = NDArray[np.float64]



def _to_float_array(values: Iterable[float], name: str) -> FloatArray:
    """Convert an input sequence to a 1D float NumPy array.

    Parameters
    ----------
    values : Iterable[float]
        Numeric input sequence.
    name : str
        Parameter name used in error messages.

    Returns
    -------
    np.ndarray
        One-dimensional floating-point array.

    Raises
    ------
    ValueError
        Raised when the input is empty, not 1D, or contains NaN/Inf.
    """
    arr = np.asarray(list(values), dtype=np.float64)
    if arr.size == 0:
        raise ValueError(f"`{name}` must not be empty")
    if arr.ndim != 1:
        raise ValueError(f"`{name}` must be a one-dimensional sequence")
    if not np.isfinite(arr).all():
        raise ValueError(f"`{name}` must not contain NaN or Inf")
    return arr


def _sort_and_merge_duplicate_x(x: FloatArray, y: FloatArray) -> Tuple[FloatArray, FloatArray]:
    """Sort by x in ascending order and merge duplicate x values by mean y.

    Parameters
    ----------
    x : FloatArray
        X-axis values.
    y : FloatArray
        Y-axis values.

    Returns
    -------
    Tuple[FloatArray, FloatArray]
        ``(x_sorted_unique, y_merged)`` arrays with duplicate x entries
        collapsed to their mean y.
    """
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]

    uniq_x = []
    merged_y = []
    start = 0
    for idx in range(1, len(x_sorted) + 1):
        if idx == len(x_sorted) or x_sorted[idx] != x_sorted[start]:
            uniq_x.append(float(x_sorted[start]))
            merged_y.append(float(np.mean(y_sorted[start:idx])))
            start = idx

    return np.asarray(uniq_x, dtype=np.float64), np.asarray(merged_y, dtype=np.float64)


def smooth_depth_sensitivity_curve(
    depth: Iterable[float],
    sensitivity: Iterable[float],
    method: SmoothMethod = "pchip",
    num_points: int = 300,
    spline_s: float = 0.0,
) -> Tuple[FloatArray, FloatArray]:
    """Smooth or interpolate discrete depth-sensitivity points.

    Preprocessing steps:
    1) validate input, 2) sort by depth, 3) merge duplicate depth values
    by averaging sensitivity.

    Parameters
    ----------
    depth : Sequence[float]
        X-axis values (sequencing depth).
    sensitivity : Sequence[float]
        Y-axis values (sensitivity).
    method : {'none', 'linear', 'cubic', 'pchip', 'spline', 'asymptotic'}, default='pchip'
        Curve generation method:
        - 'none'  : no smoothing, return sorted raw points only;
        - 'linear': linear interpolation;
        - 'cubic' : cubic spline (falls back to linear if < 4 points);
        - 'pchip' : shape-preserving piecewise cubic interpolation;
        - 'spline': smoothing spline controlled by `spline_s`.
        - 'asymptotic': log(1-y) transform + cubic spline, natural saturation at y=1.
    num_points : int, default=300
        Number of sampled points for the smoothed curve.
    spline_s : float, default=0.0
        Smoothing factor used only when `method='spline'`.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        `(depth_smooth, sensitivity_smooth)` arrays.

    Raises
    ------
    ValueError
        Raised when lengths mismatch, points are insufficient, or
        parameters are invalid.
    """
    x = _to_float_array(depth, "depth")
    y = _to_float_array(sensitivity, "sensitivity")

    if len(x) != len(y):
        raise ValueError("`depth` and `sensitivity` must have the same length")

    x, y = _sort_and_merge_duplicate_x(x, y)

    if len(x) < 2:
        raise ValueError("At least 2 distinct depth points are required")

    if method not in {"none", "linear", "cubic", "pchip", "spline", "asymptotic"}:
        raise ValueError(f"Unsupported method: {method}")
    if num_points < 2:
        raise ValueError("`num_points` must be >= 2")

    if method == "none":
        return np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)

    grid_n = max(int(num_points), len(x))
    x_new: FloatArray = np.asarray(
        np.linspace(float(x.min()), float(x.max()), grid_n), dtype=np.float64
    )

    if method == "linear":
        y_new: FloatArray = np.asarray(np.interp(x_new, x, y), dtype=np.float64)
    elif method == "cubic":
        if len(x) < 4:
            y_new = np.asarray(np.interp(x_new, x, y), dtype=np.float64)
        else:
            y_new = np.asarray(CubicSpline(x, y)(x_new), dtype=np.float64)
    elif method == "pchip":
        y_new = np.asarray(PchipInterpolator(x, y)(x_new), dtype=np.float64)
    elif method == "spline":
        degree = min(3, len(x) - 1)
        y_new = cast(
            FloatArray,
            np.asarray(UnivariateSpline(x, y, s=spline_s, k=degree)(x_new), dtype=np.float64),
        )
    else:  # method == 'asymptotic'
        # log(1 - y) transform: maps asymptote y→1 to -∞, making interpolation smooth.
        eps = 1e-10
        y_clamped = np.clip(y, eps, 1.0 - eps)
        z = -np.log(1.0 - y_clamped)           # z ∈ [0, +∞), monotonically increasing with y
        if len(x) < 4:
            z_new = np.interp(x_new, x, z)
        else:
            z_new = PchipInterpolator(x, z)(x_new)
        y_new = np.asarray(1.0 - np.exp(-z_new), dtype=np.float64)
        y_new = np.clip(y_new, 0.0, 1.0)

    return x_new, y_new


def plot_depth_sensitivity_curve(
    depth: Iterable[float],
    sensitivity: Iterable[float],
    Sensitivity_threshold: float = 0.95,
    method: SmoothMethod = "pchip",
    num_points: int = 300,
    spline_s: float = 0.0,
    x_label: str = "Sequencing Depth",
    y_label: str = "Sensitivity",
    title: str = "Depth-Sensitivity Curve",
    output_path: Union[str, Path, None] = None,
    dpi: int = 300,
    show: bool = False,
):
    """Plot depth versus sensitivity with optional smoothing.

    Parameters
    ----------
    depth : Sequence[float]
        X-axis values.
    sensitivity : Sequence[float]
        Y-axis values.
    method : {'none', 'linear', 'cubic', 'pchip', 'spline'}, default='pchip'
        Smoothing/interpolation method. See
        `smooth_depth_sensitivity_curve`.
    num_points : int, default=300
        Number of sampled points for the smoothed curve.
    spline_s : float, default=0.0
        Smoothing spline factor, used only when `method='spline'`.
    x_label : str, default='Sequencing Depth'
        X-axis label.
    y_label : str, default='Sensitivity'
        Y-axis label.
    title : str, default='Depth-Sensitivity Curve'
        Figure title.
    output_path : Union[str, Path, None], default=None
        Output image path. If None, no file is saved.
    dpi : int, default=180
        Figure resolution when saving.
    show : bool, default=False
        Whether to display the figure window.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object for further customization by callers.

    Notes
    -----
    - Raw input points are shown as scatter points.
    - The processed result is shown as a line curve.
    - This function is API-oriented and can be reused in other scripts.
    """

    outdir = os.path.dirname(output_path) if output_path is not None else None
    os.makedirs(outdir, exist_ok=True) if outdir else None
    x_raw = _to_float_array(depth, "depth")
    y_raw = _to_float_array(sensitivity, "sensitivity")

    x_smooth, y_smooth = smooth_depth_sensitivity_curve(
        x_raw,
        y_raw,
        method=method,
        num_points=num_points,
        spline_s=spline_s,
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.scatter(x_raw, y_raw, color="#1f77b4", alpha=0.85, s=22, edgecolors="#1f77b4", linewidths=0.3, label="Raw points")

    if method == "none":
        order = np.argsort(x_smooth)
        ax.plot(
            x_smooth[order],
            y_smooth[order],
            color="#ff7f0e",
            linewidth=2.0,
            label="Connected line",
        )
    else:
        ax.plot(
            x_smooth,
            y_smooth,
            color="#ff7f0e",
            linewidth=2.2,
            label=f"Smoothed curve ({method})",
        )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(alpha=0.3, linestyle="--")

    threshold_idx = int(np.searchsorted(y_smooth, Sensitivity_threshold, side="left"))
    threshold_idx = min(max(threshold_idx, 0), len(x_smooth) - 1)
    threshold_x = float(x_smooth[threshold_idx])

    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    y_bottom, _ = ax.get_ylim()

    ax.vlines(
        threshold_x,
        y_bottom,
        Sensitivity_threshold,
        colors=["red"],
        linestyles="dashed",
        label=f"threshold={Sensitivity_threshold:.2f}",
    )
    ax.scatter([threshold_x], [Sensitivity_threshold], color="red", s=26, zorder=4, edgecolors="red", linewidths=0.3)

    # Leader line: always show threshold depth below x-axis
    ax.annotate(
        f"{threshold_x:.0f}x",
        xy=(threshold_x, y_bottom),
        xytext=(8, -28),
        textcoords="offset points",
        ha="left", va="top", fontsize=8, fontweight="bold", color="red",
        arrowprops=dict(arrowstyle="-", color="red", lw=0.8,
                        connectionstyle="angle,angleA=-90,angleB=180,rad=0.2"),
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="red", lw=0.5, alpha=0.9),
    )

    # Method annotation
    if method != "none":
        ax.text(0.98, 0.02, f"fit: {method}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="gray",
                style="italic", alpha=0.7)

    yticks = np.unique(np.append(ax.get_yticks(), Sensitivity_threshold))
    ytick_candidates = np.sort(yticks[(yticks >= y_bottom) & (yticks <= 1.0)])
    if ytick_candidates.size == 0:
        ytick_candidates = np.asarray([1.0], dtype=np.float64)
    ax.set_yticks(ytick_candidates)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=True)
    fig.tight_layout(rect=(0.0, 0.08, 1.0, 1.0))

    if output_path is not None:
        fig.savefig(str(Path(output_path)), dpi=dpi)

    if show:
        plt.show()

    return fig


def plot_multiple_depth_sensitivity_curves(
    curves: Sequence[Tuple[Iterable[float], Iterable[float]]],
    sensitivity_threshold: float = 0.95,
    method: SmoothMethod = "pchip",
    num_points: int = 300,
    spline_s: float = 0.0,
    labels: Union[Sequence[str], None] = None,
    x_label: str = "Sequencing Depth",
    y_label: str = "Sensitivity",
    title: str = "",
    output_path: Union[str, Path, None] = None,
    dpi: int = 300,
    show: bool = False,
    x_tick_rotation: float = 30.0,
    is_smoothed: bool = True,
    emphasize_steep_region: bool = True,
    x_log_scale: bool = False
):
    """Plot multiple depth-sensitivity curves on the same axes.

    Overlays several curves for visual comparison, with optional x-axis
    warping to emphasize steep regions, log-scale display, and per-curve
    threshold annotations.

    Parameters
    ----------
    curves : Sequence[Tuple[Iterable[float], Iterable[float]]]
        Sequence of ``(depth, sensitivity)`` pairs, one per curve.
    sensitivity_threshold : float, default=0.95
        Horizontal reference line and threshold for effective depth.
    method : {'none', 'linear', 'cubic', 'pchip', 'spline', 'asymptotic'}, default='pchip'
        Smoothing/interpolation method. See
        ``smooth_depth_sensitivity_curve`` for details.
    num_points : int, default=300
        Number of sampled points for each smoothed curve.
    spline_s : float, default=0.0
        Smoothing factor, used only when ``method='spline'``.
    labels : Sequence[str] | None, default=None
        Display name for each curve. Auto-generated if ``None``.
    x_label : str, default='Sequencing Depth'
        X-axis label.
    y_label : str, default='Sensitivity'
        Y-axis label.
    title : str, default=''
        Figure title. Omitted when empty.
    output_path : str | Path | None, default=None
        Output image path. If ``None``, no file is saved.
    dpi : int, default=300
        Figure resolution when saving.
    show : bool, default=False
        Whether to display the figure window.
    x_tick_rotation : float, default=30.0
        Rotation angle of x-axis tick labels in degrees.
    is_smoothed : bool, default=True
        Whether to apply smoothing. When ``False``, raw points are
        connected directly.
    emphasize_steep_region : bool, default=True
        Widen high-slope regions on the x-axis for visual emphasis.
        Ignored when ``x_log_scale=True``.
    x_log_scale : bool, default=False
        Display x-axis in logarithmic scale.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object for further customization by callers.

    Raises
    ------
    ValueError
        When fewer than 2 curves are provided or ``labels`` length
        does not match the number of curves.
    """
    if len(curves) < 2:
        raise ValueError("At least two curves are required for overlay plotting")

    if labels is None:
        labels = [f"Curve {idx + 1}" for idx in range(len(curves))]
    if len(labels) != len(curves):
        raise ValueError("`labels` length must match the number of curves")

    fig, ax = plt.subplots(figsize=(12, 6))
    line_colors = [
        colorsys.hsv_to_rgb(i / max(len(curves), 1), 0.75, 0.85)
        for i in range(len(curves))
    ]

    raw_curves: list[Tuple[FloatArray, FloatArray]] = []
    smooth_curves: list[Tuple[FloatArray, FloatArray]] = []
    threshold_x_values: list[float] = []
    threshold_colors: list[tuple[float, float, float]] = []

    for idx, (depth_values, sens_values) in enumerate(curves):
        x_raw = _to_float_array(depth_values, f"depth_{idx + 1}")
        y_raw = _to_float_array(sens_values, f"sensitivity_{idx + 1}")
        if is_smoothed:
            x_smooth, y_smooth = smooth_depth_sensitivity_curve(
                x_raw,
                y_raw,
                method=method,
                num_points=num_points,
                spline_s=spline_s,
            )
        else:
            x_smooth, y_smooth = x_raw, y_raw
        raw_curves.append((x_raw, y_raw))
        smooth_curves.append((x_smooth, y_smooth))

        threshold_idx = int(np.searchsorted(y_smooth, sensitivity_threshold, side="left"))
        threshold_idx = min(max(threshold_idx, 0), len(x_smooth) - 1)
        threshold_x_values.append(float(x_smooth[threshold_idx]))
        threshold_colors.append(line_colors[idx % len(line_colors)])

    x_min = min(float(xs.min()) for xs, _ in smooth_curves)
    x_max = max(float(xs.max()) for xs, _ in smooth_curves)
    x_grid = np.linspace(x_min, x_max, 800, dtype=np.float64)

    slope_stack: list[FloatArray] = []
    for xs, ys in smooth_curves:
        yi = np.interp(x_grid, xs, ys)
        slope_stack.append(np.asarray(np.abs(np.gradient(yi, x_grid)), dtype=np.float64))

    mean_slope = np.mean(np.vstack(slope_stack), axis=0)
    slope_max = float(np.max(mean_slope))
    slope_norm = mean_slope / slope_max if slope_max > 0 else np.zeros_like(mean_slope)

    if x_log_scale:
        x_warp = x_grid.copy()

        def map_x(values: FloatArray) -> FloatArray:
            return np.asarray(values, dtype=np.float64)

        warped_threshold_x = [float(x) for x in threshold_x_values]
    else:
        steep_emphasis = 4.0 if emphasize_steep_region else 0.0
        weight = 1.0 + steep_emphasis * slope_norm
        x_warp = np.zeros_like(x_grid)
        x_warp[1:] = np.cumsum(np.diff(x_grid) * 0.5 * (weight[:-1] + weight[1:]))

        def map_x(values: FloatArray) -> FloatArray:
            return np.asarray(np.interp(values, x_grid, x_warp), dtype=np.float64)

        warped_threshold_x = [float(np.interp(x, x_grid, x_warp)) for x in threshold_x_values]

    for idx, ((x_raw, y_raw), (x_smooth, y_smooth)) in enumerate(zip(raw_curves, smooth_curves)):
        color = line_colors[idx % len(line_colors)]
        ax.scatter(map_x(x_raw), y_raw, color=color, alpha=0.75, s=22, edgecolors=color, linewidths=0.3)
        ax.plot(
            map_x(x_smooth),
            y_smooth,
            color=color,
            linewidth=2.2,
            label=f"{labels[idx]} ({threshold_x_values[idx]:.0f}x)",
        )

    y_all = np.concatenate([ys for _, ys in raw_curves])
    y_min = float(np.min(y_all))
    y_max = float(np.max(y_all))
    y_span = max(y_max - y_min, 1e-6)
    y_bottom = y_min - 0.12 * y_span
    positive_x = np.concatenate([xs[xs > 0] for xs, _ in smooth_curves])
    min_positive_x = float(np.min(positive_x)) if positive_x.size > 0 else 1.0
    if x_log_scale:
        ax.set_xscale("log")
        ax.set_xlim(left=min_positive_x, right=x_max)
    else:
        ax.set_xlim(left=max(0.0, float(x_warp.min())))
    ax.set_ylim(bottom=y_bottom)
    ax.grid(alpha=0.3, linestyle="--")

    for idx, x_threshold in enumerate(warped_threshold_x):
        ax.vlines(
            x_threshold,
            y_bottom,
            sensitivity_threshold,
            colors=[threshold_colors[idx]],
            linestyles="dashed",
            alpha=0.85,
        )
        ax.scatter([x_threshold], [sensitivity_threshold], color=threshold_colors[idx], s=24, zorder=4)

    ax.axhline(
        sensitivity_threshold,
        color="gray",
        linestyle=":",
        linewidth=1.1,
        alpha=0.9,
        label="_nolegend_",
    )
    ax.axhline(
        1.0,
        color="gray",
        linestyle=":",
        linewidth=1.1,
        alpha=0.9,
        label="_nolegend_",
    )

    # X-axis ticks: base positions only, threshold values shown via leader lines
    if x_log_scale:
        ax.xaxis.set_major_locator(mticker.LogLocator(base=10.0))
        ax.xaxis.set_minor_locator(mticker.NullLocator())
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda value, _: f"{value:g}x")
        )
    else:
        base_x_labels = np.linspace(x_min, x_max, 6)
        base_x_positions = np.asarray(np.interp(base_x_labels, x_grid, x_warp), dtype=np.float64)
        ax.set_xticks(base_x_positions)
        ax.set_xticklabels([f"{label:.0f}x" for label in base_x_labels], rotation=x_tick_rotation, ha="right")

    # Leader lines: always show threshold depth below x-axis
    for idx, x_threshold in enumerate(warped_threshold_x):
        ax.annotate(
            f"{threshold_x_values[idx]:.0f}x",
            xy=(x_threshold, y_bottom),
            xytext=(8, -28),
            textcoords="offset points",
            ha="left", va="top", fontsize=8, fontweight="bold",
            color=threshold_colors[idx],
            arrowprops=dict(arrowstyle="-", color=threshold_colors[idx], lw=0.8,
                            connectionstyle="angle,angleA=-90,angleB=180,rad=0.2"),
            bbox=dict(boxstyle="round,pad=0.15", fc="white",
                      ec=threshold_colors[idx], lw=0.5, alpha=0.9),
        )

    # Method annotation
    if method != "none":
        ax.text(0.98, 0.02, f"fit: {method}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="gray",
                style="italic", alpha=0.7)

    yticks = np.asarray(ax.get_yticks(), dtype=float)
    # Remove existing ticks too close to threshold or 1.0 to avoid label overlap
    yticks = yticks[(np.abs(yticks - sensitivity_threshold) > 0.015) & (np.abs(yticks - 1.0) > 0.015)]
    yticks = np.unique(np.append(yticks, [sensitivity_threshold, 1.0]))
    ytick_candidates = np.sort(yticks[(yticks >= y_bottom) & (yticks <= 1.05)])
    if ytick_candidates.size == 0:
        ytick_candidates = np.asarray([1.0], dtype=np.float64)
    ax.set_yticks(ytick_candidates)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    legend_handles, legend_labels = ax.get_legend_handles_labels()
    wrapped_legend_labels = [
        textwrap.fill(lbl, width=24, break_long_words=False, break_on_hyphens=False)
        for lbl in legend_labels
    ]

    ax.legend(
        legend_handles,
        wrapped_legend_labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=7.5,
        handlelength=1.4,
        handletextpad=0.4,
        columnspacing=0.8,
        borderpad=0.25,
        labelspacing=0.25,
        borderaxespad=0.0,
    )
    fig.tight_layout(rect=(0.0, 0.0, 0.76, 1.0))

    if output_path is not None:
        fig.savefig(str(Path(output_path)), dpi=dpi, bbox_inches="tight", pad_inches=0.1)

    if show:
        plt.show()

    return fig


def _effective_depth_at_sensitivity_threshold(
    depth: Iterable[float],
    sensitivity: Iterable[float],
    sensitivity_threshold: float = 0.95,
    is_smoothed: bool = True,
    method: SmoothMethod = "pchip",
    num_points: int = 300,
    spline_s: float = 0.0,
) -> float:
    """Get effective depth where sensitivity first reaches the threshold.

    Optionally smooth the input curve before searching for the first
    depth at which sensitivity meets or exceeds the threshold.

    Parameters
    ----------
    depth : Iterable[float]
        X-axis values (sequencing depth).
    sensitivity : Iterable[float]
        Y-axis values (sensitivity).
    sensitivity_threshold : float, default=0.95
        Target sensitivity level.
    is_smoothed : bool, default=True
        Whether to apply smoothing before threshold search.
    method : {'none', 'linear', 'cubic', 'pchip', 'spline', 'asymptotic'}, default='pchip'
        Smoothing method. Ignored when ``is_smoothed=False``.
    num_points : int, default=300
        Number of interpolation sample points.
    spline_s : float, default=0.0
        Smoothing factor, used only when ``method='spline'``.

    Returns
    -------
    float
        Effective sequencing depth at which sensitivity first reaches
        the threshold.  Falls back to the maximum depth if the threshold
        is never reached.
    """
    x_raw = _to_float_array(depth, "depth")
    y_raw = _to_float_array(sensitivity, "sensitivity")

    if is_smoothed:
        x_eval, y_eval = smooth_depth_sensitivity_curve(
            x_raw,
            y_raw,
            method=method,
            num_points=num_points,
            spline_s=spline_s,
        )
    else:
        x_eval, y_eval = _sort_and_merge_duplicate_x(x_raw, y_raw)

    hits = np.where(y_eval >= sensitivity_threshold)[0]
    if hits.size == 0:
        return float(x_eval[-1])
    return float(x_eval[int(hits[0])])


def plot_threshold_depth_vs_ng(
    ng_values: Iterable[float],
    effective_depth_values: Iterable[float],
    sensitivity_threshold: float = 0.95,
    x_label: str = "Input amount (ng)",
    y_label: str = "Effective sequencing depth (x)",
    title: str = "Effective depth vs input amount",
    output_path: str | Path | None = None,
    dpi: int = 300,
    show: bool = False,
):
    """Plot relationship between input DNA amount and effective depth.

    Shows how the required sequencing depth (to reach the sensitivity
    threshold) changes with different input DNA amounts.

    Parameters
    ----------
    ng_values : Iterable[float]
        Input DNA amounts in nanograms.
    effective_depth_values : Iterable[float]
        Effective sequencing depth at each input amount.
    sensitivity_threshold : float, default=0.95
        Sensitivity threshold used to compute effective depth.
    x_label : str, default='Input amount (ng)'
        X-axis label.
    y_label : str, default='Effective sequencing depth (x)'
        Y-axis label.
    title : str, default='Effective depth vs input amount'
        Figure title. Auto-generated when empty.
    output_path : str | Path | None, default=None
        Output image path. If ``None``, no file is saved.
    dpi : int, default=300
        Figure resolution when saving.
    show : bool, default=False
        Whether to display the figure window.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object for further customization by callers.

    Raises
    ------
    ValueError
        When lengths mismatch or fewer than 2 data points are provided.
    """
    ng = _to_float_array(ng_values, "ng_values")
    depth_eff = _to_float_array(effective_depth_values, "effective_depth_values")

    if len(ng) != len(depth_eff):
        raise ValueError("`ng_values` and `effective_depth_values` must have same length")
    if len(ng) < 2:
        raise ValueError("At least two points are required to draw ng-depth relation")

    ng_sorted, depth_sorted = _sort_and_merge_duplicate_x(ng, depth_eff)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(ng_sorted, depth_sorted, color="#1f77b4", s=34, label="Data points")
    ax.plot(ng_sorted, depth_sorted, color="#1f77b4", linewidth=2.0)

    if len(ng_sorted) >= 3:
        ng_dense = np.linspace(float(ng_sorted.min()), float(ng_sorted.max()), 300)
        depth_dense = PchipInterpolator(ng_sorted, depth_sorted)(ng_dense)
        ax.plot(ng_dense, depth_dense, color="#ff7f0e", linewidth=2.2, label="Smoothed trend")

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if not title:
        title = f"Effective depth at sensitivity threshold {sensitivity_threshold:.2f}"
    ax.set_title(title)
    ax.grid(alpha=0.3, linestyle="--")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.set_yticks(np.sort(np.asarray(ax.get_yticks(), dtype=float)[np.asarray(ax.get_yticks(), dtype=float) >= 0]))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=True)
    fig.tight_layout(rect=(0.0, 0.08, 1.0, 1.0))

    if output_path is not None:
        fig.savefig(str(Path(output_path)), dpi=dpi)
    if show:
        plt.show()

    return fig


def read_xy_from_csv(
    file_path: str | Path,
    x_col: str = "depth",
    y_col: str = "sensitivity",
    delimiter: str = ",",
) -> Tuple[np.ndarray, np.ndarray]:
    """Read X/Y data from a CSV or TSV file.

    Parameters
    ----------
    file_path : str | Path
        Input file path.
    x_col : str, default='depth'
        Column name for X values.
    y_col : str, default='sensitivity'
        Column name for Y values.
    delimiter : str, default=','
        Delimiter character. Use '\\t' for TSV.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        `(x, y)` arrays.

    Raises
    ------
    ValueError
        Raised when header is missing, columns are not found,
        or valid data is empty.
    """
    import csv

    resolved_delimiter = "\t" if delimiter == "\\t" else delimiter
    path = Path(file_path)

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=resolved_delimiter)
        if not reader.fieldnames:
            raise ValueError("Input file is missing a header row")
        if x_col not in reader.fieldnames or y_col not in reader.fieldnames:
            raise ValueError(
                f"Columns not found. available={reader.fieldnames}, required={x_col},{y_col}"
            )

        xs = []
        ys = []
        for row in reader:
            x_value = row.get(x_col)
            y_value = row.get(y_col)
            if x_value is None or y_value is None or x_value == "" or y_value == "":
                continue
            xs.append(float(x_value))
            ys.append(float(y_value))

    if not xs:
        raise ValueError("No valid data points were read")

    return np.asarray(xs, dtype=np.float64), np.asarray(ys, dtype=np.float64)



def run_gsva(
        count_matrix:str,
        gmt_pathway:str,
        outdir:str
):
    es = gp.gsva(
    data = count_matrix,
    gene_sets = gmt_pathway,
    outdir = outdir
    )
    return es

def plot_gsea_from_csv(
    csv_path: str,
    ranked_genes: List[str],
    out_png: str,
    top_n: Optional[int] = None,
    fdr_cutoff: float = 0.05,
    label_col: str = "Term",
    color_map: Optional[Dict[str, str]] = None,
    label_font: int = 9,
    title_font: int = 14,
    fig_size: tuple = (12, 7),
    legend_bottom: float = -0.18
):
    """
    绘制 GSEA 显著通路累积 ES 曲线（多彩、平滑、底部图例），图例列数根据通路数量和图宽自动调整。
    """

    # ===== 读取 CSV =====
    df = pd.read_csv(csv_path)

    # 定位 FDR 列
    fdr_col_candidates = [c for c in df.columns if "fdr" in c.lower()]
    if not fdr_col_candidates:
        raise ValueError("CSV 中未找到 FDR 列(例如 'FDR q-val')")
    fdr_col = fdr_col_candidates[0]

    # 只保留显著通路
    df_sig = df[df[fdr_col] < fdr_cutoff].copy()
    if df_sig.empty:
        print(f"没有显著通路 (FDR < {fdr_cutoff})，不绘图。")
        return

    df_sig = df_sig.sort_values(fdr_col)
    if top_n is not None:
        df_sig = df_sig.head(top_n)

    # ===== 绘图设置 =====
    plt.figure(figsize=fig_size)
    ax = plt.gca()
    N = len(ranked_genes)
    ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.6)

    # ===== 颜色循环 =====
    if color_map is None:
        color_map = {}
        colors_cycle = cycle(plt.get_cmap("tab20").colors)
    else:
        colors_cycle = None

    # ===== 绘制每条通路 =====
    for _, row in df_sig.iterrows():
        term = str(row[label_col])
        nes = row["NES"]
        lead_genes = str(row["Lead_genes"]).split(";")

        hits = np.array([1 if g in lead_genes else 0 for g in ranked_genes])
        nh = hits.sum()
        if nh == 0:
            continue
        no = N - nh
        running_es = np.cumsum(hits / nh - (1 - hits) / no)

        # 平滑
        win = min(len(running_es) - (1 - len(running_es) % 2), 101)
        if win >= 11:
            smooth_es = savgol_filter(running_es, window_length=win, polyorder=3)
        else:
            smooth_es = running_es

        # 获取颜色
        if term not in color_map:
            color = next(colors_cycle)
            color_map[term] = color
        else:
            color = color_map[term]

        ax.plot(smooth_es, linewidth=2, color=color, label=f"{term} (NES={nes:.2f})")

    # ===== 美化 =====
    ax.set_xlabel("Ranked Genes", fontsize=label_font)
    ax.set_ylabel("Enrichment Score (ES)", fontsize=label_font)
    # ax.set_title(f"GSEA Significant Pathways (FDR < {fdr_cutoff})", fontsize=title_font)
    ax.grid(True, alpha=0.3)

    # ===== 自动计算图例列数 =====
    n_terms = len(df_sig)
    # 根据图形宽度、字体大小、通路数量自动估算列数
    approx_char_per_col = 25  # 每列大约能放多少字符
    fig_width_inch = fig_size[0]
    legend_ncol = max(1, min(n_terms, int((fig_width_inch * 5) / approx_char_per_col)))

    # 图例底部
    ax.legend(
        fontsize=label_font,
        loc="upper center",
        bbox_to_anchor=(0.5, legend_bottom),
        ncol=legend_ncol,
        frameon=False
    )

    # 自动调整 subplots bottom
    plt.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=max(0.2, -legend_bottom + 0.05))
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"图保存到: {out_png}")
    return color_map


def plot_biax_line(
        local_alignment_file: str,
        global_alignment_file: str,
        out_png: str = "alignment_comparison_line.png"
    ):
    """
    Plot alignment comparison using dual y-axis line chart.
    """

    # 读取数据
    local_df = pd.read_csv(local_alignment_file, sep="\t")
    global_df = pd.read_csv(global_alignment_file, sep="\t")

    # 合并数据
    df = pd.merge(local_df, global_df, on="AminoAcid", suffixes=("_Local", "_Global"))

    # 排序（可选）
    df = df.sort_values("Similarity_Global", ascending=False)

    # 创建画布
    fig, ax1 = plt.subplots(figsize=(10,5))

    # 创建第二个y轴
    ax2 = ax1.twinx()

    # Local 折线
    line1 = ax1.plot(
        df["AminoAcid"],
        df["Similarity_Local"],
        marker="o",
        color="steelblue",
        label="Local Alignment"
    )

    # Global 折线
    line2 = ax2.plot(
        df["AminoAcid"],
        df["Similarity_Global"],
        marker="s",
        color="orange",
        label="Global Alignment"
    )

    # 轴标签
    ax1.set_ylabel("Max Score", color="steelblue")
    ax1.tick_params(axis='y', colors="steelblue")
    ax1.spines['left'].set_color("steelblue")

    ax2.set_ylabel("Global Similarity", color="orange")
    ax2.tick_params(axis='y', colors="orange")
    ax2.spines['right'].set_color("orange")
    ax2.spines['left'].set_color("none")  # 隐藏左侧边框

    # x轴旋转
    plt.xticks(rotation=45)

    # 合并图例
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper right")

    # 保存
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
def cor_plot(
	infile: str,
	cor_col: List[Tuple[str,str]],
	outprefix:str,
	font_size: float = 14.0
):
	"""Plot correlation scatter charts for one or multiple column pairs.

	Correlation and trendline method
	-------------------------------
	For each ``(x_col, y_col)`` pair:
	1. Both columns are converted to numeric with ``errors='coerce'``.
	2. Rows with NaN in either column are removed (pairwise complete cases).
	3. Pearson correlation coefficient is computed by
	   ``plot_df[x_col].corr(plot_df[y_col])``.
	4. A first-order linear trend line is fitted via least squares using
	   ``numpy.polyfit(x, y, 1)`` and overlaid on the scatter plot.

	Parameters
	----------
	infile : str
		Input TSV/CSV file path.
	cor_col : list[tuple[str, str]]
		Column pairs to plot, e.g. [("x1", "y1"), ("x2", "y2")].
	outprefix : str
		Output file prefix; final image is "{outprefix}_cor_scatter.png".
	font_size : float, default=14.0
		Base font size for title, axes labels, ticks, and annotations.

	Returns
	-------
	None
		Save one combined PNG figure as ``{outprefix}_cor_scatter.png``.
	"""
	if not cor_col:
		raise ValueError("cor_col is empty. Please provide at least one column pair.")
	if font_size <= 0:
		raise ValueError(f"font_size must be > 0, got {font_size}")

	# Ensure Chinese labels/titles render correctly in this plotting entrypoint.
	configure_chinese_font()

	in_path = Path(infile)
	if not in_path.exists():
		raise FileNotFoundError(f"Input file not found: {in_path}")

	sep = "\t" if in_path.suffix.lower() in {".tsv", ".txt"} else ","
	df = pd.read_csv(in_path, sep=sep)

	missing_cols = [
		col_name
		for x_col, y_col in cor_col
		for col_name in (x_col, y_col)
		if col_name not in df.columns
	]
	if missing_cols:
		missing_unique = list(dict.fromkeys(missing_cols))
		raise ValueError(f"Missing required columns: {missing_unique}")

	n_plots = len(cor_col)
	ncols = 2 if n_plots > 1 else 1
	nrows = int(np.ceil(n_plots / ncols))
	fig, axes = plt.subplots(
		nrows=nrows,
		ncols=ncols,
		figsize=(9 * ncols, 7 * nrows),
		squeeze=False,
	)
	flat_axes = axes.flatten()

	for idx, (x_col, y_col) in enumerate(cor_col):
		ax = flat_axes[idx]
		plot_df = df[[x_col, y_col]].copy()
		plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
		plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors="coerce")
		plot_df = plot_df.dropna(subset=[x_col, y_col])

		if plot_df.empty:
			ax.text(0.5, 0.5, "No numeric data", ha="center", va="center", fontsize=font_size)
			ax.set_title(f"{x_col} vs {y_col}", fontsize=font_size + 2)
			ax.set_xlabel(x_col, fontsize=font_size)
			ax.set_ylabel(y_col, fontsize=font_size)
			ax.tick_params(axis="both", labelsize=max(font_size - 1, 8))
			continue

		x_vals = plot_df[x_col].to_numpy(dtype=float)
		y_vals = plot_df[y_col].to_numpy(dtype=float)

		ax.scatter(x_vals, y_vals, alpha=0.7, s=24, edgecolors="none")
		ax.set_title(f"{x_col} vs {y_col}", fontsize=font_size + 2)
		ax.set_xlabel(x_col, fontsize=font_size)
		ax.set_ylabel(y_col, fontsize=font_size)
		ax.tick_params(axis="both", labelsize=max(font_size - 1, 8))
		ax.grid(alpha=0.25, linestyle="--")

		# Add a linear trend line when x has variability.
		if len(plot_df) >= 2 and np.unique(x_vals).size >= 2:
			slope, intercept = np.polyfit(x_vals, y_vals, 1)
			x_line = np.linspace(np.min(x_vals), np.max(x_vals), 100)
			y_line = slope * x_line + intercept
			ax.plot(x_line, y_line, color="#D62728", linewidth=1.6, linestyle="-")

		if len(plot_df) >= 2:
			corr = plot_df[x_col].corr(plot_df[y_col])
			ax.text(
				0.03,
				0.97,
				f"r = {corr:.3f}",
				transform=ax.transAxes,
				va="top",
				ha="left",
				fontsize=max(font_size - 1, 8),
				bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75},
			)

	for idx in range(n_plots, len(flat_axes)):
		flat_axes[idx].axis("off")

	fig.tight_layout()
	out_path = f"{outprefix}_cor_scatter.png"
	fig.savefig(out_path, dpi=300)
	plt.close(fig)
	print(f"Saved correlation plot: {out_path}")

if __name__ == "__main__":
    pass
    
