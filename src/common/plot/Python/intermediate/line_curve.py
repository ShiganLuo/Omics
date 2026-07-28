from pathlib import Path
from typing import Iterable, Literal, Tuple, cast
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline, PchipInterpolator, UnivariateSpline

SmoothMethod = Literal["none", "linear", "cubic", "pchip", "spline"]
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
    """Sort by x in ascending order and merge duplicate x values by mean y."""
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
    method : {'none', 'linear', 'cubic', 'pchip', 'spline'}, default='pchip'
        Curve generation method:
        - 'none'  : no smoothing, return sorted raw points only;
        - 'linear': linear interpolation;
        - 'cubic' : cubic spline (falls back to linear if < 4 points);
        - 'pchip' : shape-preserving piecewise cubic interpolation;
        - 'spline': smoothing spline controlled by `spline_s`.
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

    if method not in {"none", "linear", "cubic", "pchip", "spline"}:
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
    else:  # method == 'spline'
        degree = min(3, len(x) - 1)
        y_new = cast(
            FloatArray,
            np.asarray(UnivariateSpline(x, y, s=spline_s, k=degree)(x_new), dtype=np.float64),
        )

    return x_new, y_new


def plot_depth_sensitivity_curve(
    depth: Iterable[float],
    sensitivity: Iterable[float],
    method: SmoothMethod = "pchip",
    num_points: int = 300,
    spline_s: float = 0.0,
    x_label: str = "Sequencing Depth",
    y_label: str = "Sensitivity",
    title: str = "Depth-Sensitivity Curve",
    output_path: str | Path | None = None,
    dpi: int = 180,
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
    output_path : str | Path | None, default=None
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


    x_raw = _to_float_array(depth, "depth")
    y_raw = _to_float_array(sensitivity, "sensitivity")

    x_smooth, y_smooth = smooth_depth_sensitivity_curve(
        x_raw,
        y_raw,
        method=method,
        num_points=num_points,
        spline_s=spline_s,
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x_raw, y_raw, color="#1f77b4", alpha=0.85, label="Raw points")

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
    ax.legend(loc="best")
    fig.tight_layout()

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


def demo() -> None:
    """Minimal usage example using direct in-memory points."""
    depth = [100, 200, 300, 500, 800, 1000]
    sensitivity = [0.38, 0.56, 0.68, 0.82, 0.91, 0.95]

    plot_depth_sensitivity_curve(
        depth,
        sensitivity,
        method="pchip",
        x_label="Sequencing Depth",
        y_label="Detection Sensitivity",
        output_path="depth_sensitivity_curve.png",
        show=False,
    )


if __name__ == "__main__":
    demo()
