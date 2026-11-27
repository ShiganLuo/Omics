"""Plot fragment size distribution from histogram table."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from typing import List, Optional

PlotFormat = str

VALID_FORMATS = ("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg",
                 "jpeg", "pgf", "raw", "rgba")


def plot_fragment_size(input_file: str, out_base: str,
                       image_formats: Optional[List[PlotFormat]] = None,
                       highlight_size: int = 167) -> None:
    """Plot fragment size distribution.
    
    Args:
        input_file: Path to fragment size histogram table (from getFragmentSize.py)
        out_base: Output path without extension (``<base>.<fmt>`` per format)
        image_formats: Save formats. Default ``["png"]``.
        highlight_size: Fragment size to highlight with vertical line (default: 167)
    """
    if image_formats is None:
        image_formats = ["png"]
    bad = [f for f in image_formats if f not in VALID_FORMATS]
    if bad:
        raise ValueError(
            f"Unsupported image format(s): {', '.join(bad)}; "
            f"valid: {', '.join(VALID_FORMATS)}")
    df = pd.read_table(input_file, index_col=0)
    
    fig, ax = plt.subplots(figsize=(25, 10), dpi=300)
    
    for col in df.columns:
        total = df[col].sum()
        if total > 0:
            df[col] = df[col] / total
        ax.plot(df.index, df[col], label=col)
    
    ax.axvline(x=highlight_size, color='r', linestyle='--', alpha=0.7, label=str(highlight_size))
    
    current_xticks = ax.get_xticks()
    new_xticks = np.append(current_xticks, highlight_size)
    ax.set_xticks(new_xticks)
    
    ax.set_xlabel('Fragment size (bp)')
    ax.set_ylabel('Frequency')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    for fmt in image_formats:
        out_path = f"{out_base}.{fmt}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {out_path}")
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot fragment size distribution.")
    parser.add_argument('--input', type=str, required=True, help='Path to fragment size histogram table')
    parser.add_argument('--output', type=str, required=True,
                        help='Output path without extension (base name); '
                             'files are written as <output>.<format>')
    parser.add_argument('--highlight', type=int, default=167, help='Fragment size to highlight (default: 167)')
    parser.add_argument('-f', '--format', action='append', dest='formats', metavar='FMT',
                        help='Image format (png/pdf/svg/...). Repeatable. Default: png.')
    args = parser.parse_args()
    plot_fragment_size(args.input, args.output,
                       image_formats=args.formats if args.formats else ["png"],
                       highlight_size=args.highlight)
