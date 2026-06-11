"""Plot fragment size distribution from histogram table."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse


def plot_fragment_size(input_file: str, output_file: str, highlight_size: int = 167) -> None:
    """Plot fragment size distribution.
    
    Args:
        input_file: Path to fragment size histogram table (from getFragmentSize.py)
        output_file: Path to output plot file
        highlight_size: Fragment size to highlight with vertical line (default: 167)
    """
    df = pd.read_table(input_file, index_col=0)
    
    fig, ax = plt.subplots(figsize=(25, 10), dpi=300)
    
    for col in df.columns:
        total = df[col].sum()
        if total > 0:
            df[col] = df[col] / total
        ax.plot(df.index, df[col], label=col)
    
    # Add vertical line at highlight_size
    ax.axvline(x=highlight_size, color='r', linestyle='--', alpha=0.7, label=str(highlight_size))
    
    # Add tick for highlight_size
    current_xticks = ax.get_xticks()
    new_xticks = np.append(current_xticks, highlight_size)
    ax.set_xticks(new_xticks)
    
    ax.set_xlabel('Fragment size (bp)')
    ax.set_ylabel('Frequency')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot fragment size distribution.")
    parser.add_argument('--input', type=str, required=True, help='Path to fragment size histogram table')
    parser.add_argument('--output', type=str, required=True, help='Path to output plot file')
    parser.add_argument('--highlight', type=int, default=167, help='Fragment size to highlight (default: 167)')
    args = parser.parse_args()
    plot_fragment_size(args.input, args.output, args.highlight)
