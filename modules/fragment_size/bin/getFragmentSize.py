"""Extract fragment size distribution from samtools stats output."""
import re
import pandas as pd
import argparse
from typing import List


def fragment_hist(inlist: List[str], outfile: str) -> None:
    """Parse samtools stats files and create fragment size histogram table.
    
    Args:
        inlist: List of samtools stats output file paths
        outfile: Output file path for the histogram table
    """
    records = []
    for path in inlist:
        sample_id = re.sub(r"\.txt$", "", path.split("/")[-1])
        with open(path) as f:
            for line in f:
                if not line.startswith("IS"):
                    continue
                # IS, insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs
                fields = line.strip().split("\t")
                IS, N = fields[1], fields[3]
                records.append((sample_id, int(IS), int(N)))
    
    fs_hist_table = pd.DataFrame.from_records(records)
    fs_hist_table.columns = ["sample_id", "insertion-size", "count"]
    fs_hist = fs_hist_table.pivot(index="insertion-size", columns="sample_id", values="count")
    fs_hist.to_csv(outfile, sep="\t")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract fragment size distribution from samtools stats.")
    parser.add_argument('--input', type=str, nargs='+', required=True, help='List of samtools stats file paths')
    parser.add_argument('--out', type=str, required=True, help='Path to output file')
    args = parser.parse_args()
    fragment_hist(args.input, args.out)
