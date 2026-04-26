import gzip
import logging
from collections import Counter
from pathlib import Path
from typing import Union, Tuple, Dict

# Setup logger
logger = logging.getLogger("fastq_diversity")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def read_fastq(file_path: Union[str, Path]) -> Tuple[str, str, str, str]:
    """
    Generator function to read FASTQ file, supports gzip-compressed files.

    Parameters
    ----------
    file_path : str or Path
        Path to the FASTQ file.

    Yields
    ------
    Tuple[str, str, str, str]
        Four lines of a FASTQ record: header, sequence, plus line, quality.
    """
    file_path = Path(file_path)
    open_func = gzip.open if file_path.suffix in {".gz", ".gzip"} else open

    logger.info(f"Reading FASTQ file: {file_path}")
    with open_func(file_path, "rt") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            yield header, seq, plus, qual


def count_prefix_diversity(
    file_path: Union[str, Path],
    start: int = 0,
    end: int = 12
) -> Dict[str, int]:
    """
    Count nucleotide diversity in the prefix of FASTQ sequences.

    Parameters
    ----------
    file_path : str or Path
        Path to the FASTQ file.
    start : int, default=0
        Start position (0-based) of the prefix.
    end : int, default=12
        End position of the prefix (exclusive).

    Returns
    -------
    Dict[str, int]
        Counts of each unique prefix.
    """
    prefix_counter: Counter = Counter()
    total_reads = 0

    for _, seq, _, _ in read_fastq(file_path):
        prefix = seq[start:end]
        prefix_counter[prefix] += 1
        total_reads += 1
        if total_reads % 1_000_000 == 0:
            logger.info(f"Processed {total_reads} reads...")

    logger.info(f"Finished processing {total_reads} reads")
    return dict(prefix_counter)


def report_top_prefixes(prefix_counts: Dict[str, int], top_n: int = 20) -> None:
    """
    Print the top N most frequent prefixes.

    Parameters
    ----------
    prefix_counts : Dict[str, int]
        Dictionary of prefix counts.
    top_n : int, default=20
        Number of top prefixes to display.
    """
    logger.info(f"Reporting top {top_n} prefixes:")
    for prefix, count in sorted(prefix_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]:
        print(f"{prefix}\t{count}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="FASTQ prefix diversity counter")
    parser.add_argument("fastq", help="Path to FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("--start", type=int, default=0, help="Start position of prefix")
    parser.add_argument("--end", type=int, default=12, help="End position of prefix (exclusive)")
    parser.add_argument("--top", type=int, default=20, help="Number of top prefixes to display")
    args = parser.parse_args()

    counts = count_prefix_diversity(args.fastq, start=args.start, end=args.end)
    report_top_prefixes(counts, top_n=args.top)