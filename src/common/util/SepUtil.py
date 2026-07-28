from typing import Literal
import re
import logging
logger = logging.getLogger(__name__)
def detect_delimiter(
    file_path: str,
    sample_lines: int = 20,
) -> Literal[",", "\t", ";", "|", "whitespace"]:
    """
    Detect delimiter of a metadata / table file.

    Strategy:
    1. Test common delimiters
    2. Choose delimiter producing the most consistent column counts
    3. Fallback to whitespace detection
    """

    candidates = [",", "\t", ";", "|"]

    with open(file_path, "r", encoding="utf-8") as f:
        lines = [line for _, line in zip(range(sample_lines), f)] # not strip() to preserve potential leading/trailing delimiters

    lines = [l for l in lines if l]

    if not lines:
        raise ValueError("Metadata file is empty.")

    best_delim = None
    best_score = -1

    for delim in candidates:
        col_counts = [len(line.split(delim)) for line in lines]
        # consistency score
        unique_counts = set(col_counts)

        if len(unique_counts) == 1:
            score = col_counts[0]
        else:
            score = 0

        if score > best_score:
            best_score = score
            best_delim = delim

    if best_delim and best_score > 1:
        logger.info(f"Detected delimiter: '{best_delim}' with {best_score} columns for {file_path}")
        return best_delim

    # whitespace detection
    col_counts = [
        len(re.split(r"\s+", line))
        for line in lines
    ]

    if len(set(col_counts)) == 1 and col_counts[0] > 1:
        logger.info(f"Detected delimiter: 'whitespace' with {col_counts[0]} columns for {file_path}")
        return "whitespace"
    raise ValueError("Could not determine file delimiter.")
