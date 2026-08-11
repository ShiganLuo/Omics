from typing import Literal
import re
import logging
logger = logging.getLogger(__name__)
from collections import Counter
import re
from typing import Literal


def detect_delimiter(
    file_path: str,
    sample_lines: int = 20,
) -> Literal[",", "\t", ";", "|", "whitespace"]:
    """
    Detect delimiter of metadata / table file.

    Strategy:
    1. Test common delimiters.
    2. Calculate dominant column count ratio.
    3. Select delimiter with highest consistency.
    4. Fallback to whitespace.
    """

    candidates = [",", "\t", ";", "|"]

    with open(file_path, "r", encoding="utf-8") as f:
        lines = [
            line.rstrip("\n")
            for _, line in zip(range(sample_lines), f)
            if line.strip()
        ]

    if not lines:
        raise ValueError("Metadata file is empty.")

    best_delim = None
    best_score = 0
    best_columns = 0


    for delim in candidates:

        counts = [
            len(line.split(delim))
            for line in lines
        ]

        counter = Counter(counts)

        # 最常见列数
        dominant_cols, freq = counter.most_common(1)[0]

        ratio = freq / len(counts)

        # 至少产生2列
        if dominant_cols <= 1:
            continue

        # 综合评分
        score = ratio * dominant_cols

        if score > best_score:
            best_score = score
            best_delim = delim
            best_columns = dominant_cols


    if best_delim:
        logger.info(
            f"Detected delimiter '{best_delim}' "
            f"({best_columns} columns, consistency={best_score:.2f})"
        )
        return best_delim


    # whitespace
    counts = [
        len(re.split(r"\s+", line.strip()))
        for line in lines
    ]

    counter = Counter(counts)

    cols, freq = counter.most_common(1)[0]

    if cols > 1 and freq / len(counts) >= 0.8:
        logger.info(
            f"Detected whitespace delimiter "
            f"({cols} columns)"
        )
        return "whitespace"


    raise ValueError(
        "Could not determine file delimiter."
    )