import gzip
import numpy as np
from typing import List, Dict, Union

def detect_phred_encoding(fastq_files: Union[str, List[str]]) -> Dict[str, str]:
    """
    Detect the Phred quality encoding (Phred+33 or Phred+64) of one or more FASTQ files.

    This function reads the quality lines of FASTQ files (supports plain text and .gz),
    samples a subset of lines for efficiency, and determines the encoding based on
    the ASCII character ranges.

    Parameters
    ----------
    fastq_files : str or List[str]
        Path to a FASTQ file or a list of FASTQ file paths. 
        Both uncompressed and .gz compressed files are supported.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping each input file to detected encoding:
        - 'phred33' for Illumina 1.8+ encoding (ASCII 33–73 for Q0–40)
        - 'phred64' for Illumina 1.3–1.7 encoding (ASCII 64–104 for Q0–40)
        - 'unknown' if detection is inconclusive

    Notes
    -----
    Detection rules:
    1. Sample the first N quality lines (default N=10000 or file length if shorter).
    2. Convert quality strings to ASCII codes using ord().
    3. If any ASCII code < 59, likely Phred+33 (ASCII 33–73 covers low-Q bases)
    4. If all ASCII codes >= 59, likely Phred+64 (ASCII 64–104)
    5. If ambiguous (mixed), return 'unknown'.
    Q value region[1~40]
    
    Performance:
    - Uses numpy arrays to vectorize ASCII code conversion and min/max checks.
    - Only reads quality lines (every 4th line in FASTQ).
    """
    if isinstance(fastq_files, str):
        fastq_files = [fastq_files]

    results = {}
    for fq_file in fastq_files:
        # Open file (support gzip)
        open_func = gzip.open if fq_file.endswith(".gz") else open
        with open_func(fq_file, 'rt', encoding='utf-8', errors='ignore') as f:
            line_count = 0
            ascii_codes = []
            for i, line in enumerate(f):
                # FASTQ quality lines: every 4th line (lines 4, 8, 12, ...)
                if (i + 1) % 4 == 0:
                    codes = np.fromiter((ord(c) for c in line.strip()), dtype=np.uint8)
                    ascii_codes.append(codes)
                    line_count += 1
                    if line_count >= 10000:  # sample first 10,000 quality lines for speed
                        break
            if not ascii_codes:
                results[fq_file] = 'unknown'
                continue

            # Concatenate all sampled quality lines
            all_codes = np.concatenate(ascii_codes)
            min_code = all_codes.min()
            max_code = all_codes.max()

            # Decision rules
            # Phred+33: ASCII 33–73; Phred+64: ASCII 64–104
            if min_code < 59:  # lower than 64 minus safety margin → phred33
                results[fq_file] = 'phred33'
            elif min_code >= 59 and max_code <= 104:
                results[fq_file] = 'phred64'
            else:
                results[fq_file] = 'unknown'
    print(results)
    return results

if __name__ == "__main__":
    fastq = "/data/pub/zhousha/20260411_RNAseq/data/fastq/B1_CoD1_Mouse_WT_rep1_1.fq.gz"
    detect_phred_encoding(
        fastq_files=fastq
    )
    pass