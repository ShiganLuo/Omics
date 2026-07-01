import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from common.LogUtil import setup_logger
from common.CmdUtil import _run_cmd
from pathlib import Path
import logging
import pandas as pd
from scipy.stats import fisher_exact
from typing import List, Literal, Dict, defaultdict, Optional
import logging
from pathlib import Path
import tempfile
import shutil

logger = setup_logger("RepeatMaskerAnalysis", level=logging.INFO)


def run_te_annotation_pipeline(
    input_fasta: str, 
    output_dir: str, 
    species: str = "mus musculus",
    threads: int = 8
) -> Optional[Path]:
    """
    Execute the RepeatMasker pipeline for transposable element annotation.

    Runs RepeatMasker on the input FASTA, moves results to the output
    directory, and cleans up temporary working folders.

    Parameters
    ----------
    input_fasta : str
        Path to the input FASTA file.
    output_dir : str
        Target directory for the generated annotation results.
    species : str, optional
        Genomic repeat library species identifier. Default is "mus musculus".
    threads : int, optional
        Number of CPU cores for parallel execution. Default is 8.

    Returns
    -------
    Optional[Path]
        Path to the primary .out result file, or None if unsuccessful.
    """
    fasta_path = Path(input_fasta)
    final_out_path = Path(output_dir)
    final_out_path.mkdir(parents=True, exist_ok=True)

    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        logger.error(f"Input FASTA invalid: {input_fasta}")
        return None

    with tempfile.TemporaryDirectory(dir=final_out_path) as tmp_dir:
        tmp_path = Path(tmp_dir)
        
        rm_cmd = [
            "RepeatMasker",
            "-pa", str(threads),
            "-species", species,
            "-dir", str(tmp_path),
            "-gff",
            str(fasta_path)
        ]

        try:
            stdout = _run_cmd(rm_cmd)
            logger.info(f"RepeatMasker completed successfully:\n{stdout}")
            
            for item in tmp_path.iterdir():
                if item.is_file():
                    shutil.move(str(item), str(final_out_path / item.name))
                elif item.is_dir() and not item.name.startswith("RM_"):
                    shutil.move(str(item), str(final_out_path / item.name))

        except Exception as e:
            logger.error(f"RepeatMasker failed during execution: {e}")
            return None

    result_out = final_out_path / f"{fasta_path.name}.out"
    
    summary_file = final_out_path / f"{fasta_path.name}.tbl"
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            logger.info(f"\n{'='*20} RepeatMasker Summary {'='*20}\n{f.read()}")

    return result_out if result_out.exists() else None


class RepeatMaskerOutCompare:
    """
    Compare repeat element composition between foreground and background.
    Supports analysis at class, family, and subfamily (specific name) levels.
    """

    def __init__(self, bg_out: str, fg_out: str):
        """
        Initialize the comparator with foreground and background .out files.

        Parameters
        ----------
        bg_out : str
            Path to the background RepeatMasker .out file.
        fg_out : str
            Path to the foreground RepeatMasker .out file.
        """
        self.bg_out = bg_out
        self.fg_out = fg_out

    def parse_repeatmasker_out(
        self, 
        path: str, 
        min_score: int = 225, 
        max_div: float = 25.0,
        min_len: int = 10
    ) -> List[Dict]:
        """
        Parse a RepeatMasker .out file extracting subfamily, class, and family.

        Parameters
        ----------
        path : str
            Path to the RepeatMasker .out file.
        min_score : int, optional
            Minimum Smith-Waterman score to retain a hit. Default is 225.
        max_div : float, optional
            Maximum percent divergence allowed. Default is 25.0.
        min_len : int, optional
            Minimum fragment length in bp. Default is 10.

        Returns
        -------
        List[Dict]
            List of repeat records with keys: subfamily, class, family,
            length, div.
        """
        repeats = []

        with open(path) as f:
            for line in f:
                if line.startswith(("SW", "score", "#")) or not line.strip():
                    continue

                fields = line.strip().split()
                if len(fields) < 11:
                    continue

                try:
                    # Column Indices: 0: Score, 1: Div, 5: Start, 6: End, 9: Subfamily, 10: Class/Family
                    sw_score = int(fields[0])
                    perc_div = float(fields[1])
                    q_start = int(fields[5])
                    q_end = int(fields[6])
                    fragment_len = abs(q_end - q_start) + 1
                    
                    if sw_score < min_score or perc_div > max_div or fragment_len < min_len:
                        continue
                    
                    subfamily = fields[9]  # Specific name, e.g., L1Md_T
                    class_family = fields[10] # e.g., LINE/L1
                    
                    if "/" in class_family:
                        repeat_class, repeat_family = class_family.split("/", 1)
                    else:
                        repeat_class = class_family
                        repeat_family = "Unknown"

                    repeats.append({
                        "subfamily": subfamily,
                        "class": repeat_class,
                        "family": repeat_family,
                        "length": fragment_len,
                        "div": perc_div
                    })
                    
                except (ValueError, IndexError):
                    continue

        return repeats

    def summarize_lengths(
        self,
        repeats: List[Dict],
        level: Literal["class", "family", "subfamily"] = "class"
    ) -> Dict[str, int]:
        """
        Aggregate total base pairs per repeat element at the specified level.

        Parameters
        ----------
        repeats : List[Dict]
            List of repeat records as returned by ``parse_repeatmasker_out``.
        level : {"class", "family", "subfamily"}, optional
            Taxonomic level at which to aggregate. Default is "class".

        Returns
        -------
        Dict[str, int]
            Mapping from repeat element name to total base pairs.
        """
        length_map = defaultdict(int)
        for r in repeats:
            key = r[level]
            length_map[key] += r["length"]
        return dict(length_map)

    def enrichment_test(
        self,
        level: Literal["class", "family", "subfamily"] = "subfamily",
        min_score: int = 225,
        max_div: float = 25.0,
        min_len: int = 10
    ) -> pd.DataFrame:
        """
        Perform Fisher's exact test enrichment analysis on repeat elements.

        Compares foreground vs. background base-pair proportions for each
        repeat at the given taxonomic level.

        Parameters
        ----------
        level : {"class", "family", "subfamily"}, optional
            Taxonomic level for comparison. Default is "subfamily".
        min_score : int, optional
            Minimum Smith-Waterman score filter. Default is 225.
        max_div : float, optional
            Maximum percent divergence filter. Default is 25.0.
        min_len : int, optional
            Minimum fragment length in bp. Default is 10.

        Returns
        -------
        pd.DataFrame
            Enrichment results with columns: repeat, fg_bp, bg_bp,
            fg_ratio, bg_ratio, odds_ratio, p_value, fdr.
            Sorted by p-value ascending.
        """
        fg_raw = self.parse_repeatmasker_out(
            self.fg_out, min_score=min_score, max_div=max_div, min_len=min_len
        )
        bg_raw = self.parse_repeatmasker_out(
            self.bg_out, min_score=min_score, max_div=max_div, min_len=min_len
        )

        fg_lengths = self.summarize_lengths(fg_raw, level=level)
        bg_lengths = self.summarize_lengths(bg_raw, level=level)

        fg_total_bp = sum(fg_lengths.values())
        bg_total_bp = sum(bg_lengths.values())

        records = []
        all_keys = set(fg_lengths) | set(bg_lengths)

        for k in all_keys:
            a = fg_lengths.get(k, 0)
            b = fg_total_bp - a
            c = bg_lengths.get(k, 0)
            d = bg_total_bp - c

            if a + c == 0: continue

            odds_ratio, p_value = fisher_exact([[a, b], [c, d]])

            records.append({
                "repeat": k,
                "fg_bp": a,
                "bg_bp": c,
                "fg_ratio": a / fg_total_bp if fg_total_bp else 0,
                "bg_ratio": c / bg_total_bp if bg_total_bp else 0,
                "odds_ratio": odds_ratio,
                "p_value": p_value
            })

        if not records: return pd.DataFrame()

        df = pd.DataFrame(records)
        df["fdr"] = df["p_value"].rank(method="min") / len(df)
        return df.sort_values("p_value")


if __name__ == "__main__":
    # 执行分析
    params = {
        "bg_out":"",
        "fg_out": ""
    }
    repeatmakerMask = RepeatMaskerOutCompare(**params)