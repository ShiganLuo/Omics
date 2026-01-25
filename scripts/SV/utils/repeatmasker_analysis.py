from pathlib import Path
import logging
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
from scipy.stats import fisher_exact
from collections import Counter
from typing import List, Tuple, Literal, Dict, defaultdict
import subprocess
import logging
from pathlib import Path
import tempfile
import shutil

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)

def run_cmd(cmd:list):
        """
        执行外部命令，返回 stdout
        - 命令不存在：给出清晰提示
        - 命令执行失败：打印 stdout / stderr
        """
        cmd_str = " ".join(cmd)
        cmd_bin = cmd[0]

        logger.info(f"Running: {cmd_str}")

        # 1️⃣ 预检查：命令是否存在（比 FileNotFoundError 更友好）
        if shutil.which(cmd_bin) is None:
            logger.error(f"Command not found: '{cmd_bin}'")
            logger.error("Please make sure it is installed and in $PATH")
            raise RuntimeError(f"Command not found: {cmd_bin}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )

            if result.stdout:
                logger.info(f"Command Output:\n{result.stdout}")

            return result.stdout

        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with return code {e.returncode}")
            logger.error(f"STDOUT:\n{e.stdout or '[empty]'}")
            logger.error(f"STDERR:\n{e.stderr or '[empty]'}")
            raise RuntimeError(
                f"Command execution failed: {cmd_str}"
            ) from e


def run_te_annotation_pipeline(
    input_fasta: str, 
    output_dir: str, 
    species: str = "mus musculus",
    threads: int = 8
) -> Path:
    """
    Executes the RepeatMasker pipeline, ensuring all output files are preserved 
    in the destination directory while cleaning up temporary working folders.

    Args:
        input_fasta: Path to the input FASTA file.
        output_dir: Target directory for the generated annotation results.
        species: Genomic repeat library species identifier.
        threads: Number of CPU cores for parallel execution.

    Returns:
        Path: The path to the primary .out result file, or None if unsuccessful.
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
            run_cmd(rm_cmd)
            
            for item in tmp_path.iterdir():
                if item.is_file():
                    shutil.move(str(item), str(final_out_path / item.name))
                elif item.is_dir() and not item.name.startswith("RM_"):
                    shutil.move(str(item), str(final_out_path / item.name))

        except subprocess.CalledProcessError:
            logger.error("RepeatMasker failed during execution.")
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
        Parse .out file including subfamily (repeat name), class, and family.
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
        Sum up total base pairs for the chosen level.
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
        Perform enrichment analysis with BP lengths. 
        Set max_div low (e.g., 5.0) to find active L1 subfamilies.
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