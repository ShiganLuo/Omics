import pandas as pd
import gzip
import logging
import numpy as np
import subprocess
from pathlib import Path
from typing import Optional
import json
import re
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


def open_vcf(path: str):
    """
    Open VCF or VCF.GZ file
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_pbsv_vcf(vcf_file: str) -> pd.DataFrame:
    """
    Parse pbsv SV VCF into DataFrame
    """
    records = []

    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            chrom, pos, _, _, _, _, _, info = fields[:8]
            pos = int(pos)

            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v

            svtype = info_dict.get("SVTYPE", "NA")
            end = int(info_dict.get("END", pos))

            svlen_raw = info_dict.get("SVLEN")
            if svlen_raw is not None:
                try:
                    svlen = abs(int(svlen_raw))
                except ValueError:
                    svlen = abs(end - pos)
            else:
                svlen = abs(end - pos)

            records.append(
                {
                    "chrom": chrom,
                    "start": pos,
                    "end": end,
                    "svtype": svtype,
                    "svlen": svlen,
                }
            )

    return pd.DataFrame(records)

def generate_plot_input(df, bins=50):
    """
    对 svlen 进行对数空间分箱，计算频率占比
    :param df: 包含 svtype 和 svlen 的 DataFrame
    :param bins: 分箱数量，越多越精细（建议 30-100）
    """
    color_map = {
        "DEL": "#E41A1C", "INS": "#377EB8", 
        "DUP": "#4DAF4A", "TRA": "#984EA3", "INV": "#FF7F00"
    }
    default_color = "#999999"

    # 1. 预处理：移除长度 <= 0 的数据
    df_clean = df[df['svlen'] > 0].copy()
    
    # 2. 在对数空间定义箱边缘 (Log-spaced bins)
    min_log = np.log10(df_clean['svlen'].min())
    max_log = np.log10(df_clean['svlen'].max())
    bin_edges = np.logspace(min_log, max_log, bins + 1)
    
    # 计算箱的中点 (用于 X 轴绘图)
    bin_mids = np.sqrt(bin_edges[:-1] * bin_edges[1:]) 

    # 计算总样本量，用于计算比例
    total_count = len(df_clean)

    plot_data = []
    
    # 3. 按类型分组计算
    for sv_type, group in df_clean.groupby('svtype'):
        # 使用 np.histogram 计算每个箱子落入的频数
        counts, _ = np.histogram(group['svlen'], bins=bin_edges)
        
        # 将频数转换为占总体的比例 (0-1 之间)
        proportions = counts / total_count
        
        plot_data.append({
            "x": bin_mids.tolist(),
            "y": proportions.tolist(),
            "label": sv_type,
            "color": color_map.get(sv_type, default_color)
        })
        
    return plot_data


def run_sv_stratification(vcf_file, output_dir, run_cmd_func=None):
    """
    独立的分层统计函数
    :param vcf_file: 输入的 VCF 文件路径 (str 或 Path)
    :param output_dir: 统计结果存放目录
    :param run_cmd_func: 可选。用于执行 shell 命令的函数。如果为 None，则使用默认的 subprocess。
    """
    vcf_path = Path(vcf_file)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Generating stratification reports for: {vcf_path.name}")

    df = parse_pbsv_vcf(vcf_file)

    # 3. 长度分层 (Binning)
    bins = [0, 500, 2000, 10000, 50000, np.inf]
    labels = ['50-500bp', '500bp-2kb', '2kb-10kb', '10kb-50kb', '>50kb']
    df['size_group'] = pd.cut(df['svlen'], bins=bins, labels=labels)

    # A. 交叉矩阵 (Type vs Size)
    matrix = pd.crosstab(df['svtype'], df['size_group'])
    matrix.to_csv(out_path / "stratification_matrix.csv")

    # B. 染色体分布
    chrom_dist = df['chrom'].value_counts().reset_index()
    chrom_dist.columns = ['Chromosome', 'Count']
    chrom_dist.to_csv(out_path / "chromosome_distribution.csv", index=False)

    # C. 简报 JSON
    summary = {
        "vcf_source": str(vcf_path.absolute()),
        "total_count": len(df),
        "type_summary": df['svtype'].value_counts().to_dict(),
        "size_summary": df['size_group'].value_counts().to_dict(),
        "heavy_rearrangements": int(df['svtype'].isin(['BND', 'TRA', 'INV']).sum())
    }
    
    with open(out_path / "summary_report.json", "w") as j:
        json.dump(summary, j, indent=4)

    logger.info(f"Reports exported successfully to {out_path}")
    return matrix, summary



def extract_te_candidate_ins(
    vcf_file: str,
    output_fasta: str,
    min_len: int = 2000,
    max_len: int = 10000
) -> Optional[Path]:
    """
    Extract candidate transposable-element–associated insertion (INS) sequences
    from a VCF file and write them to a FASTA file.

    This function parses structural variants (SVs) from a VCF file and extracts
    insertion (INS) events whose estimated length falls within a specified range.
    The inserted sequences are taken from the ALT field and written to a FASTA file,
    which can be used for downstream analyses such as RepeatMasker annotation.

    Length determination strategy:
        1. Prefer the value of the INFO field "SVLEN" (absolute value of the first entry
           if multiple values are provided).
        2. If SVLEN is unavailable, estimate length as |END - POS|.

    Filtering criteria:
        - SVTYPE must be "INS"
        - Insertion length must satisfy: min_len <= length <= max_len
        - ALT must contain an explicit sequence (symbolic alleles such as "<INS>"
          are excluded)

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF file. Both plain-text ".vcf" and gzip-compressed
        ".vcf.gz" formats are supported.
    output_fasta : str
        Path to the output FASTA file for extracted insertion sequences.
        Parent directories will be created automatically if they do not exist.
    min_len : int, optional
        Minimum insertion length to retain (default: 2000 bp).
    max_len : int, optional
        Maximum insertion length to retain (default: 10000 bp).

    Returns
    -------
    pathlib.Path or None
        Path to the output FASTA file if extraction succeeds, otherwise None.

    Notes
    -----
    - This function assumes a standard VCF format with at least 8 columns.
    - For symbolic INS records (e.g. ALT="<INS>"), no actual inserted sequence
      is available and such records are skipped.
    - The FASTA header format is:
        >{ID}_{CHROM}_{POS}_len{LENGTH}
    """
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(__name__)
    ins_count = 0

    # Regular expressions for parsing INFO fields
    svlen_re = re.compile(r"SVLEN=([^;]+)")
    svtype_re = re.compile(r"SVTYPE=([^;]+)")
    end_re = re.compile(r"END=([^;]+)")

    # Choose open function based on file suffix
    open_func = gzip.open if str(vcf_file).endswith(".gz") else open

    try:
        with open_func(vcf_file, "rt") as f_in, open(output_fasta, "w") as f_out:
            for line in f_in:
                if line.startswith("#"):
                    continue

                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue

                chrom = cols[0]
                pos = int(cols[1])
                svid = cols[2]
                alt = cols[4]
                info = cols[7]

                # 1. Parse SVTYPE
                st_match = svtype_re.search(info)
                svtype = st_match.group(1) if st_match else None

                # 2. Determine insertion length
                length = None

                sl_match = svlen_re.search(info)
                if sl_match:
                    try:
                        # SVLEN may be negative or comma-separated
                        length = abs(int(float(sl_match.group(1).split(",")[0])))
                    except ValueError:
                        pass

                if length is None:
                    e_match = end_re.search(info)
                    if e_match:
                        try:
                            length = abs(int(e_match.group(1)) - pos)
                        except ValueError:
                            pass

                # 3. Apply filters
                if svtype != "INS":
                    continue
                if length is None or not (min_len <= length <= max_len):
                    continue

                # Skip symbolic ALT alleles
                if alt.startswith("<"):
                    continue

                # Write FASTA entry
                header = f">{svid}_{chrom}_{pos}_len{length}"
                f_out.write(f"{header}\n{alt}\n")
                ins_count += 1

        if ins_count > 0:
            logger.info(f"Extracted {ins_count} INS sequences to: {output_fasta}")
        else:
            logger.warning("No INS sequences matched the filtering criteria.")

    except Exception as e:
        logger.error(f"Failed to extract INS sequences from VCF: {e}")
        return None

    return output_fasta