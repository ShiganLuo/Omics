import sys
import os
import pandas as pd
import gzip
import logging
import numpy as np
from pathlib import Path
from typing import Optional
import json
import re
try:
    from .common import setup_logger
except ImportError:
    from common import setup_logger
logger = setup_logger("SV_TYPE", level=logging.INFO)



def open_vcf(path: str):
    """
    Open a VCF or VCF.GZ file for reading.

    Parameters
    ----------
    path : str
        Path to a VCF file. Gzip-compressed (``.vcf.gz``) files are handled
        automatically.

    Returns
    -------
    file object
        An open text-mode file handle.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_pbsv_vcf(vcf_file: str) -> pd.DataFrame:
    """
    Parse a pbsv SV VCF file into a DataFrame.

    Extracts chrom, start, end, svtype, and svlen from each variant record.

    Parameters
    ----------
    vcf_file : str
        Path to the pbsv VCF file (plain or gzip-compressed).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``'chrom'``, ``'start'``, ``'end'``,
        ``'svtype'``, and ``'svlen'``.
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
    Generate plot-ready data by log-space binning of SV lengths.

    Computes per-SV-type proportional distributions over log-spaced bins
    suitable for smooth curve plotting.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing ``'svtype'`` and ``'svlen'`` columns.
    bins : int, optional
        Number of logarithmic bins. Higher values give finer resolution.
        Default is ``50``.

    Returns
    -------
    list of dict
        Each dict contains keys ``'x'`` (bin midpoints), ``'y'`` (proportions),
        ``'label'`` (SV type), and ``'color'`` (hex color string).
    """
    color_map = {
        "DEL": "#E41A1C", "INS": "#377EB8", 
        "DUP": "#4DAF4A", "TRA": "#984EA3", "INV": "#FF7F00"
    }
    default_color = "#999999"

    df_clean = df[df['svlen'] > 0].copy()
    
    min_log = np.log10(df_clean['svlen'].min())
    max_log = np.log10(df_clean['svlen'].max())
    bin_edges = np.logspace(min_log, max_log, bins + 1)
    
    bin_mids = np.sqrt(bin_edges[:-1] * bin_edges[1:]) 

    total_count = len(df_clean)

    plot_data = []
    
    for sv_type, group in df_clean.groupby('svtype'):
        counts, _ = np.histogram(group['svlen'], bins=bin_edges)
        
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
    Run SV stratification analysis and generate summary reports.

    Parses a pbsv VCF, bins SVs by size and type, and exports a
    cross-tabulation matrix, chromosome distribution, and JSON summary.

    Parameters
    ----------
    vcf_file : str or Path
        Path to the input VCF file.
    output_dir : str or Path
        Directory where output reports will be saved.
    run_cmd_func : callable, optional
        Custom function for executing shell commands. If ``None``, the
        default ``subprocess`` runner is used. Currently unused.

    Returns
    -------
    tuple of (pd.DataFrame, dict)
        A cross-tabulation matrix of SV type vs. size group, and a summary
        dictionary with total count, type/size breakdowns, and
        heavy-rearrangement count.
    """
    vcf_path = Path(vcf_file)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Generating stratification reports for: {vcf_path.name}")

    df = parse_pbsv_vcf(vcf_file)

    bins = [0, 500, 2000, 10000, 50000, np.inf]
    labels = ['50-500bp', '500bp-2kb', '2kb-10kb', '10kb-50kb', '>50kb']
    df['size_group'] = pd.cut(df['svlen'], bins=bins, labels=labels)

    matrix = pd.crosstab(df['svtype'], df['size_group'])
    matrix.to_csv(out_path / "stratification_matrix.csv")

    chrom_dist = df['chrom'].value_counts().reset_index()
    chrom_dist.columns = ['Chromosome', 'Count']
    chrom_dist.to_csv(out_path / "chromosome_distribution.csv", index=False)

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
    output_fasta_str: str,
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
    output_fasta = Path(output_fasta_str)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(__name__)
    ins_count = 0

    svlen_re = re.compile(r"SVLEN=([^;]+)")
    svtype_re = re.compile(r"SVTYPE=([^;]+)")
    end_re = re.compile(r"END=([^;]+)")

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

                st_match = svtype_re.search(info)
                svtype = st_match.group(1) if st_match else None

                length = None

                sl_match = svlen_re.search(info)
                if sl_match:
                    try:
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

                if svtype != "INS":
                    continue
                if length is None or not (min_len <= length <= max_len):
                    continue

                if alt.startswith("<"):
                    continue

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