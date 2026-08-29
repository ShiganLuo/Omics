"""Quantify transposable element expression from single-cell BAM using scTE.

Wraps scTE CLI: reads BAM with cell barcode (CB) and UMI (UB) tags,
produces TE count matrix per cell. Output converted to h5ad for Scanpy.

Usage:
    python scTE_quantify.py \
        --input sample.bam \
        --output sample_TE.h5ad \
        --index /path/to/hg38.exclusive.idx \
        --cb-tag CB --umi-tag UB \
        --threads 20
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

# Add project root to sys.path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import anndata as ad
import pandas as pd
import scipy.sparse as sp

from src.common.util.LogUtil import setup_logger


def _count_reads(bam, threads=1):
    """Count reads in a BAM file using samtools."""
    try:
        result = subprocess.run(
            ["samtools", "view", "-c", "-@", str(threads), bam],
            capture_output=True, text=True, check=True
        )
        return int(result.stdout.strip())
    except Exception:
        return -1


def run_scTE(bam, outdir, index, cb_tag, umi_tag, threads, scte_bin, logger, cwd=None):
    """Run scTE quantify on a BAM file.

    Pre-filters BAM to only include reads with both CB and UMI tags.
    Cell Ranger BAMs contain reads without CB/UB (failed barcode correction)
    and scTE's check requires >= 100/100 of the first reads to have tags.

    Filter: samtools -e '[CB] && [UB]' — tag-level precision, single pass.

    Note: scTE reads the BAM twice (check_cb_umi + bam2bed), so a FIFO
    cannot be used. The filtered BAM is written to a temp file on disk.
    """
    filtered_bam = bam
    tmp_bam = None

    tags_to_filter = []
    if cb_tag and cb_tag != "False":
        tags_to_filter.append(cb_tag)
    if umi_tag and umi_tag != "False":
        tags_to_filter.append(umi_tag)

    if tags_to_filter:
        tmp_bam = os.path.join(outdir, f"_scTE_filtered_{os.getpid()}.bam")

        tag_expr = " && ".join(f"[{t}]" for t in tags_to_filter)
        filter_cmd = [
            "samtools", "view", "-e", tag_expr,
            "-@", str(threads),
            "-b", "-o", tmp_bam, bam,
        ]
        logger.info(f"BAM filter: {bam}")
        logger.info(f"  Filter expression: {tag_expr}")
        logger.info(f"  Command: {' '.join(filter_cmd)}")
        subprocess.check_call(filter_cmd)

        n_after = _count_reads(tmp_bam, threads)
        logger.info(f"  Reads after filtering: {n_after}")

        filtered_bam = tmp_bam

    cmd = [
        scte_bin, "-i", filtered_bam, "-o", outdir,
        "-x", index,
        "-p", str(threads),
        "-CB", cb_tag, "-UMI", umi_tag,
    ]
    logger.info(f"Running: {' '.join(cmd)}")
    try:
        subprocess.check_call(cmd, cwd=cwd)
    finally:
        if tmp_bam and os.path.exists(tmp_bam):
            os.remove(tmp_bam)


def scTE_csv_to_h5ad(csv_path, sample_id=""):
    """Convert scTE CSV output to h5ad AnnData object.

    scTE outputs a CSV with genes as rows and cells as columns.
    Transpose to standard AnnData format (obs=cells, var=genes).
    """
    data = pd.read_csv(csv_path, index_col=0, header=0)
    data.index = data.index.astype(str).str.replace("/", "_", regex=False)
    if data.index.name and "/" in str(data.index.name):
        data.index.name = str(data.index.name).replace("/", "_")
    # scTE: rows=genes, cols=cells → transpose to rows=cells, cols=genes
    data = data.T
    X = sp.csr_matrix(data.values.astype("float32"))
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame(index=data.index),
        var=pd.DataFrame(index=data.columns),
    )
    if sample_id:
        adata.obs["sample_id"] = sample_id
    return adata


def main():
    parser = argparse.ArgumentParser(
        description="scTE TE quantification wrapper"
    )
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output h5ad path")
    parser.add_argument("--index", required=True, help="scTE genome index (.exclusive.idx)")
    parser.add_argument("--cb-tag", default="CB", help="Cell barcode BAM tag (default: CB)")
    parser.add_argument("--umi-tag", default="UB", help="UMI BAM tag (default: UB)")
    parser.add_argument("--threads", type=int, default=1, help="Threads")
    parser.add_argument("--scte-bin", default="scTE", help="scTE binary path")
    parser.add_argument("--sample-id", default="", help="Sample ID for .obs annotation")
    args = parser.parse_args()

    logger = setup_logger("scTE_quantify")

    # Run scTE in a temp directory (scTE writes output files to -o dir)
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = os.path.join(tmpdir, "scTE_out")
        os.makedirs(outdir)
        run_scTE(args.input, outdir, args.index, args.cb_tag, args.umi_tag, args.threads, args.scte_bin, logger, cwd=tmpdir)

        # Find the CSV output (scTE names it based on input BAM name)
        csv_files = [f for f in os.listdir(outdir) if f.endswith((".csv", ".csv.gz"))]
        if not csv_files:
            csv_files = [f for f in os.listdir(tmpdir) if f.endswith((".csv", ".csv.gz"))]
            if csv_files:
                csv_path = os.path.join(tmpdir, csv_files[0])
            else:
                raise FileNotFoundError(f"No CSV output found in {outdir} or {tmpdir}")
        else:
            csv_path = os.path.join(outdir, csv_files[0])

        adata = scTE_csv_to_h5ad(csv_path, args.sample_id)
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        # Move raw scTE matrix alongside h5ad
        csv_base = csv_files[0]
        raw_csv_out = os.path.join(os.path.dirname(args.output), csv_base)
        shutil.move(csv_path, raw_csv_out)
        logger.info(f"Saved raw scTE matrix to {raw_csv_out}")
        # Move scTE_out.stat if present
        stat_src = os.path.join(outdir, "scTE_out.stat")
        if not os.path.exists(stat_src):
            stat_src = os.path.join(tmpdir, "scTE_out.stat")
        if os.path.exists(stat_src):
            stat_dst = os.path.join(os.path.dirname(args.output), "scTE_out.stat")
            shutil.move(stat_src, stat_dst)
            logger.info(f"Saved scTE_out.stat to {stat_dst}")
        adata.write_h5ad(args.output)
        logger.info(f"Wrote {adata.n_obs} cells x {adata.n_vars} TEs to {args.output}")


if __name__ == "__main__":
    main()
