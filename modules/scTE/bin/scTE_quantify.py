"""Quantify transposable element expression from single-cell BAM using scTE.

Wraps scTE CLI: reads BAM with cell barcode (CB) and UMI (UB) tags,
produces TE count matrix per cell. Output converted to h5ad for Scanpy.

Usage:
    python scTE_quantify.py \\
        --input sample.bam \\
        --output sample_TE.h5ad \\
        --index /path/to/hg38.exclusive.idx \\
        --cb-tag CB --umi-tag UB \\
        --threads 20
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import anndata as ad
import pandas as pd
import scipy.sparse as sp


def run_scTE(bam, outdir, index, cb_tag, umi_tag, threads, scte_bin, cwd=None):
    """Run scTE quantify on a BAM file."""
    cmd = [
        scte_bin, "-i", bam, "-o", outdir,
        "-x", index,
        "-p", str(threads),
        "-CB", cb_tag, "-UMI", umi_tag,
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.check_call(cmd, cwd=cwd)


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

    # Run scTE in a temp directory (scTE writes output files to -o dir)
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = os.path.join(tmpdir, "scTE_out")
        os.makedirs(outdir)
        run_scTE(args.input, outdir, args.index, args.cb_tag, args.umi_tag, args.threads, args.scte_bin, cwd=tmpdir)

        # Find the CSV output (scTE names it based on input BAM name)
        csv_files = [f for f in os.listdir(outdir) if f.endswith((".csv", ".csv.gz"))]
        if not csv_files:
            # scTE may write directly to the output dir with a fixed name
            csv_files = [f for f in os.listdir(tmpdir) if f.endswith((".csv", ".csv.gz"))]
            if csv_files:
                csv_path = os.path.join(tmpdir, csv_files[0])
            else:
                raise FileNotFoundError(f"No CSV output found in {outdir} or {tmpdir}")
        else:
            csv_path = os.path.join(outdir, csv_files[0])

        adata = scTE_csv_to_h5ad(csv_path, args.sample_id)
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        # Save raw scTE matrix alongside h5ad
        csv_base = csv_files[0]
        raw_csv_out = os.path.join(os.path.dirname(args.output), csv_base)
        shutil.copy2(csv_path, raw_csv_out)
        print(f"Saved raw scTE matrix to {raw_csv_out}")
        # Save scTE_out.stat if present
        stat_src = os.path.join(outdir, "scTE_out.stat")
        if not os.path.exists(stat_src):
            stat_src = os.path.join(tmpdir, "scTE_out.stat")
        if os.path.exists(stat_src):
            stat_dst = os.path.join(os.path.dirname(args.output), "scTE_out.stat")
            shutil.copy2(stat_src, stat_dst)
            print(f"Saved scTE_out.stat to {stat_dst}")
        adata.write_h5ad(args.output)
        print(f"Wrote {adata.n_obs} cells x {adata.n_vars} TEs to {args.output}")


if __name__ == "__main__":
    main()
