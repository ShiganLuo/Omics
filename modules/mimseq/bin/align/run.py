#!/usr/bin/env python3
"""Wrapper for mimseq alignment module.

Aligns reads using GSNAP and returns BAM files and coverage data.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

# Add mimseq source to path (parent of mimseq package)
MIMSEQ_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
if MIMSEQ_PARENT not in sys.path:
    sys.path.insert(0, MIMSEQ_PARENT)

# Add bin/ to path for serialize module
BIN_DIR = os.path.dirname(os.path.dirname(__file__))
if BIN_DIR not in sys.path:
    sys.path.insert(0, BIN_DIR)

from serialize import load_tRNA_state, save_align_state, load_json
# Import mimseq as a package
import mimseq.tRNAmap as tRNAmap_module

log = logging.getLogger(__name__)

# Plot formats forwarded to the vendored mimseq plotting code through the
# MIMSEQ_PLOT_FORMATS env var (read around ggsave()/savefig in mimseq/*.R and
# mimseq/tRNAmap.py; default "pdf" preserves upstream behaviour). Valid set =
# the R ggsave device set (every entry is also a matplotlib savefig format).
PLOT_FORMATS = ("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg")


def export_plot_formats(formats):
    """Publish requested plot formats to the vendored plotting code.

    Save formats. Default ``["pdf"]``.
    """
    fmts = list(formats) if formats else ["pdf"]
    bad = [f for f in fmts if f not in PLOT_FORMATS]
    if bad:
        raise ValueError(
            f"Unsupported image format(s): {', '.join(bad)}; "
            f"valid: {', '.join(PLOT_FORMATS)}")
    os.environ["MIMSEQ_PLOT_FORMATS"] = ",".join(fmts)


def run_align(
    sample_data: str,
    name: str,
    out: str,
    threads: int,
    keep_temp: bool,
    mismatches: float,
    remap: bool,
) -> None:
    """Run alignment using GSNAP."""

    # Ensure out has trailing /
    if not out.endswith("/"):
        out = out + "/"
    os.makedirs(out, exist_ok=True)

    # Set up logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)-5.5s] %(message)s",
        level=logging.INFO,
        handlers=[
            logging.FileHandler(os.path.join(out, "align.log")),
            logging.StreamHandler(),
        ],
    )

    log.info("Starting alignment module")

    # Load tRNA state
    tRNA_state = load_tRNA_state(out)
    snp_tolerance = tRNA_state["snp_tolerance"]

    # Load GSNAP indices
    gsnap_indices = load_json(os.path.join(out, "state", "gsnap_indices.json"))
    genome_index_path = gsnap_indices["genome_index_path"]
    genome_index_name = gsnap_indices["genome_index_name"]
    snp_index_path = gsnap_indices["snp_index_path"]
    snp_index_name = gsnap_indices["snp_index_name"]

    # Run alignment
    map_round = 1
    bams_list, coverageData = tRNAmap_module.mainAlign(
        sample_data, name,
        genome_index_path, genome_index_name,
        snp_index_path, snp_index_name,
        out, threads, snp_tolerance, keep_temp,
        mismatches, map_round, remap,
    )

    # Save alignment state
    save_align_state(
        outdir=out,
        bams_list=bams_list,
        coverageData=coverageData,
    )

    log.info("Alignment module completed successfully")


def main():
    parser = argparse.ArgumentParser(description="mimseq alignment module")
    parser.add_argument("--sample-data", required=True, help="Sample data file")
    parser.add_argument("--name", required=True, help="Experiment name")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files")
    parser.add_argument("--mismatches", type=float, default=0.075, help="Maximum mismatches")
    parser.add_argument("--remap", action="store_true", help="Enable remapping")

    parser.add_argument("-f", "--format", action="append", dest="formats",
                        metavar="FMT",
                        help="Image format (png/pdf/svg/...). Repeatable. "
                             "Default: pdf.")

    args = parser.parse_args()
    export_plot_formats(args.formats)

    run_align(
        sample_data=args.sample_data,
        name=args.name,
        out=args.out,
        threads=args.threads,
        keep_temp=args.keep_temp,
        mismatches=args.mismatches,
        remap=args.remap,
    )


if __name__ == "__main__":
    main()
