#!/usr/bin/env python3
"""Wrapper for mimseq DESeq2 differential expression module.

Runs DESeq2 analysis for differential expression.
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

# Add mimseq source to path
MIMSEQ_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
if MIMSEQ_PARENT not in sys.path:
    sys.path.insert(0, MIMSEQ_PARENT)

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


def run_deseq(
    out: str,
    control_cond: str,
    p_adj: float,
    mito_trnas: str,
) -> None:
    """Run DESeq2 differential expression analysis."""

    # Ensure out has trailing /
    if not out.endswith("/"):
        out = out + "/"
    os.makedirs(out, exist_ok=True)

    # Set up logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)-5.5s] %(message)s",
        level=logging.INFO,
        handlers=[
            logging.FileHandler(os.path.join(out, "deseq.log")),
            logging.StreamHandler(),
        ],
    )

    log.info("Starting DESeq2 module")

    # Run DESeq2 R script
    script_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "mimseq")
    deseq_cmd = [
        "Rscript",
        os.path.join(script_path, "deseq.R"),
        out,
        control_cond,
        str(p_adj),
        str(mito_trnas),
    ]
    try:
        process = subprocess.Popen(deseq_cmd, stdout=subprocess.PIPE)
        while True:
            line = process.stdout.readline()
            if not line:
                break
            line = line.decode("utf-8")
            log.info(line.rstrip())
        exitcode = process.wait()
    except subprocess.CalledProcessError:
        log.error("Error running DESeq2")
        raise

    log.info("DESeq2 module completed successfully")


def main():
    parser = argparse.ArgumentParser(description="mimseq DESeq2 module")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--control-cond", default="", help="Control condition")
    parser.add_argument("--p-adj", type=float, default=0.05, help="Adjusted p-value threshold")
    parser.add_argument("--mito-trnas", default="", help="Mitochondrial tRNA fasta file")

    parser.add_argument("-f", "--format", action="append", dest="formats",
                        metavar="FMT",
                        help="Image format (png/pdf/svg/...). Repeatable. "
                             "Default: pdf.")

    args = parser.parse_args()
    export_plot_formats(args.formats)

    run_deseq(
        out=args.out,
        control_cond=args.control_cond,
        p_adj=args.p_adj,
        mito_trnas=args.mito_trnas,
    )


if __name__ == "__main__":
    main()
