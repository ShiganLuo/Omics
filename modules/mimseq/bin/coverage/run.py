#!/usr/bin/env python3
"""Wrapper for mimseq coverage module.

Calculates coverage and generates coverage plots.
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

# Add bin/ to path for serialize module
BIN_DIR = os.path.dirname(os.path.dirname(__file__))
if BIN_DIR not in sys.path:
    sys.path.insert(0, BIN_DIR)

from serialize import load_align_state, load_cluster_state, load_pickle
import mimseq.getCoverage as getCoverage_module
import mimseq.mmQuant as mmQuant_module

log = logging.getLogger(__name__)


def run_coverage(
    out: str,
    control_cond: str,
    cca: bool,
    double_cca: bool,
    mito_trnas: str,
    misinc_thresh: float,
    mod_sites: str,
    cons_pos_list: str,
) -> None:
    """Run coverage calculation and plotting."""

    # Ensure out has trailing /
    if not out.endswith("/"):
        out = out + "/"
    os.makedirs(out, exist_ok=True)

    # Set up logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)-5.5s] %(message)s",
        level=logging.INFO,
        handlers=[
            logging.FileHandler(os.path.join(out, "coverage.log")),
            logging.StreamHandler(),
        ],
    )

    log.info("Starting coverage module")

    # Load states
    align_state = load_align_state(out)
    coverageData = align_state["coverageData"]

    # Load filtered coverage if available
    filtered_cov_path = os.path.join(out, "state", "filtered_cov.pkl")
    filtered_cov = load_pickle(filtered_cov_path) if os.path.exists(filtered_cov_path) else []

    # Load cluster state
    cluster_state = load_cluster_state(out)
    unsplitCluster_lookup = cluster_state["unsplitCluster_lookup"]

    # Run coverage calculation
    sorted_aa = getCoverage_module.getCoverage(coverageData, out, control_cond, filtered_cov, unsplitCluster_lookup)

    # Plot coverage
    getCoverage_module.plotCoverage(out, mito_trnas, sorted_aa)

    # CCA analysis
    if cca:
        cca_dir = os.path.join(out, "CCAanalysis")
        cca_file = os.path.join(cca_dir, "AlignedDinucProportions.csv")
        if os.path.exists(cca_file):
            mmQuant_module.plotCCA(out, double_cca)
        else:
            log.warning(f"CCA data file not found ({cca_file}), skipping CCA plot")

    # Modification plots (non-fatal — mods data may not exist yet with remap=True)
    script_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "mimseq")
    mods_dir = os.path.join(out, "mods")
    if os.path.isdir(mods_dir) and os.listdir(mods_dir):
        modplot_cmd = [
            "Rscript",
            os.path.join(script_path, "modPlot.R"),
            out,
            str(mod_sites),
            str(cons_pos_list),
            str(misinc_thresh),
            str(mito_trnas),
            control_cond,
        ]
        try:
            process = subprocess.Popen(modplot_cmd, stdout=subprocess.PIPE)
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                line = line.decode("utf-8")
                log.info(line.rstrip())
            exitcode = process.wait()
        except subprocess.CalledProcessError:
            log.warning("Error plotting modifications (non-fatal)")
    else:
        log.warning(f"Mods directory empty ({mods_dir}), skipping modification plots")

    log.info("Coverage module completed successfully")


def main():
    parser = argparse.ArgumentParser(description="mimseq coverage module")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--control-cond", default="", help="Control condition")
    parser.add_argument("--cca", action="store_true", help="Enable CCA analysis")
    parser.add_argument("--double-cca", action="store_true", help="Enable double CCA analysis")
    parser.add_argument("--mito-trnas", default="", help="Mitochondrial tRNA fasta file")
    parser.add_argument("--misinc-thresh", type=float, default=0.1, help="Misincorporation threshold")
    parser.add_argument("--mod-sites", default="", help="Modification sites")
    parser.add_argument("--cons-pos-list", default="", help="Conservation positions list")

    args = parser.parse_args()

    run_coverage(
        out=args.out,
        control_cond=args.control_cond,
        cca=args.cca,
        double_cca=args.double_cca,
        mito_trnas=args.mito_trnas,
        misinc_thresh=args.misinc_thresh,
        mod_sites=args.mod_sites,
        cons_pos_list=args.cons_pos_list,
    )


if __name__ == "__main__":
    main()
