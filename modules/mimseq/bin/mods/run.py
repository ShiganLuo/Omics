#!/usr/bin/env python3
"""Wrapper for mimseq modification quantification module.

Generates modification tables and misincorporation data.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
from pathlib import Path

MIMSEQ_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
if MIMSEQ_PARENT not in sys.path:
    sys.path.insert(0, MIMSEQ_PARENT)

BIN_DIR = os.path.dirname(os.path.dirname(__file__))
if BIN_DIR not in sys.path:
    sys.path.insert(0, BIN_DIR)

from serialize import load_tRNA_state, load_align_state, load_cluster_state, save_pickle
import mimseq.mmQuant as mmQuant_module
import mimseq.crosstalks as crosstalks_module
import mimseq.ssAlign as ssAlign_module

log = logging.getLogger(__name__)

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


def run_mods(
    out: str,
    name: str,
    threads: int,
    min_cov: float,
    misinc_thresh: float,
    cca: bool,
    remap: bool,
    crosstalks: bool,
) -> None:
    """Run modification quantification."""

    if not out.endswith("/"):
        out = out + "/"
    os.makedirs(out, exist_ok=True)

    logging.basicConfig(
        format="%(asctime)s [%(levelname)-5.5s] %(message)s",
        level=logging.INFO,
        handlers=[
            logging.FileHandler(os.path.join(out, "mods.log")),
            logging.StreamHandler(),
        ],
    )

    log.info("Starting modification quantification module")

    tRNA_state = load_tRNA_state(out)
    align_state = load_align_state(out)
    cluster_state = load_cluster_state(out)

    stk_files = glob.glob(os.path.join(out, "*_align.stk"))
    if not stk_files:
        raise FileNotFoundError(f"No *_align.stk file found in {out}")
    ssAlign_module.stkname = stk_files[0]
    log.info(f"Set stkname: {stk_files[0]}")

    new_mods, new_Inosines, filtered_cov, filter_warning, unsplitCluster_lookup, readRef_unsplit_newNames = \
        mmQuant_module.generateModsTable(
            align_state["coverageData"],
            out,
            name,
            threads,
            min_cov,
            tRNA_state["mismatch_dict"],
            tRNA_state["insert_dict"],
            tRNA_state["del_dict"],
            tRNA_state["cluster_dict"],
            cca,
            remap,
            misinc_thresh,
            tRNA_state["mod_lists"],
            tRNA_state["Inosine_lists"],
            tRNA_state["tRNA_dict"],
            tRNA_state["Inosine_clusters"],
            cluster_state["unique_isodecoderMMs_new"],
            cluster_state["splitBool_new"],
            cluster_state["isodecoder_sizes"],
            cluster_state["unsplitCluster_lookup"],
            True,
            crosstalks,
        )

    if crosstalks:
        log.info("Analyzing crosstalks between pairs of modifications and modification-charging...")
        crosstalks_module.crosstalks_wrapper(os.path.join(out, "single_read_data"), misinc_thresh, threads)

    save_pickle(new_mods, os.path.join(out, "state", "new_mods.pkl"))
    save_pickle(new_Inosines, os.path.join(out, "state", "new_Inosines.pkl"))
    save_pickle(filtered_cov, os.path.join(out, "state", "filtered_cov.pkl"))
    save_pickle(readRef_unsplit_newNames, os.path.join(out, "state", "readRef_unsplit_newNames.pkl"))

    log.info("Modification quantification module completed successfully")


def main():
    parser = argparse.ArgumentParser(description="mimseq modification quantification module")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--name", required=True, help="Experiment name")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--min-cov", type=float, default=0.0005, help="Minimum coverage")
    parser.add_argument("--misinc-thresh", type=float, default=0.1, help="Misincorporation threshold")
    parser.add_argument("--cca", action="store_true", help="Enable CCA analysis")
    parser.add_argument("--remap", action="store_true", help="Enable remapping")
    parser.add_argument("--crosstalks", action="store_true", help="Enable crosstalks analysis")

    parser.add_argument("-f", "--format", action="append", dest="formats",
                        metavar="FMT",
                        help="Image format (png/pdf/svg/...). Repeatable. "
                             "Default: pdf.")

    args = parser.parse_args()
    export_plot_formats(args.formats)

    run_mods(
        out=args.out,
        name=args.name,
        threads=args.threads,
        min_cov=args.min_cov,
        misinc_thresh=args.misinc_thresh,
        cca=args.cca,
        remap=args.remap,
        crosstalks=args.crosstalks,
    )


if __name__ == "__main__":
    main()
