#!/usr/bin/env python3
"""Wrapper for mimseq clusters module.

Splits isodecoders and performs cluster deconvolution.
"""

from __future__ import annotations

import argparse
import copy
import glob
import logging
import os
import sys
from collections import defaultdict

# Add mimseq source to path (parent of mimseq package)
MIMSEQ_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
if MIMSEQ_PARENT not in sys.path:
    sys.path.insert(0, MIMSEQ_PARENT)

# Add bin/ to path for serialize module
BIN_DIR = os.path.dirname(os.path.dirname(__file__))
if BIN_DIR not in sys.path:
    sys.path.insert(0, BIN_DIR)

from serialize import load_tRNA_state, load_align_state, save_cluster_state
import mimseq.splitClusters as splitClusters_module
import mimseq.ssAlign as ssAlign_module

log = logging.getLogger(__name__)


def run_clusters(
    out: str,
    name: str,
    cluster: bool,
    cluster_id: float,
    cov_diff: float,
    threads: int,
) -> None:
    """Run cluster splitting and deconvolution."""

    # Ensure out has trailing /
    if not out.endswith("/"):
        out = out + "/"
    os.makedirs(out, exist_ok=True)

    # Set up logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)-5.5s] %(message)s",
        level=logging.INFO,
        handlers=[
            logging.FileHandler(os.path.join(out, "clusters.log")),
            logging.StreamHandler(),
        ],
    )

    log.info("Starting clusters module")

    # Set ssAlign global stkname (required by tRNAclassifier -> structureParser)
    stk_files = glob.glob(os.path.join(out, "*_align.stk"))
    if stk_files:
        ssAlign_module.stkname = stk_files[0]
        log.info(f"Set stkname: {stk_files[0]}")
    else:
        log.warning("No .stk file found, structureParser may fail")

    # Load tRNA state
    tRNA_state = load_tRNA_state(out)
    tRNA_dict = tRNA_state["tRNA_dict"]
    cluster_dict = tRNA_state["cluster_dict"]
    coverage_bed = tRNA_state["coverage_bed"]
    cluster_perPos_mismatchMembers = tRNA_state["cluster_perPos_mismatchMembers"]

    # Load alignment state
    align_state = load_align_state(out)
    coverageData = align_state["coverageData"]

    # Initialize variables
    unique_isodecoderMMs = defaultdict(dict)
    unique_isodecoderMMs_new = defaultdict(dict)
    splitBool = defaultdict(set)
    splitBool_new = defaultdict(set)
    unsplitCluster_lookup = defaultdict()
    notSplit_cov_posInfo = defaultdict(set)
    notSplit_mods_posInfo = defaultdict(set)
    isodecoder_sizes = {}

    # Run cluster splitting
    if cluster and cluster_id != 1:
        cluster_dict2 = copy.deepcopy(cluster_dict)
        unique_isodecoderMMs, splitBool, notSplit_mods_posInfo = splitClusters_module.splitIsodecoder(
            cluster_perPos_mismatchMembers,
            tRNA_state["insert_dict"],
            tRNA_state["del_dict"],
            tRNA_dict,
            cluster_dict2,
            out,
            name,
        )
        splitBool_new, unique_isodecoderMMs_new, notSplit_cov_posInfo = splitClusters_module.unsplitClustersCov(
            coverageData,
            coverage_bed,
            unique_isodecoderMMs,
            splitBool,
            threads,
            1,  # map_round
            cov_diff,
        )
        isodecoder_sizes, unsplitCluster_lookup = splitClusters_module.getDeconvSizes(
            splitBool_new,
            tRNA_dict,
            cluster_dict,
            unique_isodecoderMMs_new,
        )
        splitClusters_module.writeDeconvTranscripts(out, name, tRNA_dict, isodecoder_sizes, cluster)
    elif cluster and cluster_id == 1:
        isodecoder_sizes = {iso: len(members) for iso, members in cluster_dict.items()}
        splitClusters_module.writeIsodecoderTranscripts(out, name, cluster_dict, tRNA_dict)
    elif not cluster:
        isodecoder_sizes = splitClusters_module.getIsodecoderSizes(out, name, tRNA_dict)

    # Save cluster state
    save_cluster_state(
        outdir=out,
        splitBool_new=splitBool_new,
        unique_isodecoderMMs_new=unique_isodecoderMMs_new,
        notSplit_cov_posInfo=notSplit_cov_posInfo,
        notSplit_mods_posInfo=notSplit_mods_posInfo,
        isodecoder_sizes=isodecoder_sizes,
        unsplitCluster_lookup=unsplitCluster_lookup,
    )

    log.info("Clusters module completed successfully")


def main():
    parser = argparse.ArgumentParser(description="mimseq clusters module")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--name", required=True, help="Experiment name")
    parser.add_argument("--cluster", action="store_true", help="Enable clustering")
    parser.add_argument("--cluster-id", type=float, default=0.97, help="Cluster identity threshold")
    parser.add_argument("--cov-diff", type=float, default=0.5, help="Coverage difference threshold")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")

    args = parser.parse_args()

    run_clusters(
        out=args.out,
        name=args.name,
        cluster=args.cluster,
        cluster_id=args.cluster_id,
        cov_diff=args.cov_diff,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
