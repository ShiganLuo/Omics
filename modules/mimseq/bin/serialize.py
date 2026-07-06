#!/usr/bin/env python3
"""Serialization utilities for mimseq intermediate state.

Saves/loads Python objects (dict, defaultdict, etc.) to/from pickle files
so that independent Snakemake modules can communicate state.
"""

from __future__ import annotations

import json
import logging
import os
import pickle
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Optional

log = logging.getLogger(__name__)


def save_pickle(obj: Any, path: str) -> None:
    """Save a Python object to a pickle file."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    log.info(f"Saved pickle: {path} ({os.path.getsize(path)} bytes)")


def load_pickle(path: str) -> Any:
    """Load a Python object from a pickle file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Pickle file not found: {path}")
    with open(path, "rb") as f:
        obj = pickle.load(f)
    log.info(f"Loaded pickle: {path}")
    return obj


def save_json(obj: Any, path: str) -> None:
    """Save a JSON-serializable object to a JSON file."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
    log.info(f"Saved JSON: {path}")


def load_json(path: str) -> Any:
    """Load a JSON object from a file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"JSON file not found: {path}")
    with open(path, "r") as f:
        obj = json.load(f)
    log.info(f"Loaded JSON: {path}")
    return obj


# Specific save/load functions for mimseq state

def save_tRNA_state(
    outdir: str,
    tRNA_dict: dict,
    cluster_dict: dict,
    coverage_bed: str,
    snp_tolerance: bool,
    mismatch_dict: dict,
    insert_dict: dict,
    del_dict: dict,
    mod_lists: dict,
    Inosine_lists: dict,
    Inosine_clusters: dict,
    cluster_perPos_mismatchMembers: dict,
) -> None:
    """Save tRNA tools state to pickle files."""
    state_dir = os.path.join(outdir, "state")
    os.makedirs(state_dir, exist_ok=True)

    save_pickle(tRNA_dict, os.path.join(state_dir, "tRNA_dict.pkl"))
    save_pickle(cluster_dict, os.path.join(state_dir, "cluster_dict.pkl"))
    save_pickle(mismatch_dict, os.path.join(state_dir, "mismatch_dict.pkl"))
    save_pickle(insert_dict, os.path.join(state_dir, "insert_dict.pkl"))
    save_pickle(del_dict, os.path.join(state_dir, "del_dict.pkl"))
    save_pickle(mod_lists, os.path.join(state_dir, "mod_lists.pkl"))
    save_pickle(Inosine_lists, os.path.join(state_dir, "Inosine_lists.pkl"))
    save_pickle(Inosine_clusters, os.path.join(state_dir, "Inosine_clusters.pkl"))
    save_pickle(cluster_perPos_mismatchMembers, os.path.join(state_dir, "cluster_perPos_mismatchMembers.pkl"))

    # Save simple values as JSON
    meta = {
        "coverage_bed": coverage_bed,
        "snp_tolerance": snp_tolerance,
    }
    save_json(meta, os.path.join(state_dir, "tRNA_meta.json"))


def load_tRNA_state(outdir: str) -> dict:
    """Load tRNA tools state from pickle files."""
    state_dir = os.path.join(outdir, "state")

    meta = load_json(os.path.join(state_dir, "tRNA_meta.json"))

    return {
        "tRNA_dict": load_pickle(os.path.join(state_dir, "tRNA_dict.pkl")),
        "cluster_dict": load_pickle(os.path.join(state_dir, "cluster_dict.pkl")),
        "mismatch_dict": load_pickle(os.path.join(state_dir, "mismatch_dict.pkl")),
        "insert_dict": load_pickle(os.path.join(state_dir, "insert_dict.pkl")),
        "del_dict": load_pickle(os.path.join(state_dir, "del_dict.pkl")),
        "mod_lists": load_pickle(os.path.join(state_dir, "mod_lists.pkl")),
        "Inosine_lists": load_pickle(os.path.join(state_dir, "Inosine_lists.pkl")),
        "Inosine_clusters": load_pickle(os.path.join(state_dir, "Inosine_clusters.pkl")),
        "cluster_perPos_mismatchMembers": load_pickle(os.path.join(state_dir, "cluster_perPos_mismatchMembers.pkl")),
        "coverage_bed": meta["coverage_bed"],
        "snp_tolerance": meta["snp_tolerance"],
    }


def save_align_state(
    outdir: str,
    bams_list: list,
    coverageData: dict,
) -> None:
    """Save alignment state to pickle files."""
    state_dir = os.path.join(outdir, "state")
    os.makedirs(state_dir, exist_ok=True)

    save_pickle(bams_list, os.path.join(state_dir, "bams_list.pkl"))
    save_pickle(coverageData, os.path.join(state_dir, "coverageData.pkl"))


def load_align_state(outdir: str) -> dict:
    """Load alignment state from pickle files."""
    state_dir = os.path.join(outdir, "state")

    return {
        "bams_list": load_pickle(os.path.join(state_dir, "bams_list.pkl")),
        "coverageData": load_pickle(os.path.join(state_dir, "coverageData.pkl")),
    }


def save_cluster_state(
    outdir: str,
    splitBool_new: dict,
    unique_isodecoderMMs_new: dict,
    notSplit_cov_posInfo: dict,
    notSplit_mods_posInfo: dict,
    isodecoder_sizes: dict,
    unsplitCluster_lookup: dict,
) -> None:
    """Save cluster state to pickle files."""
    state_dir = os.path.join(outdir, "state")
    os.makedirs(state_dir, exist_ok=True)

    save_pickle(splitBool_new, os.path.join(state_dir, "splitBool_new.pkl"))
    save_pickle(unique_isodecoderMMs_new, os.path.join(state_dir, "unique_isodecoderMMs_new.pkl"))
    save_pickle(notSplit_cov_posInfo, os.path.join(state_dir, "notSplit_cov_posInfo.pkl"))
    save_pickle(notSplit_mods_posInfo, os.path.join(state_dir, "notSplit_mods_posInfo.pkl"))
    save_pickle(isodecoder_sizes, os.path.join(state_dir, "isodecoder_sizes.pkl"))
    save_pickle(unsplitCluster_lookup, os.path.join(state_dir, "unsplitCluster_lookup.pkl"))


def load_cluster_state(outdir: str) -> dict:
    """Load cluster state from pickle files."""
    state_dir = os.path.join(outdir, "state")

    return {
        "splitBool_new": load_pickle(os.path.join(state_dir, "splitBool_new.pkl")),
        "unique_isodecoderMMs_new": load_pickle(os.path.join(state_dir, "unique_isodecoderMMs_new.pkl")),
        "notSplit_cov_posInfo": load_pickle(os.path.join(state_dir, "notSplit_cov_posInfo.pkl")),
        "notSplit_mods_posInfo": load_pickle(os.path.join(state_dir, "notSplit_mods_posInfo.pkl")),
        "isodecoder_sizes": load_pickle(os.path.join(state_dir, "isodecoder_sizes.pkl")),
        "unsplitCluster_lookup": load_pickle(os.path.join(state_dir, "unsplitCluster_lookup.pkl")),
    }


def get_sample_dir(out_dir: str, bam_path: str) -> str:
    """Return per-sample subdirectory under out_dir/samples/<sample>/.

    Creates the directory if it doesn't exist. The sample name is derived
    from the BAM filename by stripping the .single.fq.gz.uniq.bam suffix.
    """
    import re
    bam_name = os.path.basename(bam_path)
    sample_name = re.sub(r'\.single\.fq\.gz\..*$', '', bam_name)
    sample_dir = os.path.join(out_dir, "samples", sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    return sample_dir
