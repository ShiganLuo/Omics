#!/usr/bin/env python3
"""Wrapper for mimseq tRNAtools module.

Parses tRNA sequences and modifications, generates SNP index and GSNAP indices.
Saves intermediate state to pickle files for downstream modules.
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

from serialize import save_tRNA_state
# Import mimseq as a package
import mimseq.tRNAtools as tRNAtools_module
import mimseq.ssAlign as ssAlign_module

log = logging.getLogger(__name__)


def run_tRNAtools(
    trnas: str,
    trnaout: str,
    mito_trnas: str,
    plastid_trnas: str,
    modifications: str,
    name: str,
    out: str,
    double_cca: bool,
    threads: int,
    snp_tolerance: bool,
    cluster: bool,
    cluster_id: float,
    posttrans: bool,
    pretrnas: bool,
    local_mod: bool,
    species: str,
) -> None:
    """Run tRNA tools: parse tRNA, generate SNP index and GSNAP indices."""

    # Ensure out has trailing /
    if not out.endswith("/"):
        out = out + "/"
    os.makedirs(out, exist_ok=True)

    # Set up logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)-5.5s] %(message)s",
        level=logging.INFO,
        handlers=[
            logging.FileHandler(os.path.join(out, "tRNAtools.log")),
            logging.StreamHandler(),
        ],
    )

    log.info("Starting tRNAtools module")

    # If species is specified, set reference paths from mimseq data
    if species:
        mimseq_data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "mimseq", "data")
        species_refs = {
            "Hsap": ("hg38-eColitK/hg38-tRNAs-all.fa", "hg38-eColitK/hg38-tRNAs-detailed.out", "hg38-eColitK/hg38-mitotRNAs.fa"),
            "Hsap19": ("hg19-eColitK/hg19-tRNAs-all.fa", "hg19-eColitK/hg19_eschColi-tRNAs.out", "hg19-eColitK/hg19-mitotRNAs.fa"),
            "Mmus": ("mm39-eColitK/mm39-tRNAs-all.fa", "mm39-eColitK/mm39-tRNAs-detailed.out", "mm39-eColitK/mm39-mitotRNAs.fa"),
            "Scer": ("sacCer3-eColitK/sacCer3_eschColitK.fa", "sacCer3-eColitK/sacCer3_eschColi-tRNAs.out", "sacCer3-eColitK/sacCer3-mitotRNAs.fa"),
            "Spom": ("schiPomb-eColitK/schiPomb_eschColitK.fa", "schiPomb-eColitK/schiPomb_eschColi-tRNAs.out", "schiPomb-eColitK/schiPomb-mitotRNAs.fa"),
            "Dmel": ("dm6-eColitK/dm6-tRNAs-all.fa", "dm6-eColitK/dm6_eschColi-tRNAs.out", "dm6-eColitK/dm6-mitotRNAs.fa"),
            "Drer": ("danRer11-eColitK/danRer11-tRNAs-all.fa", "danRer11-eColitK/danRer11_eschColi-tRNAs.out", "danRer11-eColitK/danRer11-mitotRNAs.fa"),
            "Cele": ("ce11-eColitK/ce11-tRNAs-all.fa", "ce11-eColitK/ce11-tRNAs-detailed.out", "ce11-eColitK/ce11-mitotRNAs.fa"),
            "Ecol": ("eschColi-K_12_MG1655-tRNAs/eschColi_K_12_MG1655-tRNAs.fa", "eschColi-K_12_MG1655-tRNAs/eschColi_K_12_MG1655-tRNAs.out", ""),
        }
        if species in species_refs:
            trnas_ref, trnaout_ref, mito_ref = species_refs[species]
            trnas = os.path.join(mimseq_data_dir, trnas_ref)
            trnaout = os.path.join(mimseq_data_dir, trnaout_ref)
            if mito_ref and not mito_trnas:
                mito_trnas = os.path.join(mimseq_data_dir, mito_ref)
            log.info(f"Using built-in references for {species}: {trnas}")
        else:
            raise ValueError(f"Unknown species: {species}")

    # Parse tRNA and modifications, generate SNP index
    (
        coverage_bed,
        snp_tolerance_out,
        mismatch_dict,
        insert_dict,
        del_dict,
        mod_lists,
        Inosine_lists,
        Inosine_clusters,
        tRNA_dict,
        cluster_dict,
        cluster_perPos_mismatchMembers,
    ) = tRNAtools_module.modsToSNPIndex(
        trnas, trnaout, mito_trnas, plastid_trnas,
        modifications, name, out, double_cca, threads,
        snp_tolerance, cluster, cluster_id, posttrans, pretrnas, local_mod,
    )

    # Parse tRNA structure
    ssAlign_module.structureParser()

    # Generate GSNAP indices
    map_round = 1
    genome_index_path, genome_index_name, snp_index_path, snp_index_name = \
        tRNAtools_module.generateGSNAPIndices(species, name, out, map_round, snp_tolerance, cluster)

    # Save state for downstream modules
    save_tRNA_state(
        outdir=out,
        tRNA_dict=tRNA_dict,
        cluster_dict=cluster_dict,
        coverage_bed=coverage_bed,
        snp_tolerance=snp_tolerance_out,
        mismatch_dict=mismatch_dict,
        insert_dict=insert_dict,
        del_dict=del_dict,
        mod_lists=mod_lists,
        Inosine_lists=Inosine_lists,
        Inosine_clusters=Inosine_clusters,
        cluster_perPos_mismatchMembers=cluster_perPos_mismatchMembers,
    )

    # Save GSNAP index paths
    from serialize import save_json
    save_json({
        "genome_index_path": genome_index_path,
        "genome_index_name": genome_index_name,
        "snp_index_path": snp_index_path,
        "snp_index_name": snp_index_name,
    }, os.path.join(out, "state", "gsnap_indices.json"))

    log.info("tRNAtools module completed successfully")


def main():
    parser = argparse.ArgumentParser(description="mimseq tRNAtools module")
    parser.add_argument("--trnas", default="", help="tRNA fasta file")
    parser.add_argument("--trnaout", default="", help="tRNAscan-SE output file")
    parser.add_argument("--mito-trnas", default="", help="Mitochondrial tRNA fasta file")
    parser.add_argument("--plastid-trnas", default="", help="Plastid tRNA fasta file")
    parser.add_argument("--modifications", required=True, help="Modifications table")
    parser.add_argument("--name", required=True, help="Experiment name")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--double-cca", action="store_true", help="Enable double CCA analysis")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--snp-tolerance", action="store_true", help="Enable SNP tolerance")
    parser.add_argument("--cluster", action="store_true", help="Enable clustering")
    parser.add_argument("--cluster-id", type=float, default=0.97, help="Cluster identity threshold")
    parser.add_argument("--posttrans", action="store_true", help="Disable post-transcriptional modification")
    parser.add_argument("--pretrnas", action="store_true", help="Use pre-tRNAs")
    parser.add_argument("--local-mod", action="store_true", help="Use local Modomics data")
    parser.add_argument("--species", required=True, help="Species name")

    args = parser.parse_args()

    run_tRNAtools(
        trnas=args.trnas,
        trnaout=args.trnaout,
        mito_trnas=args.mito_trnas,
        plastid_trnas=args.plastid_trnas,
        modifications=args.modifications,
        name=args.name,
        out=args.out,
        double_cca=args.double_cca,
        threads=args.threads,
        snp_tolerance=args.snp_tolerance,
        cluster=args.cluster,
        cluster_id=args.cluster_id,
        posttrans=args.posttrans,
        pretrnas=args.pretrnas,
        local_mod=args.local_mod,
        species=args.species,
    )


if __name__ == "__main__":
    main()
