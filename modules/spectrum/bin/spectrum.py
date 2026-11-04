#!/usr/bin/env python
"""Mutation spectrum analysis: trinucleotide-context counts from somatic VCFs.

Loads a reference genome, classifies each SNV into its pyrimidine-strand
trinucleotide context, aggregates counts per sample group, and renders a
stacked bar plot, a heatmap, and a counts table.

Usage:
    python spectrum.py \
        --vcf-map experiment_vcf_map.tsv \
        --group-map experiment_group_map.tsv \
        --ref-fasta genome.fa \
        --output-prefix somatic_spectrum
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

_BIN_DIR = Path(__file__).resolve().parent
_ROOT_DIR = _BIN_DIR.parent.parent.parent
for _p in (str(_BIN_DIR), str(_ROOT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from src.common.util.LogUtil import setup_logger  # noqa: E402

import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402
from Bio import SeqIO  # noqa: E402

logger = setup_logger("Spectrum")

# Canonical order of the six pyrimidine-strand base substitutions.
BASE_MUT_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def load_reference_genome(fasta_file):
    """Load a reference genome into a dict of contig name -> uppercase sequence."""
    logger.info("Loading reference genome...")
    ref_seq = {record.id: str(record.seq).upper() for record in SeqIO.parse(fasta_file, "fasta")}
    logger.info(f"Reference genome loaded: {len(ref_seq)} contigs.")
    return ref_seq


def get_trinucleotide(ref_seq, chrom, pos, ref, alt):
    """Return the trinucleotide context of a SNV, normalized to the pyrimidine strand.

    The context is formatted as ``<5'base>[<ref>><alt>]<3'base>``; purine
    reference alleles are reverse-complemented so all contexts are reported
    on the pyrimidine (C/T) strand. Returns None when the flanking context
    is unavailable.
    """
    seq = ref_seq.get(chrom)
    if not seq:
        return None
    pos0 = pos - 1
    if pos0 < 1 or pos0 + 1 >= len(seq):
        return None
    tri = seq[pos0 - 1:pos0 + 2]
    ref_base = ref.upper()
    alt_base = alt.upper()
    # Normalize to the pyrimidine strand.
    if ref_base in ["A", "G"]:
        complement = str.maketrans("ACGT", "TGCA")
        tri = tri.translate(complement)[::-1]
        ref_base = ref_base.translate(complement)
        alt_base = alt_base.translate(complement)
    return f"{tri[0]}[{ref_base}>{alt_base}]{tri[2]}"


def process_vcf(vcf_file, ref_seq):
    """Parse a VCF (first 6 columns) and count trinucleotide contexts of SNVs."""
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL"]
    df = pd.read_csv(vcf_file, sep="\t", comment="#", header=None, usecols=range(6), names=cols)
    df = df[df["ALT"].str.len() == 1]
    df["trinuc"] = df.apply(lambda x: get_trinucleotide(ref_seq, x["CHROM"], int(x["POS"]), x["REF"], x["ALT"]), axis=1)
    df = df.dropna(subset=["trinuc"])
    return Counter(df["trinuc"])


def build_group_dataframe(all_counts, sample_groups):
    """Aggregate per-sample trinucleotide counters into a group-by-context DataFrame."""
    group_counts = {}
    for sample, counts in all_counts.items():
        group = sample_groups[sample]
        if group not in group_counts:
            group_counts[group] = Counter()
        group_counts[group] += counts

    # Order rows by the six base substitution classes.
    tri_list = []
    for base_mut in BASE_MUT_ORDER:
        for tri in sorted(set().union(*[c.keys() for c in group_counts.values()])):
            if f"[{base_mut}]" in tri:
                tri_list.append(tri)

    tri_df = pd.DataFrame(0, index=tri_list, columns=list(group_counts.keys()))
    for group, counts in group_counts.items():
        for tri, num in counts.items():
            if tri in tri_df.index:
                tri_df.loc[tri, group] = num
    return tri_df


def plot_stacked_bar(tri_df, output_file):
    """Render the per-group trinucleotide mutation spectrum as a stacked bar plot."""
    tri_df.plot(kind="bar", stacked=True, figsize=(20, 6))
    plt.ylabel("Mutation count")
    plt.xlabel("Trinucleotide context")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_heatmap(tri_df, output_file):
    """Render the group-by-context counts as a heatmap."""
    plt.figure(figsize=(12, 8))
    sns.heatmap(tri_df, cmap="Reds", annot=False)
    plt.ylabel("Trinucleotide context")
    plt.xlabel("Group")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def mutation_spectrum_analysis(vcf_files, sample_groups, ref_fasta, output_prefix):
    """Run the full spectrum analysis and write csv + stacked bar + heatmap."""
    ref_seq = load_reference_genome(ref_fasta)
    all_counts = {}
    for sample, vcf in vcf_files.items():
        logger.info(f"Processing {sample} ...")
        all_counts[sample] = process_vcf(vcf, ref_seq)

    tri_df = build_group_dataframe(all_counts, sample_groups)
    tri_df.to_csv(f"{output_prefix}.csv")
    plot_stacked_bar(tri_df, f"{output_prefix}_stacked_bar.png")
    plot_heatmap(tri_df, f"{output_prefix}_heatmap.png")
    logger.info(f"Wrote {output_prefix}.csv, {output_prefix}_stacked_bar.png, {output_prefix}_heatmap.png")
    return tri_df


def read_mapping_tsv(file_path, key_col, value_col):
    """Read a two-column TSV mapping (header required) into a dict."""
    mapping = {}
    with open(file_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if key_col not in reader.fieldnames or value_col not in reader.fieldnames:
            raise ValueError(f"Missing required columns: {key_col}, {value_col}")
        for row in reader:
            key = row[key_col].strip()
            value = row[value_col].strip()
            if key:
                mapping[key] = value
    return mapping


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Mutation spectrum analysis")
    parser.add_argument("--vcf-map", required=True, help="TSV with columns: experiment_sample_id, vcf")
    parser.add_argument("--group-map", required=True, help="TSV with columns: experiment_sample_id, group")
    parser.add_argument("--ref-fasta", required=True, help="Reference FASTA")
    parser.add_argument("--output-prefix", required=True, help="Output prefix")
    parser.add_argument("--log", default=None, help="Optional log file (in addition to stdout)")
    return parser.parse_args()


def main():
    """Entry point: read maps, run the analysis, write outputs."""
    args = parse_args()
    if args.log:
        _log_path = Path(args.log)
        _log_path.parent.mkdir(parents=True, exist_ok=True)
        _lg = setup_logger("Spectrum", log_file=str(_log_path))
        logger.handlers = _lg.handlers
    vcf_files = read_mapping_tsv(args.vcf_map, "experiment_sample_id", "vcf")
    sample_groups = read_mapping_tsv(args.group_map, "experiment_sample_id", "group")
    mutation_spectrum_analysis(vcf_files, sample_groups, args.ref_fasta, args.output_prefix)


if __name__ == "__main__":
    main()
