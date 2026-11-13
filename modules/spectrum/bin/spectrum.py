#!/usr/bin/env python
"""Mutation spectrum analysis: trinucleotide-context counts from somatic VCFs.

Loads a reference genome, classifies each SNV into its pyrimidine-strand
trinucleotide context, and writes both per-sample and per-group outputs:

  <prefix>_per_sample.csv    96 contexts x samples, raw counts (SBS96 matrix;
                             the input contract for signature_refit.py and
                             signature_nmf.py)
  <prefix>_sample_stats.csv  per-sample SNV counts and drop accounting
  <prefix>.csv               96 contexts x groups (group-aggregated view)
  <prefix>_stacked_bar.png   group spectra stacked bar plot
  <prefix>_heatmap.png       group x context heatmap

Filtering policy (drop reasons are counted in the sample-stats table):
  - FILTER scope: only records passing the VCF filters are counted by
    default (``--filter PASS``; ``PASS`` and ``.`` both count as passing).
    Somatic spectrum claims must rest on the same high-confidence callset
    as burden statistics; ``--filter all`` keeps every record when the
    caller's unfiltered output is genuinely wanted.
  - non-SNV records (REF or any ALT allele not a single A/C/G/T base),
  - records on contigs absent from the reference FASTA,
  - contexts whose flanks contain non-ACGT bases (N stretches, edges).

Usage:
    python spectrum.py \
        --vcf-map experiment_vcf_map.tsv \
        --group-map experiment_group_map.tsv \
        --ref-fasta genome.fa \
        --output-prefix somatic_spectrum
"""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, TextIO

# setup_logger copied verbatim from src/common/util/LogUtil.py (module bin
# scripts are self-contained; see modules.md 脚本撰写规范 #0).
def setup_logger(
    logger_name: str,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    stream: Optional[TextIO] = None,
) -> logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


logger = setup_logger("Spectrum", stream=sys.stdout)

# Canonical order of the six pyrimidine-strand base substitutions.
BASE_MUT_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical_trinuc_list() -> List[str]:
    """Return the canonical SBS96 channel labels in COSMIC order.

    Channels are grouped by the six pyrimidine-strand substitutions and
    ordered by 5' then 3' flank (A, C, G, T), matching the row order of
    COSMIC signature matrices.
    """
    return [
        f"{up}[{base_mut}]{down}"
        for base_mut in BASE_MUT_ORDER
        for up in "ACGT"
        for down in "ACGT"
    ]


def load_reference_genome(fasta_file: str) -> Dict[str, str]:
    """Load a reference genome into a dict of contig name -> uppercase sequence."""
    from Bio import SeqIO

    logger.info("Loading reference genome...")
    ref_seq = {record.id: str(record.seq).upper() for record in SeqIO.parse(fasta_file, "fasta")}
    logger.info(f"Reference genome loaded: {len(ref_seq)} contigs.")
    return ref_seq


def get_trinucleotide(ref_seq: Dict[str, str], chrom: str, pos: int, ref: str, alt: str) -> Optional[str]:
    """Return the trinucleotide context of a SNV, normalized to the pyrimidine strand.

    The context is formatted as ``<5'base>[<ref>><alt>]<3'base>``; purine
    reference alleles are reverse-complemented so all contexts are reported
    on the pyrimidine (C/T) strand. Returns None when the contig is unknown,
    the position is at a contig edge, or a flank base is not A/C/G/T.
    """
    seq = ref_seq.get(chrom)
    if not seq:
        return None
    pos0 = pos - 1
    if pos0 < 1 or pos0 + 1 >= len(seq):
        return None
    tri = seq[pos0 - 1:pos0 + 2]
    if any(base not in "ACGT" for base in tri):
        return None
    ref_base = ref.upper()
    alt_base = alt.upper()
    # Normalize to the pyrimidine strand.
    if ref_base in ("A", "G"):
        tri = tri.translate(_COMPLEMENT)[::-1]
        ref_base = ref_base.translate(_COMPLEMENT)
        alt_base = alt_base.translate(_COMPLEMENT)
    return f"{tri[0]}[{ref_base}>{alt_base}]{tri[2]}"


def process_vcf(vcf_file: str, ref_seq: Dict[str, str], sample_id: str,
                pass_only: bool = True) -> Dict[str, object]:
    """Count SBS96 contexts of one sample's VCF and account for dropped records.

    Multi-allelic ALT fields (``G,T``) are split and each SNV allele counted
    separately. With ``pass_only``, records whose FILTER is neither ``PASS``
    nor ``.`` are dropped before counting. Records are further dropped as
    non-SNV (REF/ALT allele is not a single A/C/G/T base), unknown contig,
    or context-unavailable (N flank / edge).

    Returns
    -------
    dict with keys: counts (Counter of context label), n_records,
    n_filter_fail, n_snv_alleles, n_counted, n_non_snv, n_unknown_contig,
    n_no_context.
    """
    import pandas as pd

    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
    df = pd.read_csv(vcf_file, sep="\t", comment="#", header=None, usecols=range(7), names=cols)
    n_records = len(df)
    n_filter_fail = 0
    if pass_only:
        keep = df["FILTER"].isin(("PASS", "."))
        n_filter_fail = int((~keep).sum())
        df = df[keep]
    n_snv_alleles = 0
    n_non_snv = 0
    n_unknown_contig = 0
    n_no_context = 0
    counts: Counter = Counter()
    for row in df.itertuples(index=False):
        ref = str(row.REF).upper()
        for alt in str(row.ALT).upper().split(","):
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                n_non_snv += 1
                continue
            n_snv_alleles += 1
            if row.CHROM not in ref_seq:
                n_unknown_contig += 1
                continue
            tri = get_trinucleotide(ref_seq, row.CHROM, int(row.POS), ref, alt)
            if tri is None:
                n_no_context += 1
                continue
            counts[tri] += 1
    n_counted = sum(counts.values())
    logger.info(
        f"{sample_id}: {n_records} records -> {n_snv_alleles} SNV alleles, "
        f"{n_counted} counted (dropped: {n_filter_fail} FILTER, "
        f"{n_non_snv} non-SNV, {n_unknown_contig} unknown contig, "
        f"{n_no_context} no context)"
    )
    return {
        "counts": counts,
        "n_records": n_records,
        "n_filter_fail": n_filter_fail,
        "n_snv_alleles": n_snv_alleles,
        "n_counted": n_counted,
        "n_non_snv": n_non_snv,
        "n_unknown_contig": n_unknown_contig,
        "n_no_context": n_no_context,
    }


def build_per_sample_matrix(all_counts: Dict[str, Counter], sample_order: List[str]):
    """Build the canonical 96 x samples count matrix (all 96 rows always present)."""
    import pandas as pd

    tri_list = canonical_trinuc_list()
    mat = pd.DataFrame(0, index=tri_list, columns=sample_order, dtype=int)
    for sample in sample_order:
        for tri, num in all_counts[sample].items():
            mat.loc[tri, sample] = num
    return mat


def build_group_dataframe(per_sample_mat, sample_groups: Dict[str, str], sample_order: List[str]):
    """Aggregate the per-sample matrix into a context-by-group DataFrame."""
    import pandas as pd

    group_order: List[str] = []
    for sample in sample_order:
        group = sample_groups[sample]
        if group not in group_order:
            group_order.append(group)
    tri_df = pd.DataFrame(0, index=per_sample_mat.index, columns=group_order, dtype=int)
    for sample in sample_order:
        tri_df[sample_groups[sample]] += per_sample_mat[sample]
    return tri_df


def plot_stacked_bar(tri_df, output_file: str) -> None:
    """Render the per-group trinucleotide mutation spectrum as a stacked bar plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tri_df.plot(kind="bar", stacked=True, figsize=(20, 6))
    plt.ylabel("Mutation count")
    plt.xlabel("Trinucleotide context")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_heatmap(tri_df, output_file: str) -> None:
    """Render the context-by-group counts as a heatmap."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.figure(figsize=(12, 8))
    sns.heatmap(tri_df, cmap="Reds", annot=False)
    plt.ylabel("Trinucleotide context")
    plt.xlabel("Group")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def write_sample_stats(stats: Dict[str, Dict[str, object]], sample_groups: Dict[str, str],
                       sample_order: List[str], output_file: str) -> None:
    """Write the per-sample SNV accounting table."""
    fields = [
        "sample_id", "group", "n_records", "n_filter_fail", "n_snv_alleles",
        "n_counted", "n_non_snv", "n_unknown_contig", "n_no_context",
    ]
    with open(output_file, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for sample in sample_order:
            s = stats[sample]
            writer.writerow([
                sample, sample_groups.get(sample, ""),
                s["n_records"], s["n_filter_fail"], s["n_snv_alleles"],
                s["n_counted"], s["n_non_snv"], s["n_unknown_contig"],
                s["n_no_context"],
            ])


def mutation_spectrum_analysis(vcf_files: Dict[str, str], sample_groups: Dict[str, str],
                               ref_fasta: str, output_prefix: str,
                               pass_only: bool = True) -> None:
    """Run the full spectrum analysis and write per-sample + per-group outputs."""
    ref_seq = load_reference_genome(ref_fasta)
    sample_order = list(vcf_files.keys())
    all_counts: Dict[str, Counter] = {}
    stats: Dict[str, Dict[str, object]] = {}
    for sample in sample_order:
        logger.info(f"Processing {sample} ...")
        res = process_vcf(vcf_files[sample], ref_seq, sample, pass_only=pass_only)
        all_counts[sample] = res["counts"]
        stats[sample] = res

    per_sample_mat = build_per_sample_matrix(all_counts, sample_order)
    per_sample_mat.to_csv(f"{output_prefix}_per_sample.csv")
    write_sample_stats(stats, sample_groups, sample_order, f"{output_prefix}_sample_stats.csv")

    tri_df = build_group_dataframe(per_sample_mat, sample_groups, sample_order)
    tri_df.to_csv(f"{output_prefix}.csv")
    plot_stacked_bar(tri_df, f"{output_prefix}_stacked_bar.png")
    plot_heatmap(tri_df, f"{output_prefix}_heatmap.png")
    logger.info(
        f"Wrote {output_prefix}_per_sample.csv, {output_prefix}_sample_stats.csv, "
        f"{output_prefix}.csv, {output_prefix}_stacked_bar.png, {output_prefix}_heatmap.png"
    )


def read_mapping_tsv(file_path: str, key_col: str, value_col: str) -> Dict[str, str]:
    """Read a two-column TSV mapping (header required) into a dict."""
    mapping: Dict[str, str] = {}
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
    parser.add_argument("-v", "--vcf-map", required=True,
                        help="TSV with columns: experiment_sample_id, vcf")
    parser.add_argument("-g", "--group-map", required=True,
                        help="TSV with columns: experiment_sample_id, group")
    parser.add_argument("-r", "--ref-fasta", required=True,
                        help="Reference FASTA")
    parser.add_argument("-o", "--output-prefix", required=True,
                        help="Output prefix")
    parser.add_argument("--filter", choices=["PASS", "all"], default="PASS",
                        help="FILTER scope: PASS counts only records passing "
                             "the VCF filters (PASS/'.'); all keeps every "
                             "record (default: PASS)")
    parser.add_argument("-l", "--log", default=None,
                        help="Optional log file (in addition to stdout)")
    return parser.parse_args()


def main():
    """Entry point: read maps, run the analysis, write outputs."""
    args = parse_args()
    if args.log:
        _log_path = Path(args.log)
        _log_path.parent.mkdir(parents=True, exist_ok=True)
        _lg = setup_logger("Spectrum", log_file=str(_log_path), stream=sys.stdout)
        logger.handlers = _lg.handlers
    vcf_files = read_mapping_tsv(args.vcf_map, "experiment_sample_id", "vcf")
    sample_groups = read_mapping_tsv(args.group_map, "experiment_sample_id", "group")
    missing_groups = [s for s in vcf_files if s not in sample_groups]
    if missing_groups:
        raise ValueError(
            f"Samples missing from the group map: {missing_groups}; "
            f"every sample in {args.vcf_map} needs a row in {args.group_map}."
        )
    mutation_spectrum_analysis(vcf_files, sample_groups, args.ref_fasta,
                               args.output_prefix, pass_only=(args.filter == "PASS"))


if __name__ == "__main__":
    main()
