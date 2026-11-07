#!/usr/bin/env python
"""Mutation-burden calculations from somatic variant-caller outputs.

Three burden metrics from three caller outputs, plus a summary merger:

* ``snv-indel``: tumor mutational burden (TMB) and indel burden from a
  filtered somatic SNV/InDel VCF (e.g. GATK Mutect2 + FilterMutectCalls).
* ``sv``: structural-variant burden from a somatic SV VCF (e.g. Manta
  ``somaticSV.vcf.gz``).
* ``cnv``: copy-number burden from a CNVkit ``.cns`` segment table.
* ``aggregate``: merge per-sample burden tables into one summary TSV.

The ``snv-indel`` step annotates the input VCF with Ensembl VEP (via
``bin/vep_ann.py``, cache-based offline) and counts only meaningful —
protein-impacting or core-splice — variants (VEP consequence whitelist).
This annotation + filtering is switchable and ON by default: ``--no-annotate``
skips VEP and counts every passing allele without consequence filtering.

Denominator for the per-megabase rates (``snv-indel`` only), resolved in
precedence order: ``--targets-bed`` (captured territory) > ``--gtf``
(merged exon/CDS length) > ``--size`` (fixed coding size in Mb). Variant
counting and territory computation share the same chromosome universe
(canonical human/mouse autosomes, optionally plus X/Y).

Usage:
    python burden.py snv-indel -i pair.filtered.vcf.gz -s TUMOR \\
        -o TUMOR.snv_indel_burden.tsv -b targets.bed -v 0.05
    python burden.py snv-indel -i pair.filtered.vcf.gz -s TUMOR \\
        -o TUMOR.raw_burden.tsv --no-annotate
    python burden.py sv -i somaticSV.vcf.gz -s TUMOR -o TUMOR.sv_burden.tsv
    python burden.py cnv -i TUMOR.cns -s TUMOR -o TUMOR.cnv_burden.tsv
    python burden.py aggregate --snv *.snv_indel_burden.tsv \\
        --sv *.sv_burden.tsv --cnv *.cnv_burden.tsv -o burden_summary.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, NamedTuple, Optional, Sequence

import pandas as pd

from LogUtil import setup_logger
from vep_ann import VEP_Annotation

# Protein-impacting / core-splice consequences counted toward TMB when the
# VCF carries VEP CSQ annotations (Ensembl VEP predicted consequences:
# https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html).
INCLUDE_CONSEQUENCES = {
    "transcript_ablation",
    "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "feature_elongation", "feature_truncation",
    "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant",
}

# Per-type output columns. Names are already unique across types so the
# aggregate step can outer-join every table on "sample" without renaming.
SNV_INDEL_HEADER = [
    "sample", "snv_count", "indel_count", "tmb_count", "frameshift_count",
    "tmb_per_mb", "snv_per_mb", "indel_per_mb", "territory_mb", "denominator",
]
SV_HEADER = ["sample", "sv_total", "sv_DEL", "sv_DUP", "sv_INV", "sv_INS",
             "sv_BND", "sv_OTHER"]
CNV_HEADER = [
    "sample", "cnv_total_segments", "cnv_gain_segments", "cnv_loss_segments",
    "cnv_altered_length_mb", "cnv_altered_gene_count",
]
SV_TYPES = ["DEL", "DUP", "INV", "INS", "BND"]

# Single source of truth for CLI defaults and programmatic calls alike.
DEFAULT_MIN_VAF = 0.05
DEFAULT_SIZE_MB = 34.0
DEFAULT_GTF_FEATURE = "exon"
DEFAULT_MIN_ABS_LOG2 = 0.3
DEFAULT_VEP_CACHE = "~/.vep"
DEFAULT_VEP_SPECIES = "mus_musculus"
DEFAULT_VEP_ASSEMBLY = "GRCm39"


# Module-level logger; reconfigured in main() with the CLI log file (same
# singleton object, stream on stdout).
logger = setup_logger("Burden", stream=sys.stdout)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _open_text(path: str | Path) -> IO[str]:
    """Open a plain or gzip-compressed text file for reading.

    Parameters
    ----------
    path : str | Path
        File path ending in ``.gz`` (gzip) or not (plain text).

    Returns
    -------
    IO[str]
        Text-mode file handle.
    """
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def _chrom_universe(include_sex: bool) -> set[str]:
    """Build the set of chromosome names eligible for counting.

    Canonical human (1-22) and mouse (1-19) autosomes, with and without the
    ``chr`` prefix; X/Y are only included when ``include_sex`` is set.
    Non-canonical scaffolds are always excluded.

    Parameters
    ----------
    include_sex : bool
        Include X/Y alongside the autosomes.

    Returns
    -------
    set[str]
        Eligible chromosome names.
    """
    chroms = {str(i) for i in range(1, 23)} | {f"chr{i}" for i in range(1, 23)}
    chroms |= {str(i) for i in range(1, 20)} | {f"chr{i}" for i in range(1, 20)}
    if include_sex:
        chroms |= {"X", "Y", "chrX", "chrY"}
    return chroms


def _merged_length(intervals: Sequence[tuple[int, int]]) -> int:
    """Merge 1-based inclusive ``(start, end)`` intervals and sum their length.

    Parameters
    ----------
    intervals : Sequence[tuple[int, int]]
        Intervals as 1-based inclusive coordinate pairs.

    Returns
    -------
    int
        Total length in bp of the union of the intervals (0 when empty).
    """
    if not intervals:
        return 0
    intervals = sorted(intervals, key=lambda x: x[0])
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    total += cur_end - cur_start + 1
    return total


def _bed_territory_bp(bed_file: str | Path, universe: set[str]) -> int:
    """Sum the merged interval length of a BED target file.

    Parameters
    ----------
    bed_file : str | Path
        BED file (0-based half-open coordinates; ``track``/``browser``
        header lines ignored).

    universe : set[str]
        Chromosome names eligible for counting (see `_chrom_universe`).

    Returns
    -------
    int
        Merged territory length in bp, restricted to ``universe``.
    """
    per_chrom: dict[str, list[tuple[int, int]]] = {}
    with open(bed_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            if chrom not in universe:
                continue
            start = int(fields[1]) + 1  # BED is 0-based half-open -> 1-based
            end = int(fields[2])
            if end >= start:
                per_chrom.setdefault(chrom, []).append((start, end))
    return sum(_merged_length(v) for v in per_chrom.values())


def _gtf_territory_bp(gtf_file: str | Path, universe: set[str],
                      feature: str) -> int:
    """Sum the merged length of one GTF feature type (``exon`` or ``CDS``).

    Parameters
    ----------
    gtf_file : str | Path
        GTF annotation file (comment lines ignored).

    universe : set[str]
        Chromosome names eligible for counting (see `_chrom_universe`).

    feature : str
        GTF feature type whose intervals are merged, e.g. ``exon`` or
        ``CDS``.

    Returns
    -------
    int
        Merged feature length in bp, restricted to ``universe``.
    """
    per_chrom: dict[str, list[tuple[int, int]]] = {}
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[2] != feature:
                continue
            chrom = fields[0]
            if chrom not in universe:
                continue
            per_chrom.setdefault(chrom, []).append((int(fields[3]), int(fields[4])))
    return sum(_merged_length(v) for v in per_chrom.values())


class Territory(NamedTuple):
    """Denominator for per-megabase burden rates.

    Attributes
    ----------
    mb : float
        Territory size in megabases.
    denominator : str
        Provenance label: ``targets_bed``, ``gtf_<feature>`` or
        ``fixed_size``.
    """

    mb: float
    denominator: str


def resolve_territory(size_mb: float,
                      targets_bed: Optional[str] = None,
                      gtf: Optional[str] = None,
                      gtf_feature: str = DEFAULT_GTF_FEATURE,
                      universe: Optional[set[str]] = None) -> Territory:
    """Resolve the per-Mb denominator in precedence order.

    Precedence: ``targets_bed`` (captured territory) > ``gtf`` (merged
    feature length) > ``size_mb`` (fixed coding size). A configured source
    that yields 0 bp on the counting chromosomes falls back to ``size_mb``
    with a warning.

    Parameters
    ----------
    size_mb : float
        Fixed territory in Mb used when neither BED nor GTF is given (or
        both yield 0 bp).
    targets_bed : Optional[str], optional
        Capture targets BED path (highest precedence).
    gtf : Optional[str], optional
        Gene annotation GTF path (second precedence).
    gtf_feature : str, optional
        GTF feature used for the GTF denominator, by default ``exon``.
    universe : Optional[set[str]], optional
        Chromosome universe; derived with ``include_sex=False`` when None.

    Returns
    -------
    Territory
        Resolved size (Mb) plus its provenance label.
    """
    if universe is None:
        universe = _chrom_universe(include_sex=False)
    territory_bp = 0
    denominator = "fixed_size"
    if targets_bed:
        _require_input(targets_bed, "targets BED")
        territory_bp = _bed_territory_bp(targets_bed, universe)
        denominator = "targets_bed"
    elif gtf:
        _require_input(gtf, "gene annotation GTF")
        territory_bp = _gtf_territory_bp(gtf, universe, gtf_feature)
        denominator = f"gtf_{gtf_feature}"
    if territory_bp > 0:
        return Territory(territory_bp / 1e6, denominator)
    if targets_bed or gtf:
        logger.warning(
            "Configured territory source yielded 0 bp on the counting "
            "chromosomes; falling back to --size.")
    return Territory(float(size_mb), "fixed_size")


def _write_tsv(path: str | Path, header: Sequence[str],
               rows: Sequence[Sequence[object]]) -> None:
    """Atomically write a TSV (header + rows); always emits every column.

    Parameters
    ----------
    path : str | Path
        Output TSV path (parents created; written via a ``.tmp`` file and
        ``os.replace``).
    header : Sequence[str]
        Column names.
    rows : Sequence[Sequence[object]]
        Data rows, one sequence of cells per row.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    os.replace(tmp, path)
    logger.info(f"Wrote {path} ({len(rows)} data row(s))")


def _require_input(path: Optional[str], what: str) -> None:
    """Fail with a clear message when a required input file is missing.

    Parameters
    ----------
    path : Optional[str]
        Input path as given on the CLI / config.
    what : str
        Human-readable description of the input, used in the error message.

    Raises
    ------
    FileNotFoundError
        When ``path`` is empty or does not exist.
    """
    if not path or not os.path.exists(path):
        raise FileNotFoundError(
            f"Missing {what} input: {path!r} (set the module config key that "
            f"provides it, then rerun)."
        )
    logger.debug(f"{what} input OK: {path}")


# ---------------------------------------------------------------------------
# snv-indel: TMB / indel burden
# ---------------------------------------------------------------------------

def _parse_csq(info_field: str) -> list[str]:
    """Return the list of VEP CSQ entries (allele|Consequence|...) from INFO.

    Parameters
    ----------
    info_field : str
        VCF INFO column value.

    Returns
    -------
    list[str]
        CSQ entry strings; empty when the record carries no CSQ.
    """
    for item in info_field.split(";"):
        if item.startswith("CSQ="):
            return item.split("=", 1)[1].split(",")
    return []


def _allele_matches(csq_allele: str, alt: str, ref: str) -> bool:
    """Match a CSQ allele token against a VCF REF/ALT pair.

    VEP encodes deletions as ``-``; treat them as matching shorter ALT
    alleles.

    Parameters
    ----------
    csq_allele : str
        Allele token taken from the first CSQ field.
    alt : str
        VCF ALT allele.
    ref : str
        VCF REF allele.

    Returns
    -------
    bool
        True when the CSQ entry refers to this ALT allele.
    """
    if csq_allele == alt:
        return True
    return csq_allele == "-" and len(alt) < len(ref)


def _sample_allele_fraction(format_field: str, sample_field: str,
                            alt_idx: int) -> Optional[float]:
    """Extract one ALT allele's fraction in the tumor sample column.

    Uses the FORMAT ``AF`` field (comma-separated per ALT); falls back to a
    fraction computed from ``AD`` (ref + alt depths) when AF is absent.

    Parameters
    ----------
    format_field : str
        VCF FORMAT column value (``AF:AD:...``).
    sample_field : str
        Tumor sample column value.
    alt_idx : int
        Zero-based index of the ALT allele within the record.

    Returns
    -------
    Optional[float]
        Allele fraction, or None when not computable from AF/AD.
    """
    keys = format_field.split(":")
    parts = sample_field.split(":")

    def _get(key: str) -> Optional[str]:
        if key not in keys:
            return None
        idx = keys.index(key)
        return parts[idx] if idx < len(parts) else None

    af = _get("AF")
    if af is not None:
        values = af.split(",")
        if alt_idx < len(values):
            try:
                return float(values[alt_idx])
            except ValueError:
                pass
    ad = _get("AD")
    if ad is not None:
        try:
            depths = [float(t) for t in ad.split(",")]
        except ValueError:
            return None
        if alt_idx + 1 < len(depths):
            total = depths[0] + depths[alt_idx + 1]
            if total > 0:
                return depths[alt_idx + 1] / total
    return None


@dataclass
class VariantCounts:
    """Allele counts from one somatic SNV/InDel VCF.

    Attributes
    ----------
    snv : int
        SNV alleles passing all filters.
    indel : int
        InDel alleles passing all filters.
    frameshift : int
        Indels annotated as frameshift_variant (annotate mode only).
    skipped_no_af : int
        Alleles skipped because AF/AD yielded no fraction.
    skipped_no_csq : int
        Records skipped for carrying no CSQ annotation (annotate mode).
    csq_seen : bool
        True once any record carried a CSQ field.
    """

    snv: int = 0
    indel: int = 0
    frameshift: int = 0
    skipped_no_af: int = 0
    skipped_no_csq: int = 0
    csq_seen: bool = field(default=False)

    @property
    def tmb(self) -> int:
        """int: Total counted alleles (SNV + InDel)."""
        return self.snv + self.indel


def count_snv_indel(vcf: str | Path,
                    sample: str,
                    min_vaf: float = DEFAULT_MIN_VAF,
                    annotate: bool = True,
                    include_filtered: bool = False,
                    universe: Optional[set[str]] = None) -> VariantCounts:
    """Count SNVs, indels and frameshifts in one somatic VCF.

    In annotate mode only alleles whose VEP CSQ consequences intersect the
    protein-impacting / core-splice whitelist are counted, and records
    without any CSQ are skipped. Without annotation every allele passing
    the record FILTER (PASS/``.`` unless ``include_filtered``) is counted.

    Parameters
    ----------
    vcf : str | Path
        Filtered somatic SNV/InDel VCF (plain or gzipped), already carrying
        VEP CSQ annotations when ``annotate`` is True.
    sample : str
        Tumor sample name; must match a column in the ``#CHROM`` header.
    min_vaf : float, optional
        Minimum tumor allele fraction, by default 0.05.
    annotate : bool, optional
        Require VEP CSQ annotations and filter by consequence whitelist,
        by default True.
    include_filtered : bool, optional
        Count non-PASS records too, by default False (PASS only).
    universe : Optional[set[str]], optional
        Chromosome universe; autosomes-only when None.

    Returns
    -------
    VariantCounts
        Allele counts plus skip statistics.

    Raises
    ------
    ValueError
        When the tumor sample column is absent, the file has no
        ``#CHROM`` header, or annotate mode saw no CSQ anywhere.
    """
    if universe is None:
        universe = _chrom_universe(include_sex=False)
    counts = VariantCounts()
    with _open_text(vcf) as fh:
        tumor_idx: Optional[int] = None
        for line in fh:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                columns = fields
                if sample not in columns:
                    raise ValueError(
                        f"Tumor sample column {sample!r} not found in "
                        f"{vcf}; available columns: {columns[9:]}")
                tumor_idx = columns.index(sample)
                continue
            if tumor_idx is None:
                raise ValueError(f"{vcf} has no #CHROM header line")
            chrom, _pos, _id, ref, alt_field = (fields[0], fields[1],
                                                fields[2], fields[3],
                                                fields[4])
            if chrom not in universe:
                continue
            flt = fields[6]
            if not include_filtered and flt not in (".", "PASS"):
                continue
            info = fields[7]
            fmt = fields[8] if len(fields) > 8 else ""
            tumor = fields[tumor_idx] if tumor_idx < len(fields) else ""

            csq_entries = _parse_csq(info)
            has_csq = bool(csq_entries)
            counts.csq_seen = counts.csq_seen or has_csq
            if annotate and not has_csq:
                counts.skipped_no_csq += 1
                continue

            for alt_idx, alt in enumerate(alt_field.split(",")):
                if alt in (".", "*", "") or alt.startswith("<"):
                    continue
                if has_csq and annotate:
                    allele_csq = [
                        entry for entry in csq_entries
                        if _allele_matches(entry.split("|", 1)[0], alt, ref)
                    ]
                    if not allele_csq:
                        continue
                    tokens = set()
                    for entry in allele_csq:
                        fields_csq = entry.split("|")
                        if len(fields_csq) > 1 and fields_csq[1]:
                            tokens.update(fields_csq[1].split("&"))
                    if not (tokens & INCLUDE_CONSEQUENCES):
                        continue
                    is_frameshift = "frameshift_variant" in tokens
                else:
                    is_frameshift = False
                af = _sample_allele_fraction(fmt, tumor, alt_idx)
                if af is None:
                    counts.skipped_no_af += 1
                    continue
                if af < min_vaf:
                    continue
                if len(ref) == 1 and len(alt) == 1:
                    counts.snv += 1
                else:
                    counts.indel += 1
                    if is_frameshift:
                        counts.frameshift += 1

    if annotate and not counts.csq_seen:
        raise ValueError(
            f"VEP annotation produced no CSQ field in {vcf}; the annotated "
            f"VCF is unusable for burden counting. Check the VEP cache/"
            f"species/assembly settings, or rerun with --no-annotate.")
    return counts


def run_snv_indel(vcf: str,
                  sample: str,
                  out: str,
                  min_vaf: float = DEFAULT_MIN_VAF,
                  size_mb: float = DEFAULT_SIZE_MB,
                  targets_bed: Optional[str] = None,
                  gtf: Optional[str] = None,
                  gtf_feature: str = DEFAULT_GTF_FEATURE,
                  annotate: bool = True,
                  vep_cache: str = DEFAULT_VEP_CACHE,
                  vep_species: str = DEFAULT_VEP_SPECIES,
                  vep_assembly: str = DEFAULT_VEP_ASSEMBLY,
                  include_sex: bool = False,
                  include_filtered: bool = False) -> None:
    """Compute SNV/InDel burden and TMB for one filtered somatic VCF.

    With ``annotate`` (default) the VCF is first annotated with Ensembl VEP
    (``vep_ann.VEP_Annotation``; the annotated copy is written beside the
    output as ``<sample>.vep.vcf``) and only variants carrying VEP CSQ
    consequences from the protein-impacting / core-splice whitelist are
    counted — annotation is mandatory in this mode and an unannotatable VCF
    raises. With ``annotate=False`` every allele passing the record FILTER
    (PASS/``.`` unless ``include_filtered``) is counted without consequence
    filtering. Both modes apply the VAF floor and the chromosome universe;
    SNVs and indels split by allele length.

    Writes one row even for zero-count inputs (header + zero values).

    Parameters
    ----------
    vcf : str
        Filtered somatic SNV/InDel VCF path (``.vcf``/``.vcf.gz``).
    sample : str
        Sample label; also the tumor sample-column name in the VCF.
    out : str
        Output TSV path.
    min_vaf : float, optional
        Minimum tumor allele fraction, by default 0.05.
    size_mb : float, optional
        Fixed coding territory in Mb when no BED/GTF is given, by default
        34.0.
    targets_bed : Optional[str], optional
        Capture targets BED; summed territory becomes the per-Mb
        denominator (highest precedence).
    gtf : Optional[str], optional
        Gene annotation GTF; merged feature length becomes the per-Mb
        denominator.
    gtf_feature : str, optional
        GTF feature used for the GTF denominator, by default ``exon``.
    annotate : bool, optional
        Annotate with VEP and count only protein-impacting / core-splice
        variants, by default True.
    vep_cache : str, optional
        VEP cache root directory, by default ``~/.vep``.
    vep_species : str, optional
        VEP species name, by default ``mus_musculus``.
    vep_assembly : str, optional
        VEP assembly name, by default ``GRCm39``.
    include_sex : bool, optional
        Count X/Y in addition to autosomes, by default False.
    include_filtered : bool, optional
        Count non-PASS records too, by default False.

    Raises
    ------
    FileNotFoundError
        When the VCF (or a configured BED/GTF) is missing.
    ValueError
        When the VCF header lacks the tumor sample column, or annotate
        mode produced no CSQ anywhere.
    """
    _require_input(vcf, "somatic SNV/InDel VCF")
    universe = _chrom_universe(include_sex)

    # Annotation + meaningful-variant filtering (default on): annotate the
    # filtered VCF with VEP first, then count from the annotated copy.
    source_vcf = vcf
    if annotate:
        vep_out = str(Path(out).parent / f"{sample}.vep.vcf")
        Path(vep_out).parent.mkdir(parents=True, exist_ok=True)
        analysis = VEP_Annotation(
            vep_cache_dir=vep_cache,
            species=vep_species,
            assembly=vep_assembly,
        )
        analysis.annotate_sv_vep(vcf, vep_out)
        source_vcf = vep_out

    # Denominator: captured territory > GTF feature length > fixed size.
    territory = resolve_territory(size_mb=size_mb, targets_bed=targets_bed,
                                   gtf=gtf, gtf_feature=gtf_feature,
                                   universe=universe)

    logger.info(
        f"snv-indel: sample={sample} vcf={source_vcf} annotate={annotate} "
        f"min_vaf={min_vaf} territory_mb={territory.mb:.4f} "
        f"denominator={territory.denominator}")

    counts = count_snv_indel(vcf=source_vcf, sample=sample, min_vaf=min_vaf,
                             annotate=annotate,
                             include_filtered=include_filtered,
                             universe=universe)

    if counts.skipped_no_csq:
        logger.warning(
            f"{counts.skipped_no_csq} record(s) skipped: no CSQ annotation.")
    if not annotate:
        logger.info(
            "Annotation filtering disabled (--no-annotate): counting all "
            "passing alleles without consequence filtering.")
    if counts.skipped_no_af:
        logger.warning(
            f"{counts.skipped_no_af} allele(s) skipped: allele fraction not "
            f"computable from FORMAT AF/AD.")

    _write_tsv(out, SNV_INDEL_HEADER, [[
        sample, counts.snv, counts.indel, counts.tmb, counts.frameshift,
        f"{counts.tmb / territory.mb:.4f}",
        f"{counts.snv / territory.mb:.4f}",
        f"{counts.indel / territory.mb:.4f}",
        f"{territory.mb:.4f}", territory.denominator,
    ]])


# ---------------------------------------------------------------------------
# sv: structural-variant burden
# ---------------------------------------------------------------------------

def count_sv(vcf: str | Path) -> dict[str, int]:
    """Count somatic structural variants by INFO/SVTYPE bucket.

    Parameters
    ----------
    vcf : str | Path
        Somatic SV VCF (plain or gzipped), Manta ``somaticSV`` convention.

    Returns
    -------
    dict[str, int]
        Counts keyed by ``total``, each of DEL/DUP/INV/INS/BND, and
        ``other``.
    """
    counts = {sv_type: 0 for sv_type in SV_TYPES}
    counts["other"] = 0
    counts["total"] = 0
    with _open_text(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            counts["total"] += 1
            sv_type = None
            for item in fields[7].split(";"):
                if item.startswith("SVTYPE="):
                    sv_type = item.split("=", 1)[1]
                    break
            if sv_type in SV_TYPES:
                counts[sv_type] += 1
            else:
                counts["other"] += 1
    return counts


def run_sv(vcf: str, sample: str, out: str) -> None:
    """Count somatic structural variants by SVTYPE for one SV VCF.

    Reads INFO/SVTYPE per record and buckets DEL/DUP/INV/INS/BND;
    everything else lands in ``sv_OTHER``. Writes one row even for
    zero-count inputs.

    Parameters
    ----------
    vcf : str
        Somatic SV VCF path (``.vcf``/``.vcf.gz``), e.g. Manta
        ``somaticSV.vcf.gz``.
    sample : str
        Sample label.
    out : str
        Output TSV path.

    Raises
    ------
    FileNotFoundError
        When the SV VCF is missing.
    """
    _require_input(vcf, "somatic SV VCF")
    counts = count_sv(vcf)
    logger.info(f"sv: sample={sample} vcf={vcf} total={counts['total']}")
    _write_tsv(out, SV_HEADER, [[
        sample, counts["total"], counts["DEL"], counts["DUP"],
        counts["INV"], counts["INS"], counts["BND"], counts["other"],
    ]])


# ---------------------------------------------------------------------------
# cnv: copy-number burden
# ---------------------------------------------------------------------------

@dataclass
class CnvSummary:
    """Altered-segment summary of one CNVkit ``.cns`` table.

    Attributes
    ----------
    total : int
        Total segment rows read.
    gain : int
        Altered segments called as gain.
    loss : int
        Altered segments called as loss.
    altered_bp : int
        Total length (bp) of altered segments with parseable coordinates.
    genes : set[str]
        Unique gene symbols inside altered segments.
    """

    total: int = 0
    gain: int = 0
    loss: int = 0
    altered_bp: int = 0
    genes: set[str] = field(default_factory=set)


def summarize_cnv(cns: str | Path,
                  min_abs_log2: float = DEFAULT_MIN_ABS_LOG2) -> CnvSummary:
    """Summarize altered copy-number segments from a CNVkit ``.cns`` table.

    A segment counts as altered when the ``call`` column is present and
    non-zero, otherwise when ``abs(log2) >= min_abs_log2``.

    Parameters
    ----------
    cns : str | Path
        CNVkit ``.cns`` segment table (TSV with header).
    min_abs_log2 : float, optional
        Altered threshold on ``|log2|`` when the table has no usable
        ``call`` column, by default 0.3.

    Returns
    -------
    CnvSummary
        Segment counts, altered length and altered-gene set.

    Raises
    ------
    ValueError
        When the table lacks the required columns (start/end plus log2 or
        call).
    """
    summary = CnvSummary()
    with open(cns) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = reader.fieldnames or []
        colmap = {c.lower(): c for c in columns}

        def _first(*names: str) -> Optional[str]:
            for name in names:
                if name in colmap:
                    return colmap[name]
            return None

        start_col = _first("start")
        end_col = _first("end")
        log2_col = _first("log2")
        call_col = _first("call")
        gene_col = _first("gene", "genes", "gene_names", "name")
        if not (start_col and end_col and (log2_col or call_col)):
            raise ValueError(
                f"{cns} lacks required columns (start/end/log2 or call); "
                f"found: {columns}")
        for row in reader:
            summary.total += 1
            if call_col and row.get(call_col, "").strip() not in ("", "."):
                call = row[call_col].strip()
                altered = call not in ("0", "neutral")
                is_gain = call in ("+", "1", "+1", "gain")
                is_loss = call in ("-", "-1", "loss")
            else:
                try:
                    log2 = float(row[log2_col])
                except (KeyError, TypeError, ValueError):
                    continue
                altered = abs(log2) >= min_abs_log2
                is_gain = log2 >= min_abs_log2
                is_loss = log2 <= -min_abs_log2
            if not altered:
                continue
            if is_gain:
                summary.gain += 1
            elif is_loss:
                summary.loss += 1
            try:
                summary.altered_bp += int(row[end_col]) - int(row[start_col])
            except (KeyError, TypeError, ValueError):
                pass
            if gene_col:
                for gene in row.get(gene_col, "").split(","):
                    gene = gene.strip()
                    if gene and gene not in (".", "-"):
                        summary.genes.add(gene)
    return summary


def run_cnv(cns: str, sample: str, out: str,
            min_abs_log2: float = DEFAULT_MIN_ABS_LOG2) -> None:
    """Summarize altered copy-number segments for one CNVkit ``.cns`` table.

    Reports gain/loss segment counts, total altered length (Mb) and the
    number of unique gene symbols inside altered segments. Writes one row
    even for zero altered segments.

    Parameters
    ----------
    cns : str
        CNVkit ``.cns`` segment table path.
    sample : str
        Sample label.
    out : str
        Output TSV path.
    min_abs_log2 : float, optional
        Altered threshold on ``|log2|`` when the table has no usable
        ``call`` column, by default 0.3.

    Raises
    ------
    FileNotFoundError
        When the ``.cns`` table is missing.
    ValueError
        When the table lacks the required columns.
    """
    _require_input(cns, "CNVkit .cns table")
    summary = summarize_cnv(cns, min_abs_log2=min_abs_log2)
    logger.info(
        f"cnv: sample={sample} cns={cns} total_segments={summary.total} "
        f"altered={summary.gain + summary.loss}")
    _write_tsv(out, CNV_HEADER, [[
        sample, summary.total, summary.gain, summary.loss,
        f"{summary.altered_bp / 1e6:.4f}", len(summary.genes),
    ]])


# ---------------------------------------------------------------------------
# aggregate: merge per-sample tables
# ---------------------------------------------------------------------------

def run_aggregate(out: str,
                  snv_tables: Optional[Sequence[str]] = None,
                  sv_tables: Optional[Sequence[str]] = None,
                  cnv_tables: Optional[Sequence[str]] = None) -> None:
    """Merge the per-sample burden tables into one summary TSV.

    Tables of the SAME type are concatenated row-wise (one row per sample,
    identical columns); the type frames are then outer-joined on the
    ``sample`` column — per-type column names are already unique across
    types. A sample present in only some types keeps empty cells for the
    others. With no input tables a header-only TSV is written.

    Parameters
    ----------
    out : str
        Output summary TSV path.
    snv_tables : Optional[Sequence[str]], optional
        Per-sample ``snv_indel`` burden TSV paths.
    sv_tables : Optional[Sequence[str]], optional
        Per-sample ``sv`` burden TSV paths.
    cnv_tables : Optional[Sequence[str]], optional
        Per-sample ``cnv`` burden TSV paths.
    """
    type_frames = []
    for paths, label in ((snv_tables or [], "snv_indel"),
                         (sv_tables or [], "sv"),
                         (cnv_tables or [], "cnv")):
        tables = []
        for path in paths:
            _require_input(path, f"{label} burden table")
            tables.append(pd.read_csv(path, sep="\t"))
        logger.info(f"aggregate: {label}: {len(tables)} table(s)")
        if tables:
            frame = pd.concat(tables, ignore_index=True)
            # One row per sample per type is the contract; keep the last on
            # accidental duplicate submissions instead of exploding the join.
            dup = frame["sample"].duplicated().sum()
            if dup:
                logger.warning(f"{label}: {dup} duplicate sample row(s) dropped")
                frame = frame.drop_duplicates("sample", keep="last")
            type_frames.append(frame)
    if type_frames:
        merged = type_frames[0]
        for frame in type_frames[1:]:
            merged = merged.merge(frame, on="sample", how="outer")
        merged = merged.sort_values("sample").reset_index(drop=True)
    else:
        merged = pd.DataFrame(columns=["sample"])
    _write_tsv(out, list(merged.columns),
               merged.fillna("").astype(str).values.tolist())


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_common(parser: argparse.ArgumentParser) -> None:
    """Attach the options shared by every per-sample subcommand.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Subparser receiving the common ``-s``/``-o``/``-l``/``-w``
        options.
    """
    parser.add_argument("-s", "--sample", required=True,
                        help="sample label (tumor sample; also the VCF "
                             "sample-column name for snv-indel)")
    parser.add_argument("-o", "--out", required=True,
                        help="output TSV path (always written, even when all "
                             "counts are zero)")
    parser.add_argument("-l", "--log", default=None, type=Path,
                        help="optional log file path (default: stdout only)")
    parser.add_argument("-w", "--overwrite", action="store_true",
                        help="rewrite the output when it already exists "
                             "(default: keep it and exit)")


def _run_command(args: argparse.Namespace) -> None:
    """Shared skip-existing gate, dispatch and failure wrap.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments; ``args.command`` selects the runner and
        ``args.func`` receives the matching explicit keyword arguments.
    """
    if os.path.exists(args.out) and not args.overwrite:
        logger.info(f"{args.out} exists; keeping it (use -w/--overwrite to redo)")
        return
    subject = (f" for sample {args.sample}"
               if args.command != "aggregate" else "")
    try:
        args.func(args)
    except Exception as exc:
        logger.error(f"burden {args.command} failed{subject}: {exc}")
        raise RuntimeError(
            f"burden {args.command} failed{subject}: {exc}; "
            f"see the log for details.")
    logger.info(f"burden {args.command} completed{subject}")


def _run_snv_indel_cli(args: argparse.Namespace) -> None:
    """Adapt the snv-indel CLI namespace to `run_snv_indel` parameters."""
    run_snv_indel(vcf=args.vcf, sample=args.sample, out=args.out,
                  min_vaf=args.min_vaf, size_mb=args.size,
                  targets_bed=args.targets_bed, gtf=args.gtf,
                  gtf_feature=args.gtf_feature, annotate=args.annotate,
                  vep_cache=args.vep_cache, vep_species=args.vep_species,
                  vep_assembly=args.vep_assembly, include_sex=args.include_sex,
                  include_filtered=args.include_filtered)


def _run_sv_cli(args: argparse.Namespace) -> None:
    """Adapt the sv CLI namespace to `run_sv` parameters."""
    run_sv(vcf=args.vcf, sample=args.sample, out=args.out)


def _run_cnv_cli(args: argparse.Namespace) -> None:
    """Adapt the cnv CLI namespace to `run_cnv` parameters."""
    run_cnv(cns=args.cns, sample=args.sample, out=args.out,
            min_abs_log2=args.min_abs_log2)


def _run_aggregate_cli(args: argparse.Namespace) -> None:
    """Adapt the aggregate CLI namespace to `run_aggregate` parameters."""
    run_aggregate(out=args.out, snv_tables=args.snv, sv_tables=args.sv,
                  cnv_tables=args.cnv)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments for the burden subcommands.

    Parameters
    ----------
    argv : Optional[Sequence[str]], optional
        Argument list (defaults to ``sys.argv[1:]``).

    Returns
    -------
    argparse.Namespace
        Parsed arguments with a ``func`` dispatcher per subcommand.
    """
    parser = argparse.ArgumentParser(
        description="Mutation-burden calculations (SNV/InDel TMB, SV, CNV).")
    sub = parser.add_subparsers(dest="command", required=True)

    p_snv = sub.add_parser(
        "snv-indel", help="TMB and indel burden from a filtered somatic VCF")
    _add_common(p_snv)
    p_snv.add_argument("-i", "--vcf", required=True,
                       help="filtered somatic SNV/InDel VCF (.vcf/.vcf.gz)")
    p_snv.add_argument("-v", "--min-vaf", type=float, default=DEFAULT_MIN_VAF,
                       help=f"minimum tumor allele fraction (default: "
                            f"{DEFAULT_MIN_VAF})")
    p_snv.add_argument("-m", "--size", type=float, default=DEFAULT_SIZE_MB,
                       help=f"fixed coding territory in Mb when no BED/GTF "
                            f"is given (default: {DEFAULT_SIZE_MB})")
    p_snv.add_argument("-b", "--targets-bed", default=None,
                       help="capture targets BED; summed territory becomes "
                            "the per-Mb denominator (highest precedence)")
    p_snv.add_argument("-g", "--gtf", default=None,
                       help="gene annotation GTF; merged feature length "
                            "becomes the per-Mb denominator")
    p_snv.add_argument("-e", "--gtf-feature", default=DEFAULT_GTF_FEATURE,
                       choices=["exon", "CDS"],
                       help=f"GTF feature used for the GTF denominator "
                            f"(default: {DEFAULT_GTF_FEATURE})")
    p_snv.add_argument("-a", "--annotate", action=argparse.BooleanOptionalAction,
                       default=True,
                       help="annotate with VEP and count only protein-"
                            "impacting / core-splice variants (default: "
                            "enabled; --no-annotate counts every passing "
                            "allele without consequence filtering)")
    p_snv.add_argument("-c", "--vep-cache", default=DEFAULT_VEP_CACHE,
                       help=f"VEP cache root directory (default: "
                            f"{DEFAULT_VEP_CACHE})")
    p_snv.add_argument("-p", "--vep-species", default=DEFAULT_VEP_SPECIES,
                       help=f"VEP species name (default: {DEFAULT_VEP_SPECIES})")
    p_snv.add_argument("-y", "--vep-assembly", default=DEFAULT_VEP_ASSEMBLY,
                       help=f"VEP assembly name (default: {DEFAULT_VEP_ASSEMBLY})")
    p_snv.add_argument("-x", "--include-sex", action="store_true",
                       help="count X/Y in addition to autosomes")
    p_snv.add_argument("-F", "--include-filtered", action="store_true",
                       help="count non-PASS records too (default: PASS only)")
    p_snv.set_defaults(func=_run_snv_indel_cli)

    p_sv = sub.add_parser(
        "sv", help="structural-variant burden from a somatic SV VCF")
    _add_common(p_sv)
    p_sv.add_argument("-i", "--vcf", required=True,
                      help="somatic SV VCF (.vcf/.vcf.gz), e.g. Manta "
                           "somaticSV.vcf.gz")
    p_sv.set_defaults(func=_run_sv_cli)

    p_cnv = sub.add_parser(
        "cnv", help="copy-number burden from a CNVkit .cns table")
    _add_common(p_cnv)
    p_cnv.add_argument("-i", "--cns", required=True,
                       help="CNVkit .cns segment table")
    p_cnv.add_argument("-m", "--min-abs-log2", type=float,
                       default=DEFAULT_MIN_ABS_LOG2,
                       help=f"altered threshold on |log2| when the table "
                            f"has no call column (default: "
                            f"{DEFAULT_MIN_ABS_LOG2})")
    p_cnv.set_defaults(func=_run_cnv_cli)

    p_agg = sub.add_parser(
        "aggregate", help="merge per-sample burden tables into one summary")
    p_agg.add_argument("-n", "--snv", nargs="*", default=[],
                       help="per-sample snv_indel burden TSVs")
    p_agg.add_argument("-v", "--sv", nargs="*", default=[],
                       help="per-sample sv burden TSVs")
    p_agg.add_argument("-c", "--cnv", nargs="*", default=[],
                       help="per-sample cnv burden TSVs")
    p_agg.add_argument("-o", "--out", required=True, help="summary TSV path")
    p_agg.add_argument("-l", "--log", default=None, type=Path,
                       help="optional log file path (default: stdout only)")
    p_agg.add_argument("-w", "--overwrite", action="store_true",
                       help="rewrite the output when it already exists")
    p_agg.set_defaults(sample="(all)", func=_run_aggregate_cli)

    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Parse CLI arguments and run the selected burden calculation.

    Parameters
    ----------
    argv : Optional[Sequence[str]], optional
        Argument list (defaults to ``sys.argv[1:]``).
    """
    global logger
    args = parse_args(argv)
    if args.log is not None:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        logger = setup_logger("Burden", log_file=str(args.log), stream=sys.stdout)
    _run_command(args)


if __name__ == "__main__":
    main()
