from __future__ import annotations
from typing import List, Tuple, Dict, Optional, Iterable
import pysam
import logging
from dataclasses import dataclass

try:
    from .LogUtil import setup_logger
except Exception:
    from LogUtil import setup_logger


# =========================
# Data structures
# =========================

@dataclass
class Breakpoint:
    chrom: str
    pos: int


@dataclass
class TRARecord:
    id: str
    chr1: str
    pos1: int
    chr2: str
    pos2: int


@dataclass
class TERange:
    chrom: str
    start: int
    end: int
    te_name: str


# =========================
# VCF parsing
# =========================

def parse_tra_from_vcf(vcf_path: str, logger: logging.Logger) -> List[TRARecord]:
    """
    Parse TRA (translocation) records from a VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to VCF file.
    logger : logging.Logger
        Logger instance.

    Returns
    -------
    List[TRARecord]
        List of parsed TRA records.
    """
    vcf = pysam.VariantFile(vcf_path)
    results: List[TRARecord] = []

    for record in vcf:
        svtype = record.info.get("SVTYPE", None)
        if svtype != "TRA":
            continue

        try:
            chr1 = record.chrom
            pos1 = record.pos

            # standard BND format parsing
            alt = record.alts[0]
            chr2, pos2 = _parse_bnd_alt(alt)

            results.append(
                TRARecord(
                    id=record.id,
                    chr1=chr1,
                    pos1=pos1,
                    chr2=chr2,
                    pos2=pos2,
                )
            )
        except Exception as e:
            logger.warning(f"Failed parsing record {record.id}: {e}")

    logger.info(f"Parsed {len(results)} TRA records")
    return results


def _parse_bnd_alt(alt: str) -> Tuple[str, int]:
    """
    Parse BND ALT field to extract partner breakpoint.

    Example:
        N]chr2:12345] -> (chr2, 12345)
    """
    import re
    match = re.search(r'[\[\]](.+):(\d+)[\[\]]', alt)
    if not match:
        raise ValueError(f"Invalid BND ALT: {alt}")
    return match.group(1), int(match.group(2))


# =========================
# TE annotation
# =========================

def load_te_bed(bed_path: str) -> List[TERange]:
    """
    Load TE annotation BED.

    Expected format:
        chrom start end te_name
    """
    te_list: List[TERange] = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, start, end, name = line.strip().split()[:4]
            te_list.append(
                TERange(chrom, int(start), int(end), name)
            )
    return te_list


def query_te_overlap(
    chrom: str,
    pos: int,
    te_list: List[TERange],
    window: int = 200
) -> List[TERange]:
    """
    Query TE overlaps within a window around breakpoint.
    """
    result: List[TERange] = []
    for te in te_list:
        if te.chrom != chrom:
            continue
        if te.start <= pos + window and te.end >= pos - window:
            result.append(te)
    return result


# =========================
# BAM evidence
# =========================

def extract_reads(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos: int,
    window: int = 200
) -> Iterable[pysam.AlignedSegment]:
    """
    Fetch reads around breakpoint.
    """
    return bam.fetch(chrom, max(0, pos - window), pos + window)


def count_support_reads(
    reads: Iterable[pysam.AlignedSegment]
) -> Dict[str, int]:
    """
    Count split and discordant reads.
    """
    split = 0
    discordant = 0

    for r in reads:
        if r.is_unmapped:
            continue

        # split read (soft clip)
        if any(op in (4, 5) for op, _ in (r.cigartuples or [])):
            split += 1

        # discordant pair
        if not r.is_proper_pair:
            discordant += 1

    return {
        "split": split,
        "discordant": discordant
    }


# =========================
# Core logic
# =========================

def is_te_mediated_tra(
    tra: TRARecord,
    bam_path: str,
    te_list: List[TERange],
    logger: logging.Logger,
    window: int = 200,
    min_split: int = 2
) -> bool:
    """
    Determine whether a TRA event is likely TE-mediated.

    Parameters
    ----------
    tra : TRARecord
        Translocation record.
    bam_path : str
        Path to BAM file.
    te_list : List[TERange]
        TE annotation.
    logger : logging.Logger
        Logger.
    window : int, optional
        Window size around breakpoint.
    min_split : int, optional
        Minimum split reads threshold.

    Returns
    -------
    bool
        True if likely TE-mediated.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # TE overlap
    te1 = query_te_overlap(tra.chr1, tra.pos1, te_list, window)
    te2 = query_te_overlap(tra.chr2, tra.pos2, te_list, window)

    # BAM evidence
    reads1 = extract_reads(bam, tra.chr1, tra.pos1, window)
    reads2 = extract_reads(bam, tra.chr2, tra.pos2, window)

    support1 = count_support_reads(reads1)
    support2 = count_support_reads(reads2)

    logger.info(
        f"{tra.id}: TE1={len(te1)}, TE2={len(te2)}, "
        f"split=({support1['split']},{support2['split']})"
    )

    # decision rules
    if len(te1) > 0 and len(te2) > 0:
        if support1["split"] >= min_split or support2["split"] >= min_split:
            return True

    return False


# =========================
# Batch processing
# =========================

def annotate_tra_te(
    vcf_path: str,
    bam_path: str,
    te_bed: str,
    logger: Optional[logging.Logger] = None
) -> Dict[str, bool]:
    """
    Annotate all TRA events for TE mediation.

    Returns
    -------
    Dict[str, bool]
        Mapping from TRA ID to TE-mediated flag.
    """
    if logger is None:
        logger = setup_logger(__name__)

    te_list = load_te_bed(te_bed)
    tra_list = parse_tra_from_vcf(vcf_path, logger)

    results: Dict[str, bool] = {}

    for tra in tra_list:
        flag = is_te_mediated_tra(
            tra, bam_path, te_list, logger
        )
        results[tra.id] = flag

    logger.info("Annotation completed")
    return results