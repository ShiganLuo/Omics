#!/usr/bin/env python3
"""Per-fragment RNA/background exonic-density ratio for ecDNA fragments.

ratio = (RNA reads over fragment / fragment∩exonic kb)
      / (RNA reads over linear background / background∩exonic kb)

RNA-seq signal only lands on exons, so a raw interval-RPK ratio conflates
transcription with exon density (project data: fragments 0-14.4% exonic,
flank windows 1.0-10.9%). Both denominators here are therefore the EXONIC
intersection lengths, not interval lengths. Exon intervals come from a
GENCODE GTF (all exon features, merged per chromosome).

Background = [start-flank, end+flank] on the same chromosome MINUS every
ecDNA interval of the same sample (all circles in the decoil filtered bed).
Masking is required: for short fragments the naive left-flank window is
mostly other fragments of the same circle, which inflates the background
with circle transcription and falsely calls active fragments silent
(2026-10 recheck: homolog-ring M/S fragments 0.2-0.3x -> 1.4-3.0x).

Reads are deduplicated by query_name per window; absolute library size
cancels in the ratio.

Evidence flags (column ``evidence``):
    flank_empty   background window has zero RNA reads
    no_exonic     fragment or background has zero exonic intersection
                  (non-coding / intergenic circle content, not missing data)
    low_conf      background exonic intersection < --min-exonic-kb (poisson
                  noise dominates; excluded from strength ranking)
    rna_above_flank / rna_equal_flank / rna_below_flank

This is supporting evidence only — see METHOD.md (证据 3) for the decision
tree and the exon-normalization rationale.

References (METHOD.md):
    Gao et al., Front Oncol 2025, doi:10.3389/fonc.2025.1701989
    Turner et al., Nature 2017, doi:10.1038/nature21356
"""
from __future__ import annotations

import argparse
import bisect
import logging
import os
import sys

try:
    import pysam
except ImportError:
    pysam = None  # user installs pysam themselves (repo convention)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("rna_dna_ratio")

#: Background windows whose exonic intersection falls below this many kb are
#: flagged low_conf (default 20 kb; too few exonic bases -> noisy ratio).
DEFAULT_MIN_EXONIC_KB = 5.0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Build the CLI parser.

    Parameters
    ----------
    argv : list of str, optional
        Argument tokens; ``None`` uses ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        Parsed arguments (ecdna_bed, rna_bam, gtf, flank, min_mapq,
        min_exonic_kb, out).
    """
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--ecdna-bed", required=True,
                        help="decoil reconstruct.ecDNA.filtered.bed (coverage in col 7)")
    parser.add_argument("--rna-bam", required=True,
                        help="RNA-seq BAM aligned to the linear genome")
    parser.add_argument("--gtf", required=True,
                        help="GENCODE GTF (same assembly as the BAM); exon "
                             "features are merged per chromosome to define "
                             "exonic space")
    parser.add_argument("--flank", type=int, default=50000,
                        help="Flanking window (bp) each side for linear background "
                             "(default: %(default)s)")
    parser.add_argument("--min-mapq", type=int, default=10,
                        help="Minimum MAPQ for RNA reads (default: %(default)s)")
    parser.add_argument("--min-exonic-kb", type=float,
                        default=DEFAULT_MIN_EXONIC_KB,
                        help="Flag windows whose exonic intersection is below "
                             "this many kb as low_conf (default: %(default)s)")
    parser.add_argument("--out", required=True, help="Output TSV path")
    return parser.parse_args(argv)


def load_fragments(bed_path: str) -> list[dict]:
    """Parse the decoil filtered bed into fragment dicts.

    Parameters
    ----------
    bed_path : str
        Path to ``reconstruct.ecDNA.filtered.bed`` (coverage in column 7).

    Returns
    -------
    list of dict
        One dict per row with keys ``chrom`` (str), ``start`` (int, 0-based),
        ``end`` (int), ``circ_id`` (str), ``fragment_id`` (str),
        ``dna_coverage`` (float).
    """
    fragments = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            fragments.append({
                "chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                "circ_id": f[3], "fragment_id": f[4],
                "dna_coverage": float(f[6]) if f[6] else 0.0,
            })
    return fragments


def load_exons(gtf_path: str) -> dict[str, list[tuple[int, int]]]:
    """Extract all exon intervals from a GENCODE GTF, merged per chromosome.

    Parameters
    ----------
    gtf_path : str
        Path to a GENCODE GTF (GTF feature column must contain ``exon``
        rows; coordinates are 1-based inclusive and converted to 0-based
        half-open).

    Returns
    -------
    dict
        ``chrom -> sorted merged [(start, end), ...]`` (0-based half-open,
        non-overlapping).
    """
    per_chrom: dict[str, list[tuple[int, int]]] = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[2] != "exon":
                continue
            per_chrom.setdefault(f[0], []).append((int(f[3]) - 1, int(f[4])))
    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, ivs in per_chrom.items():
        ivs.sort()
        out: list[tuple[int, int]] = []
        for s, e in ivs:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    logger.info("Loaded exons for %d chromosomes from %s", len(merged), gtf_path)
    return merged


def exonic_length(exons: dict[str, list[tuple[int, int]]], chrom: str,
                  start: int, end: int) -> int:
    """Compute the exonic intersection length of [start, end) on a chromosome.

    Parameters
    ----------
    exons : dict
        Output of :func:`load_exons`.
    chrom : str
        Chromosome name (must match the GTF/BAM naming).
    start, end : int
        Query interval (0-based half-open).

    Returns
    -------
    int
        Sum of merged-exon bases inside the query interval (0 if the
        chromosome has no exons).
    """
    ivs = exons.get(chrom, [])
    if not ivs:
        return 0
    i = bisect.bisect_right([s for s, _ in ivs], start) - 1
    if i < 0:
        i = 0
    total = 0
    for s, e in ivs[i:]:
        if s >= end:
            break
        total += max(0, min(e, end) - max(s, start))
    return total


def count_rna_reads(bam, chrom: str, start: int, end: int, min_mapq: int,
                    names: set | None = None) -> int:
    """Count RNA fragments overlapping [start, end) with MAPQ >= min_mapq.

    Deduplicates by query_name so a fragment with two mates inside the window
    counts once. Pass an existing ``names`` set to accumulate across multiple
    windows (returns the running unique count).

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Opened RNA-seq BAM.
    chrom : str
        Chromosome name.
    start, end : int
        Window (0-based half-open).
    min_mapq : int
        Minimum mapping quality.
    names : set of str, optional
        Running set of query names shared across windows.

    Returns
    -------
    int
        Unique fragment count so far.
    """
    if names is None:
        names = set()
    for aln in bam.fetch(chrom, max(0, start), end):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if (aln.mapping_quality or 0) < min_mapq:
            continue
        names.add(aln.query_name)
    return len(names)


def background_windows(frag: dict, all_fragments: list[dict],
                       flank: int) -> list[tuple[int, int]]:
    """Left+right flank windows of ``frag`` minus every ecDNA interval.

    Subtracts the merged spans of ALL circles of the sample on the same
    chromosome (including the fragment's own circle), so only genuine linear
    genome remains.

    Parameters
    ----------
    frag : dict
        The fragment in question (``chrom``/``start``/``end`` keys).
    all_fragments : list of dict
        Every fragment of the sample (masking source).
    flank : int
        Flanking window size in bp (each side).

    Returns
    -------
    list of tuple(int, int)
        Non-overlapping masked windows (0-based half-open).
    """
    lo = max(0, frag["start"] - flank)
    hi = frag["end"] + flank
    blocked = sorted((f["start"], f["end"]) for f in all_fragments
                     if f["chrom"] == frag["chrom"])
    pieces: list[tuple[int, int]] = []
    cur = lo
    for bs, be in blocked:
        if be <= cur or bs >= hi:
            continue
        if bs > cur:
            pieces.append((cur, min(bs, hi)))
        cur = max(cur, be)
        if cur >= hi:
            break
    if cur < hi:
        pieces.append((cur, hi))
    return pieces


def main(argv: list[str] | None = None) -> int:
    """Run the per-fragment exonic-ratio computation.

    For every fragment: counts RNA reads and exonic bases over the fragment
    and over the circle-masked background windows, then writes one TSV row
    with both raw-interval and exonic-normalized ratios plus an evidence
    flag (see module docstring).

    Parameters
    ----------
    argv : list of str, optional
        Argument tokens; ``None`` uses ``sys.argv[1:]``.

    Returns
    -------
    int
        Process exit code (0 on success, 1 when pysam is missing).
    """
    args = parse_args(argv)
    if pysam is None:
        logger.error("pysam is required: pip install pysam (user-managed, no auto-install)")
        return 1

    fragments = load_fragments(args.ecdna_bed)
    logger.info("Loaded %d fragments from %s", len(fragments), args.ecdna_bed)
    exons = load_exons(args.gtf)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with pysam.AlignmentFile(args.rna_bam, "rb") as bam, open(args.out, "w") as fh:
        fh.write("circ_id\tfragment_id\tfragment\tdna_coverage\trna_reads\t"
                 "frag_exonic_kb\trna_per_exonic_kb\t"
                 "flank_len_kb\tflank_exonic_kb\tflank_reads\tflank_per_exonic_kb\t"
                 "ratio_vs_flank\tevidence\n")
        for frag in fragments:
            length = frag["end"] - frag["start"]
            rna = count_rna_reads(bam, frag["chrom"], frag["start"], frag["end"],
                                  args.min_mapq)
            frag_ex_kb = exonic_length(exons, frag["chrom"],
                                       frag["start"], frag["end"]) / 1000.0
            # Exonic density (reads per exonic kb): the RNA-seq-valid
            # denominator. Interval length is kept only as a reference column.
            frag_e_rpk = rna / frag_ex_kb if frag_ex_kb > 0 else 0.0
            # Linear background: both flanks, circle intervals masked out.
            names: set = set()
            flank_len = 0
            flank_ex_len = 0
            for s, e in background_windows(frag, fragments, args.flank):
                flank_len += e - s
                flank_ex_len += exonic_length(exons, frag["chrom"], s, e)
                count_rna_reads(bam, frag["chrom"], s, e, args.min_mapq,
                                names=names)
            flank_ex_kb = flank_ex_len / 1000.0
            flank_e_rpk = len(names) / flank_ex_kb if flank_ex_kb > 0 else 0.0
            ratio = (frag_e_rpk / flank_e_rpk) if flank_e_rpk > 0 else float("inf")
            if flank_e_rpk <= 0 and len(names) == 0:
                evidence = "flank_empty"
            elif frag_ex_kb <= 0 or flank_ex_kb <= 0:
                evidence = "no_exonic"
            elif flank_ex_kb < args.min_exonic_kb:
                evidence = "low_conf"
            else:
                evidence = ("rna_above_flank" if ratio >= 2 else
                            "rna_equal_flank" if ratio >= 1 else
                            "rna_below_flank")
            fh.write(f"{frag['circ_id']}\t{frag['fragment_id']}\t"
                     f"{frag['chrom']}:{frag['start']}-{frag['end']}\t"
                     f"{frag['dna_coverage']}\t{rna}\t"
                     f"{frag_ex_kb:.2f}\t{frag_e_rpk:.2f}\t"
                     f"{flank_len / 1000:.1f}\t{flank_ex_kb:.2f}\t"
                     f"{len(names)}\t{flank_e_rpk:.2f}\t"
                     f"{ratio:.2f}\t{evidence}\n")

    logger.info("Wrote %s", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
