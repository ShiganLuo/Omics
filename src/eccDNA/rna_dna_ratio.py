#!/usr/bin/env python3
"""Per-fragment RNA/DNA coverage ratio for ecDNA fragments (evidence level 3).

ratio = RNA-seq RPKM over the fragment / decoil DNA coverage of the fragment.

High ecDNA copy number inflates both coverages proportionally, so a ratio
clearly above the fragment's own flanking linear background (± flank bp)
suggests active transcription rather than copy-number artifact. This is
supporting evidence only — see METHOD.md for the decision tree.

References (METHOD.md):
    Gao et al., Front Oncol 2025, doi:10.3389/fonc.2025.1701989
    Turner et al., Nature 2017, doi:10.1038/nature21356
"""
from __future__ import annotations

import argparse
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


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--ecdna-bed", required=True,
                        help="decoil reconstruct.ecDNA.filtered.bed (coverage in col 7)")
    parser.add_argument("--rna-bam", required=True,
                        help="RNA-seq BAM aligned to the linear genome")
    parser.add_argument("--flank", type=int, default=50000,
                        help="Flanking window (bp) for linear background (default: 50000)")
    parser.add_argument("--min-mapq", type=int, default=10,
                        help="Minimum MAPQ for RNA reads (default: 10)")
    parser.add_argument("--out", required=True, help="Output TSV path")
    return parser.parse_args(argv)


def load_fragments(bed_path: str) -> list[dict]:
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


def count_rna_reads(bam, chrom: str, start: int, end: int, min_mapq: int) -> int:
    """Count RNA fragments overlapping [start, end) with MAPQ >= min_mapq.

    Uses fetch(until_eof=False) + dedup by query_name so a fragment with two
    mates inside the window counts once.
    """
    names = set()
    for aln in bam.fetch(chrom, max(0, start), end):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if (aln.mapping_quality or 0) < min_mapq:
            continue
        names.add(aln.query_name)
    return len(names)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if pysam is None:
        logger.error("pysam is required: pip install pysam (user-managed, no auto-install)")
        return 1

    fragments = load_fragments(args.ecdna_bed)
    logger.info("Loaded %d fragments from %s", len(fragments), args.ecdna_bed)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with pysam.AlignmentFile(args.rna_bam, "rb") as bam, open(args.out, "w") as fh:
        fh.write("circ_id\tfragment_id\tfragment\tdna_coverage\trna_reads\t"
                 "rna_reads_per_kb\tflank_reads_per_kb\tratio_vs_flank\tevidence\n")
        for frag in fragments:
            length = frag["end"] - frag["start"]
            rna = count_rna_reads(bam, frag["chrom"], frag["start"], frag["end"],
                                  args.min_mapq)
            # Reads-per-kb (rpk) proxy; the flank window is the internal
            # control, so absolute library size cancels in the ratio.
            rna_rpk = rna / (length / 1000) if length else 0.0
            fstart = max(0, frag["start"] - args.flank)
            fend = frag["start"]  # left flank only: avoids neighboring fragment on circle
            flank_reads = count_rna_reads(bam, frag["chrom"], fstart, fend,
                                          args.min_mapq)
            flank_len = fend - fstart
            flank_rpk = flank_reads / (flank_len / 1000) if flank_len else 0.0
            ratio = (rna_rpk / flank_rpk) if flank_rpk > 0 else float("inf")
            # DNA-side comparison: circle coverage vs (unknown) linear coverage
            # is not available here; report the RNA ratio and flag.
            evidence = ("rna_above_flank" if ratio >= 2 else
                        "rna_equal_flank" if ratio >= 1 else "rna_below_flank")
            fh.write(f"{frag['circ_id']}\t{frag['fragment_id']}\t"
                     f"{frag['chrom']}:{frag['start']}-{frag['end']}\t"
                     f"{frag['dna_coverage']}\t{rna}\t{rna_rpk:.2f}\t"
                     f"{flank_rpk:.2f}\t{ratio:.2f}\t{evidence}\n")

    logger.info("Wrote %s", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
