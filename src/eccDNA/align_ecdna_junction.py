#!/usr/bin/env python3
"""Count junction-spanning RNA reads on reconstructed ecDNA reference sequences.

Evidence level 2 in METHOD.md — the core method. RNA-seq reads are aligned to
the decoil-reconstructed ecDNA FASTA (each record is a linearized circle: the
fragments concatenated head-to-tail). Reads that cross the linearization point
(the circular junction) WITHOUT a spliced alignment are direct evidence that
the ecDNA molecule is transcribed; no read from the linear genome can produce
this pattern.

Read categories per ecDNA (record) and per internal junction between fragments:
    junction_spanning_nonsplice  -- continuous alignment across the circle
        junction (no N in CIGAR). Direct evidence of ecDNA transcription.
    junction_spanning_backsplice -- alignment crossing the junction with a
        canonical GT-AG splice. Likely a circRNA back-splicing artifact, NOT
        evidence of ecDNA transcription. Reported separately.
    internal_reads               -- reads fully inside the record (background,
        used for normalization).

Prerequisites (see METHOD.md):
    STAR --runMode genomeGenerate --genomeDir ecdna_index \
         --genomeFastaFiles reconstruct.ecDNA.filtered.fasta \
         --sjdbOverhang 99 --genomeSAindexNbases 8
    STAR --genomeDir ecdna_index --readFilesIn R1 R2 --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 \
         --outFilterMultimapNmax 100 --outFilterMismatchNoverReadLmax 0.04

Junction localization: for record of length L, the circular junction sits at
position 1|L. A read spans the junction when its alignment wraps the
linearization point. The script derives junction positions from the ecDNA BED
(fragment boundaries concatenated) rather than assuming a single junction per
record, so multi-fragment circles report per-fragment-junction counts too.

References (METHOD.md):
    Zhu et al., Nat Genet 2023, doi:10.1038/s41588-023-01386-y
    Møller et al., Nucleic Acids Res 2018, doi:10.1093/nar/gky352
    Dobin et al., Bioinformatics 2013, doi:10.1093/bioinformatics/bts635
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from collections import defaultdict

try:
    import pysam
except ImportError:
    pysam = None  # user installs pysam themselves (repo convention: no auto-install)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("align_ecdna_junction")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--ecdna-fasta", required=True,
                        help="decoil reconstruct.ecDNA.filtered.fasta (records = linearized circles)")
    parser.add_argument("--ecdna-bed", required=True,
                        help="decoil reconstruct.ecDNA.filtered.bed (fragment layout per circ_id)")
    parser.add_argument("--bam", required=True,
                        help="RNA-seq BAM aligned to the ecDNA index")
    parser.add_argument("--min-mapq", type=int, default=10,
                        help="Minimum MAPQ for a read to be counted (default: 10)")
    parser.add_argument("--min-anchor", type=int, default=10,
                        help="Minimum bases on each side of a junction for the read to count (default: 10)")
    parser.add_argument("--out", required=True, help="Output TSV path")
    return parser.parse_args(argv)


def load_record_lengths(fasta_path: str) -> dict[str, int]:
    """Record name -> total length, from FASTA headers and sequence lengths."""
    lengths: dict[str, int] = {}
    name = None
    total = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0]
                total = 0
            else:
                total += len(line.strip())
    if name is not None:
        lengths[name] = total
    logger.info("Loaded %d ecDNA records from %s", len(lengths), fasta_path)
    return lengths


def load_fragment_boundaries(bed_path: str) -> dict[str, list[tuple[int, int]]]:
    """circ_id -> ordered list of (fragment_length, cumulative_start) on the record.

    decoil concatenates fragments head-to-tail in BED order, so the k-th
    fragment of a record occupies [cum_start, cum_start + frag_len) and the
    inter-fragment junctions sit at each cum_start (k > 0). The circular
    junction sits at record position 0 == record_length.
    """
    by_circ: dict[str, list[int]] = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            circ_id = fields[3]
            by_circ[circ_id].append(int(fields[2]) - int(fields[1]))
    layout: dict[str, list[tuple[int, int]]] = {}
    for circ_id, frag_lens in by_circ.items():
        entries = []
        cum = 0
        for flen in frag_lens:
            entries.append((flen, cum))
            cum += flen
        layout[circ_id] = entries
    return layout


def iter_junctions(record: str, length: int,
                   layout: list[tuple[int, int]]) -> list[tuple[int, str]]:
    """Junction positions (0-based: a read alignment covering position x-1|x
    spans it) with labels. Junction 0 is the circular junction (record end ->
    start); junctions k>0 are internal fragment joints."""
    junctions = [(0, "circle_junction")]
    for i in range(1, len(layout)):
        junctions.append((layout[i][1], f"fragment_junction_{i}"))
    return junctions


def read_blocks(aln) -> list[tuple[int, int]]:
    """Aligned reference blocks of a read, [(start, end), ...] sorted.

    pysam's get_blocks() merges M/=/X runs and splits on both D and N. A
    deletion (D) is genuine genomic sequence within a fragment, so blocks
    separated by D are merged back; blocks separated by N (splice) stay apart.
    """
    blocks = sorted((s, e) for s, e in aln.get_blocks() if e > s)
    if not blocks:
        return []
    merged: list[list[int]] = []
    cigar = aln.cigartuples or []
    # walk the CIGAR to know which gaps are D (merge) vs N (keep apart)
    gaps: list[tuple[int, int, int]] = []  # (gap_start, gap_end, op)
    refpos = blocks[0][0] if blocks else 0
    for op, length in cigar:
        if op in (0, 7, 8):  # M, =, X
            refpos += length
        elif op in (2, 3):  # D or N
            gaps.append((refpos, refpos + length, op))
            refpos += length
    for s, e in blocks:
        if merged and s - merged[-1][1] > 0:
            gap = next((g for g in gaps if g[0] == merged[-1][1] and g[1] == s), None)
            if gap and gap[2] == 2:  # D gap: merge
                merged[-1][1] = e
                continue
        merged.append([s, e])
    return [(s, e) for s, e in merged]


def splice_motifs_at(aln) -> bool:
    """True if the alignment contains an N (splice) operation."""
    return any(op == 3 for op, _ in (aln.cigartuples or []))


def classify_read(aln, junctions: list[tuple[int, str]], length: int,
                  min_anchor: int) -> list[tuple[str, bool]]:
    """Return [(junction_label, is_backsplice), ...] for junctions the read spans.

    A read spans junction at position p when it has aligned bases on both sides
    with at least *min_anchor* bases anchored each side (for the circular
    junction at 0, the read must wrap: end of record -> start of record).
    is_backsplice is True when the read's CIGAR contains N at that junction.
    """
    blocks = read_blocks(aln)
    if not blocks:
        return []
    spans: list[tuple[str, bool]] = []
    has_n = splice_motifs_at(aln)
    for pos, label in junctions:
        if pos == 0:
            left = [b for b in blocks if b[1] >= length - min_anchor]
            right = [b for b in blocks if b[0] <= min_anchor]
            wraps = bool(left) and bool(right)
            if not wraps and len(blocks) >= 1:
                # a single block crossing 0 is impossible in linear coords;
                # soft-clip/supplementary split means primary covers one side
                wraps = False
            if wraps:
                spans.append((label, False))
            elif has_n:
                # N-based wrap: first block near start, last block near end
                first, last = blocks[0], blocks[-1]
                if (first[0] <= min_anchor
                        and last[1] >= length - min_anchor
                        and len(blocks) >= 2):
                    spans.append((label, True))
        else:
            # continuous: one block extends min_anchor bases to both sides
            contiguous = any(b[0] <= pos - min_anchor and b[1] >= pos + min_anchor
                             for b in blocks)
            if contiguous:
                spans.append((label, False))
            elif has_n:
                # spliced: the read's N lands on the junction — either the donor
                # (left block ends at pos) or the acceptor (right block starts
                # at pos) side, with min_anchor anchoring on the block itself.
                # split reads abut the junction: left block ends exactly at
                # pos (donor) or right block starts exactly at pos (acceptor),
                # each with >= min_anchor anchored bases.
                donor = any(b[1] == pos and pos - b[0] >= min_anchor
                            for b in blocks)
                acceptor = any(b[0] == pos and b[1] - pos >= min_anchor
                               for b in blocks)
                if donor or acceptor:
                    spans.append((label, True))
    return spans


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if pysam is None:
        logger.error("pysam is required: pip install pysam (user-managed, no auto-install)")
        return 1

    lengths = load_record_lengths(args.ecdna_fasta)
    layout = load_fragment_boundaries(args.ecdna_bed)
    junctions_by_record = {
        rec: iter_junctions(rec, lengths[rec], layout.get(rec, [(lengths[rec], 0)]))
        for rec in lengths
    }

    counts: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: {"nonsplice": 0, "backsplice": 0, "read_names": set()})
    internal: dict[str, int] = defaultdict(int)
    total_seen = 0

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for aln in bam:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if (aln.mapping_quality or 0) < args.min_mapq:
                continue
            record = aln.reference_name
            if record not in lengths:
                continue
            total_seen += 1
            spans = classify_read(aln, junctions_by_record[record],
                                  lengths[record], args.min_anchor)
            if spans:
                for label, is_backsplice in spans:
                    key = (record, label)
                    if is_backsplice:
                        counts[key]["backsplice"] += 1
                    else:
                        counts[key]["nonsplice"] += 1
                        counts[key]["read_names"].add(aln.query_name)
            else:
                internal[record] += 1

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("circ_id\tjunction\tnonsplice_reads\tbacksplice_reads\t"
                 "nonsplice_read_names\tinternal_reads\trecord_len_bp\tevidence\n")
        for record in sorted(lengths):
            rec_internal = internal.get(record, 0)
            written = False
            for pos, label in junctions_by_record[record]:
                key = (record, label)
                non = counts[key]["nonsplice"]
                back = counts[key]["backsplice"]
                if non or back:
                    evidence = ("ecDNA_transcribed" if label == "circle_junction"
                                and non >= 3 else
                                "junction_reads_detected" if non else
                                "backsplice_only_suspect_circRNA")
                    fh.write(f"{record}\t{label}\t{non}\t{back}\t"
                             f"{';'.join(sorted(counts[key]['read_names'])[:20])}\t"
                             f"{rec_internal}\t{lengths[record]}\t{evidence}\n")
                    written = True
            if not written:
                fh.write(f"{record}\t-\t0\t0\t-\t{rec_internal}\t"
                         f"{lengths[record]}\tno_junction_reads\n")

    logger.info("Aligned reads scanned: %d -> %s", total_seen, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
