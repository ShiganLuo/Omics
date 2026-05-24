#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import collections
from typing import Dict, Literal, Optional, Set
from intervaltree import IntervalTree

PositionType = Literal[
    'five_prime',
    'three_prime',
    'end_unknown',
    'internal',
    'outside'
]


def parse_gtf_attributes(attr: str) -> Dict[str, str]:
	"""
	Parse a GTF attribute string into a dict.

	Parameters
	----------
	attr : str
		Raw attribute string from a GTF line.

	Returns
	-------
	attrs : dict
		Dictionary of attribute key/value pairs.
	"""
	attrs = {}
	for field in attr.strip().strip(';').split(';'):
		field = field.strip()
		if not field:
			continue
		key, sep, value = field.partition(' ')
		if not sep:
			continue
		attrs[key] = value.strip().strip('"')
	return attrs

def parse_te_gtf(
	te_gtf: str,
	classification_level: Literal['gene_id', 'family_id', 'class_id'] = 'gene_id',
	chrom_filter: Optional[Set[str]] = None,
) -> Dict[str, IntervalTree]:
	"""
	Parse TE (transposon element) GTF file and build interval trees for fast overlap queries.

	Parameters
	----------
	te_gtf : str
		Path to the TE GTF annotation file.
	classification_level : str
		Attribute key to use for TE classification (e.g., 'gene_id', 'family_id', 'class_id').
	chrom_filter : set of str, optional
		If provided, only TE records on these chromosomes are parsed.

	Returns
	-------
	te_tree : dict
		Dictionary mapping chromosome to IntervalTree of TE regions. Each interval stores (te_type, te_id).
	"""
	te_tree = collections.defaultdict(IntervalTree)
	with open(te_gtf) as f:
		for line in f:
			if line.startswith('#') or line.strip() == '':
				continue
			fields = line.strip().split('\t')
			if len(fields) < 9:
				continue
			chrom, source, feature, start, end, score, strand, frame, attr = fields
			if chrom_filter is not None and chrom not in chrom_filter:
				continue
			start, end = int(start), int(end)
			attrs = parse_gtf_attributes(attr)
			te_id = attrs.get(classification_level, None)
			if te_id is None:
				continue
			te_tree[chrom].addi(start, end + 1, te_id)
	return te_tree

def parse_stringtie_gtf(gtf: str) -> Dict:
	"""
	Parse StringTie GTF file and extract transcript-level coordinates.

	Parameters
	----------
	gtf : str
		Path to the StringTie GTF file.

	Returns
	-------
	trans_coords : dict
		Dictionary mapping transcript_id to (chrom, start, end, strand, cov, fpkm, tpm).
	"""
	trans_coords = dict()  # transcript_id: (chrom, start, end, strand, cov, fpkm, tpm)
	with open(gtf) as f:
		for line in f:
			if line.startswith('#') or line.strip() == '':
				continue
			fields = line.strip().split('\t')
			if len(fields) < 9:
				continue
			chrom, source, feature, start, end, score, strand, frame, attr = fields
			if feature != 'transcript':
				continue
			start, end = int(start), int(end)
			attrs = parse_gtf_attributes(attr)
			tx_id = attrs.get('transcript_id', None)
			if tx_id:
				cov = attrs.get('cov', '.')
				fpkm = attrs.get('FPKM', '.')
				tpm = attrs.get('TPM', '.')
				trans_coords[tx_id] = (chrom, start, end, strand, cov, fpkm, tpm)
	return trans_coords

def interval_overlap(
        start1: int,
        end1: int,
        start2: int,
        end2: int) -> int:
    """
    Calculate overlap length between two 1-based closed intervals.

    Returns
    -------
    overlap_length : int
        Overlap length (0 means no overlap).
    """
    return max(0, min(end1, end2) - max(start1, start2) + 1)


def classify_te_position(
        trans_start: int,
        trans_end: int,
        strand: str,
        te_start: int,
        te_end: int,
        internal_distance: int) -> PositionType:
    """
    Classify TE position relative to transcript ends.

    Parameters
    ----------
    trans_start : int
        Transcript start (1-based, inclusive).

    trans_end : int
        Transcript end (1-based, inclusive).

    strand : str
        Transcript strand:
            '+' : positive strand
            '-' : negative strand
            '.' : unknown strand

    te_start : int
        TE start (1-based, inclusive).

    te_end : int
        TE end (1-based, inclusive).

    internal_distance : int
        Distance from transcript ends regarded as terminal region.

    Returns
    -------
    position : PositionType

        five_prime:
            TE overlaps transcript 5' terminal region

        three_prime:
            TE overlaps transcript 3' terminal region

        end_unknown:
            TE overlaps transcript end region but transcript strand unknown

        internal:
            TE overlaps transcript internal region only

        outside:
            TE does not overlap transcript
    """

    if trans_start > trans_end:
        trans_start, trans_end = trans_end, trans_start

    if te_start > te_end:
        te_start, te_end = te_end, te_start


    transcript_overlap = interval_overlap(
        trans_start,
        trans_end,
        te_start,
        te_end
    )

    if transcript_overlap == 0:
        return 'outside'

    # ------------------------------------------------------------------
    # transcript length
    # ------------------------------------------------------------------

    transcript_length = trans_end - trans_start + 1

    # terminal region size should not exceed half transcript length
    #
    # otherwise short transcripts may have fully overlapping
    # 5' and 3' regions
    # ------------------------------------------------------------------

    terminal_size = min(
        internal_distance,
        max(1, transcript_length // 2)
    )

    # ------------------------------------------------------------------
    # define transcript-relative terminal regions
    #
    # positive strand:
    #
    # 5' ---------------------> 3'
    #
    # negative strand:
    #
    # 3' <--------------------- 5'
    # ------------------------------------------------------------------

    if strand == '+':

        five_prime_region = (
            trans_start,
            trans_start + terminal_size - 1
        )

        three_prime_region = (
            trans_end - terminal_size + 1,
            trans_end
        )

    elif strand == '-':

        five_prime_region = (
            trans_end - terminal_size + 1,
            trans_end
        )

        three_prime_region = (
            trans_start,
            trans_start + terminal_size - 1
        )

    else:
        five_prime_region = None
        three_prime_region = None

    # ------------------------------------------------------------------
    # unknown strand
    # ------------------------------------------------------------------

    if strand not in ('+', '-'):

        start_region = (
            trans_start,
            trans_start + terminal_size - 1
        )

        end_region = (
            trans_end - terminal_size + 1,
            trans_end
        )

        if (
            interval_overlap(
                te_start,
                te_end,
                *start_region
            ) > 0
            or
            interval_overlap(
                te_start,
                te_end,
                *end_region
            ) > 0
        ):
            return 'end_unknown'

        return 'internal'

    # ------------------------------------------------------------------
    # calculate overlaps
    # ------------------------------------------------------------------

    five_prime_overlap = interval_overlap(
        te_start,
        te_end,
        *five_prime_region
    )

    three_prime_overlap = interval_overlap(
        te_start,
        te_end,
        *three_prime_region
    )

    # ------------------------------------------------------------------
    # classify
    # ------------------------------------------------------------------

    if five_prime_overlap > 0 and three_prime_overlap > 0:

        # short transcript or long TE
        #
        # choose side with larger overlap

        if five_prime_overlap > three_prime_overlap:
            return 'five_prime'

        elif three_prime_overlap > five_prime_overlap:
            return 'three_prime'

        else:
            # equal overlap
            #
            # prioritize 5' because promoter-associated
            # TE events are usually biologically more important

            return 'five_prime'

    if five_prime_overlap > 0:
        return 'five_prime'

    if three_prime_overlap > 0:
        return 'three_prime'

    return 'internal'


def main(stringtie_gtf, te_gtf, output, internal_distance: int):
	"""
	Main function to compute TE chimeric transcripts statistics.

	Parameters
	----------
	stringtie_gtf : str
		Path to the StringTie GTF file.
	te_gtf : str
		Path to the TE GTF annotation file.
	output : str
		Path to the output TSV file.

	Output
	------
	Writes a TSV file with TE counts and TE names by position.
	"""
	trans_coords = parse_stringtie_gtf(stringtie_gtf)
	chrom_filter = {chrom for chrom, _, _, _, _, _, _ in trans_coords.values()}
	te_tree = parse_te_gtf(te_gtf, classification_level='gene_id', chrom_filter=chrom_filter)
	with open(output, 'w') as out:
		out.write(
			'transcript_id\tchrom\tstart\tend\tstrand\tcov\tFPKM\tTPM\tTE_count\t'
			'TE_5p_count\tTE_3p_count\tTE_end_unknown_count\tTE_internal_count\t'
			'TE_5p_names\tTE_3p_names\tTE_end_unknown_names\tTE_internal_names\n'
		)
		for tx_id, (chrom, start, end, strand, cov, fpkm, tpm) in trans_coords.items():
			hits = te_tree[chrom].overlap(start, end + 1)
			if not hits:
				out.write(
					f'{tx_id}\t{chrom}\t{start}\t{end}\t{strand}\t{cov}\t{fpkm}\t{tpm}\t'
					'0\t0\t0\t0\t.\t.\t.\t.\n'
				)
				continue
			te_types = {iv.data for iv in hits}
			five_prime = set()
			three_prime = set()
			end_unknown = set()
			internal = set()
			for iv in hits:
				position = classify_te_position(start, end, strand, iv.begin, iv.end - 1, internal_distance)
				if position == 'five_prime':
					five_prime.add(iv.data)
				elif position == 'three_prime':
					three_prime.add(iv.data)
				elif position == 'end_unknown':
					end_unknown.add(iv.data)
				elif position == 'outside':
					continue
				else:
					internal.add(iv.data)
			out.write(
				f'{tx_id}\t{chrom}\t{start}\t{end}\t{strand}\t{cov}\t{fpkm}\t{tpm}\t{len(te_types)}\t'
				f'{len(five_prime)}\t{len(three_prime)}\t{len(end_unknown)}\t{len(internal)}\t'
				f'{";".join(sorted(five_prime)) if five_prime else "."}\t'
				f'{";".join(sorted(three_prime)) if three_prime else "."}\t'
				f'{";".join(sorted(end_unknown)) if end_unknown else "."}\t'
				f'{";".join(sorted(internal)) if internal else "."}\n'
			)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Calculate TE chimeric transcript statistics from StringTie GTF and TE GTF.'
	)
	parser.add_argument('-s', '--stringtie_gtf', required=True, help='Input StringTie assembled GTF file')
	parser.add_argument('-t', '--te_gtf', required=True, help='Input transposable element GTF annotation file')
	parser.add_argument('-o', '--output', required=True, help='Output TSV file for statistics')
	parser.add_argument(
		'-d', '--internal-distance',
		type=int,
		default=200,
		help='Distance (bp) from transcript ends to define 5\'/3\' regions; otherwise internal',
	)
	args = parser.parse_args()
	main(args.stringtie_gtf, args.te_gtf, args.output, args.internal_distance)
