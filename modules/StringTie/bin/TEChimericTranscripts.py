#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import collections
from typing import Dict, Literal, Optional, Set
from intervaltree import IntervalTree


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
) -> Dict:
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
	tx_coords : dict
		Dictionary mapping transcript_id to (chrom, start, end, strand).
	"""
	tx_coords = dict()  # transcript_id: (chrom, start, end, strand)
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
				tx_coords[tx_id] = (chrom, start, end, strand)
	return tx_coords

def classify_te_position(tx_start: int, tx_end: int, strand: str, te_start: int, te_end: int, internal_distance: int) -> str:
	"""
	Classify TE position relative to transcript ends.

	Parameters
	----------
	tx_start : int
		Transcript start (1-based, inclusive).
	tx_end : int
		Transcript end (1-based, inclusive).
	strand : str
		Transcript strand ('+' or '-').
	te_start : int
		TE start (1-based, inclusive).
	te_end : int
		TE end (1-based, inclusive).
	internal_distance : int
		Distance from transcript ends to consider as 5'/3'.

	Returns
	-------
	position : str
		One of 'five_prime', 'three_prime', or 'internal'.
	"""
	pos = (te_start + te_end) // 2
	if strand == '+':
		if pos - tx_start <= internal_distance:
			return 'five_prime'
		if tx_end - pos <= internal_distance:
			return 'three_prime'
		return 'internal'
	if tx_end - pos <= internal_distance:
		return 'five_prime'
	if pos - tx_start <= internal_distance:
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
	tx_coords = parse_stringtie_gtf(stringtie_gtf)
	chrom_filter = {chrom for chrom, _, _, _ in tx_coords.values()}
	te_tree = parse_te_gtf(te_gtf, classification_level='gene_id', chrom_filter=chrom_filter)
	with open(output, 'w') as out:
		out.write(
			'transcript_id\tTE_count\t'
			'TE_5p_count\tTE_3p_count\tTE_internal_count\t'
			'TE_5p_names\tTE_3p_names\tTE_internal_names\n'
		)
		for tx_id, (chrom, start, end, strand) in tx_coords.items():
			hits = te_tree[chrom].overlap(start, end + 1)
			if not hits:
				out.write(f'{tx_id}\t0\t.\t0\t0\t0\t.\t.\t.\n')
				continue
			te_types = {iv.data for iv in hits}
			five_prime = set()
			three_prime = set()
			internal = set()
			for iv in hits:
				position = classify_te_position(start, end, strand, iv.begin, iv.end - 1, internal_distance)
				if position == 'five_prime':
					five_prime.add(iv.data)
				elif position == 'three_prime':
					three_prime.add(iv.data)
				else:
					internal.add(iv.data)
			out.write(
				f'{tx_id}\t{len(te_types)}\t'
				f'{len(five_prime)}\t{len(three_prime)}\t{len(internal)}\t'
				f'{";".join(sorted(five_prime)) if five_prime else "."}\t'
				f'{";".join(sorted(three_prime)) if three_prime else "."}\t'
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
