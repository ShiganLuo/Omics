#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import csv
import os
from collections import OrderedDict, defaultdict
from typing import Dict, Iterable, List, Tuple


def detect_delimiter(path: str) -> str:
	"""Detect delimiter from a small file sample."""
	with open(path, 'r') as fh:
		sample = fh.read(4096)
		try:
			return csv.Sniffer().sniff(sample).delimiter
		except csv.Error:
			return '\t'


def read_group_table(path: str) -> List[Tuple[str, str]]:
	"""Read sample-group pairs from a TSV/CSV table."""
	delimiter = detect_delimiter(path)
	with open(path, 'r') as fh:
		reader = csv.DictReader(fh, delimiter=delimiter)
		if reader.fieldnames is None:
			raise ValueError('Group table has no header.')
		field_map = {name.lower(): name for name in reader.fieldnames}
		if 'sample' not in field_map or 'group' not in field_map:
			raise ValueError('Group table must contain columns: sample, group')
		sample_key = field_map['sample']
		group_key = field_map['group']
		rows = []
		for row in reader:
			sample = row.get(sample_key, '').strip()
			group = row.get(group_key, '').strip()
			if sample == '' or group == '':
				continue
			rows.append((sample, group))
		return rows


def parse_int(value: str) -> int:
	"""Parse integer values with fallback to 0."""
	value = value.strip()
	if value == '' or value == '.':
		return 0
	try:
		return int(value)
	except ValueError:
		return 0


def summarize_chimeric_file(path: str, require_te: bool) -> Dict[str, int]:
	"""Summarize 5'/3' chimeric transcript counts for one sample file."""
	counts = defaultdict(int)
	with open(path, 'r') as fh:
		reader = csv.DictReader(fh, delimiter='\t')
		for row in reader:
			te_count = parse_int(row.get('TE_count', '0'))
			if require_te and te_count == 0:
				continue
			te_5p = parse_int(row.get('TE_5p_count', '0'))
			te_3p = parse_int(row.get('TE_3p_count', '0'))

			is_5p = te_5p > 0
			is_3p = te_3p > 0
			if is_5p and not is_3p:
				counts['five_only'] += 1
			elif is_3p and not is_5p:
				counts['three_only'] += 1
			elif is_5p and is_3p:
				counts['both'] += 1
			else:
				counts['none'] += 1

			if not (is_5p or is_3p):
				counts['dominant_none'] += 1
			elif te_5p > te_3p:
				counts['dominant_five'] += 1
			elif te_3p > te_5p:
				counts['dominant_three'] += 1
			else:
				counts['dominant_tie'] += 1
	return counts


def parse_te_names(value: str) -> List[str]:
	"""Split TE name field into a list, ignoring empty placeholders."""
	value = value.strip()
	if value == '' or value == '.':
		return []
	return [item for item in value.split(';') if item]


def summarize_te_types(path: str, require_te: bool) -> Dict[str, Dict[str, int]]:
	"""Count TE types at 5' and 3' ends for one sample file."""
	counts = defaultdict(lambda: {'five': 0, 'three': 0})
	with open(path, 'r') as fh:
		reader = csv.DictReader(fh, delimiter='\t')
		for row in reader:
			te_count = parse_int(row.get('TE_count', '0'))
			if require_te and te_count == 0:
				continue
			five_names = parse_te_names(row.get('TE_5p_names', '.'))
			three_names = parse_te_names(row.get('TE_3p_names', '.'))
			for name in five_names:
				counts[name]['five'] += 1
			for name in three_names:
				counts[name]['three'] += 1
	return counts


def find_input_file(input_dir: str, sample: str, suffix: str) -> str:
	"""Find the input file for a given sample."""
	path = os.path.join(input_dir, sample, f'{sample}{suffix}')
	if os.path.exists(path):
		return path
	raise FileNotFoundError(f'Missing input file for sample {sample}: {path}')


def write_summary(path: str, rows: Iterable[Dict[str, str]], header: List[str]) -> None:
	"""Write a list of dict rows to a TSV file."""
	with open(path, 'w', newline='') as fh:
		writer = csv.DictWriter(fh, fieldnames=header, delimiter='\t')
		writer.writeheader()
		for row in rows:
			writer.writerow(row)


def build_count_ticks(values: List[int]) -> List[int]:
	"""Build integer tick marks based on the maximum value."""
	max_val = max(values) if values else 0
	if max_val <= 10:
		step = 1
	elif max_val <= 50:
		step = 5
	elif max_val <= 100:
		step = 10
	else:
		step = 20
	return list(range(0, max_val + step, step))


def plot_group_proportions(path: str, group_rows: List[Dict[str, str]]) -> None:
	"""Plot stacked proportions of 5'/3' chimeric transcripts by group."""
	try:
		import matplotlib.pyplot as plt
	except ImportError as exc:
		raise RuntimeError('matplotlib is required for plotting.') from exc

	groups = [row['group'] for row in group_rows]
	five = [float(row['five_only_prop']) for row in group_rows]
	three = [float(row['three_only_prop']) for row in group_rows]
	both = [float(row['both_prop']) for row in group_rows]

	fig, ax = plt.subplots(figsize=(8, 4))
	bottom = [0.0] * len(groups)
	bar_kwargs = dict(width=0.6)
	ax.bar(groups, five, label="5' only", color='#4E79A7', **bar_kwargs)
	bottom = [b + v for b, v in zip(bottom, five)]
	ax.bar(groups, three, bottom=bottom, label="3' only", color='#F28E2B', **bar_kwargs)
	bottom = [b + v for b, v in zip(bottom, three)]
	ax.bar(groups, both, bottom=bottom, label="5' & 3'", color='#59A14F', **bar_kwargs)

	ax.set_ylabel('Proportion of chimeric transcripts')
	ax.set_ylim(0, 1)
	ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
	ax.tick_params(axis='both', labelsize=9)
	ax.set_title("5'/3' TE chimeric transcript composition")
	ax.legend(frameon=False)
	plt.tight_layout()
	fig.savefig(path, dpi=150)


def plot_te_type_top(path: str, te_rows: List[Dict[str, str]], top_n: int) -> None:
	"""Plot overall top TE types with 5'/3' counts."""
	try:
		import matplotlib.pyplot as plt
	except ImportError as exc:
		raise RuntimeError('matplotlib is required for plotting.') from exc

	rows = sorted(
		te_rows,
		key=lambda r: int(r['five_count']) + int(r['three_count']),
		reverse=True,
	)
	rows = rows[:top_n]
	if not rows:
		return

	labels = [row['te_type'] for row in rows]
	five = [int(row['five_count']) for row in rows]
	three = [int(row['three_count']) for row in rows]

	fig_height = max(4, 0.3 * len(labels) + 1.5)
	fig, ax = plt.subplots(figsize=(8, fig_height))
	y = list(range(len(labels)))
	bar_height = 0.4
	ax.barh([v - bar_height / 2 for v in y], five, height=bar_height, color='#4E79A7', label="5' count")
	ax.barh([v + bar_height / 2 for v in y], three, height=bar_height, color='#F28E2B', label="3' count")
	ax.set_yticks(y)
	ax.set_yticklabels(labels)
	ax.invert_yaxis()
	ax.set_xlabel('Transcript count')
	ax.tick_params(axis='both', labelsize=9)
	ax.set_title("Top TE types by 5'/3' counts")
	ax.legend(frameon=False)
	plt.tight_layout()
	fig.savefig(path, dpi=150)


def plot_te_type_grouped(
	path: str,
	group_te_counts: Dict[str, Dict[str, Dict[str, int]]],
	ordered_te: List[str],
) -> None:
	"""Plot multi-group TE type 5'/3' counts with shared TE ordering."""
	try:
		import matplotlib.pyplot as plt
	except ImportError as exc:
		raise RuntimeError('matplotlib is required for plotting.') from exc

	groups = list(group_te_counts.keys())
	if not groups or not ordered_te:
		return

	fig_height = max(4, 0.3 * len(ordered_te) + 1.8)
	fig, axes = plt.subplots(
		nrows=1,
		ncols=len(groups),
		figsize=(4.0 * len(groups), fig_height),
		sharey=True,
	)
	if len(groups) == 1:
		axes = [axes]

	for ax, group in zip(axes, groups):
		te_map = group_te_counts[group]
		five = [te_map.get(te, {}).get('five', 0) for te in ordered_te]
		three = [te_map.get(te, {}).get('three', 0) for te in ordered_te]
		y = list(range(len(ordered_te)))
		bar_height = 0.4
		ax.barh([v - bar_height / 2 for v in y], five, height=bar_height, color='#4E79A7', label="5' count")
		ax.barh([v + bar_height / 2 for v in y], three, height=bar_height, color='#F28E2B', label="3' count")
		ax.set_title(group)
		ax.tick_params(axis='both', labelsize=9)
		ax.invert_yaxis()
		ax.set_yticks(y)
		ax.set_yticklabels(ordered_te)

	axes[0].legend(frameon=False, loc='lower right')
	fig.suptitle("TE types: 5'/3' counts by group", y=0.98)
	plt.tight_layout()
	fig.savefig(path, dpi=150)


def main() -> None:
	parser = argparse.ArgumentParser(
		description='Summarize TE chimeric transcripts by 5\'/3\' position and plot group composition.'
	)
	parser.add_argument('-g', '--group-table', required=True, help='TSV/CSV with columns: sample, group')
	parser.add_argument('-i', '--input-dir', required=True, help='Directory with *_TE_chimeric_transcripts.txt files')
	parser.add_argument('-s', '--suffix', default='_TE_chimeric_transcripts.txt', help='Input file suffix')
	parser.add_argument('-o', '--out-prefix', required=True, help='Output prefix for TSV/plot')
	parser.add_argument('--require-te', action='store_true', help='Only count transcripts with TE_count > 0')
	parser.add_argument('--no-plot', action='store_true', help='Skip plot output')
	parser.add_argument('--te-type-top', type=int, default=20, help='Top N TE types to plot (overall)')
	parser.add_argument('--no-group-plot', action='store_true', help='Skip multi-group TE type plot')
	args = parser.parse_args()

	group_rows = read_group_table(args.group_table)
	if not group_rows:
		raise ValueError('No valid rows in group table.')

	sample_summaries = []
	group_summary = OrderedDict()
	group_te_counts = OrderedDict()
	for sample, group in group_rows:
		input_file = find_input_file(args.input_dir, sample, args.suffix)
		counts = summarize_chimeric_file(input_file, args.require_te)
		te_counts = summarize_te_types(input_file, args.require_te)

		total = counts['five_only'] + counts['three_only'] + counts['both'] + counts['none']
		total_chimeric = counts['five_only'] + counts['three_only'] + counts['both']
		row = {
			'sample': sample,
			'group': group,
			'total_tx': str(total),
			'chimeric_tx': str(total_chimeric),
			'five_only': str(counts['five_only']),
			'three_only': str(counts['three_only']),
			'both': str(counts['both']),
			'none': str(counts['none']),
			'dominant_five': str(counts['dominant_five']),
			'dominant_three': str(counts['dominant_three']),
			'dominant_tie': str(counts['dominant_tie']),
			'dominant_none': str(counts['dominant_none']),
		}
		sample_summaries.append(row)

		if group not in group_summary:
			group_summary[group] = defaultdict(int)
		if group not in group_te_counts:
			group_te_counts[group] = defaultdict(lambda: {'five': 0, 'three': 0})
		for key in counts:
			group_summary[group][key] += counts[key]
		for te_name, te_count in te_counts.items():
			group_te_counts[group][te_name]['five'] += te_count['five']
			group_te_counts[group][te_name]['three'] += te_count['three']

	sample_header = [
		'sample', 'group', 'total_tx', 'chimeric_tx',
		'five_only', 'three_only', 'both', 'none',
		'dominant_five', 'dominant_three', 'dominant_tie', 'dominant_none',
	]
	write_summary(f'{args.out_prefix}.sample_summary.tsv', sample_summaries, sample_header)

	group_rows_out = []
	for group, counts in group_summary.items():
		total_chimeric = counts['five_only'] + counts['three_only'] + counts['both']
		denom = total_chimeric if total_chimeric > 0 else 1
		row = {
			'group': group,
			'chimeric_tx': str(total_chimeric),
			'five_only': str(counts['five_only']),
			'three_only': str(counts['three_only']),
			'both': str(counts['both']),
			'five_only_prop': f"{counts['five_only'] / denom:.6f}",
			'three_only_prop': f"{counts['three_only'] / denom:.6f}",
			'both_prop': f"{counts['both'] / denom:.6f}",
			'dominant_five': str(counts['dominant_five']),
			'dominant_three': str(counts['dominant_three']),
			'dominant_tie': str(counts['dominant_tie']),
		}
		group_rows_out.append(row)

	group_header = [
		'group', 'chimeric_tx',
		'five_only', 'three_only', 'both',
		'five_only_prop', 'three_only_prop', 'both_prop',
		'dominant_five', 'dominant_three', 'dominant_tie',
	]
	write_summary(f'{args.out_prefix}.group_summary.tsv', group_rows_out, group_header)

	te_rows_out = []
	te_overall = defaultdict(lambda: {'five': 0, 'three': 0})
	for group, te_map in group_te_counts.items():
		for te_name, counts in te_map.items():
			total = counts['five'] + counts['three']
			if total == 0:
				continue
			te_rows_out.append({
				'group': group,
				'te_type': te_name,
				'five_count': str(counts['five']),
				'three_count': str(counts['three']),
				'total_count': str(total),
			})
			te_overall[te_name]['five'] += counts['five']
			te_overall[te_name]['three'] += counts['three']

	te_rows_out.sort(key=lambda r: (r['group'], r['te_type']))
	te_header = ['group', 'te_type', 'five_count', 'three_count', 'total_count']
	write_summary(f'{args.out_prefix}.te_type_counts.tsv', te_rows_out, te_header)

	overall_rows = []
	for te_name, counts in te_overall.items():
		total = counts['five'] + counts['three']
		if total == 0:
			continue
		overall_rows.append({
			'te_type': te_name,
			'five_count': str(counts['five']),
			'three_count': str(counts['three']),
			'total_count': str(total),
		})

	overall_rows_sorted = sorted(
		overall_rows,
		key=lambda r: int(r['five_count']) + int(r['three_count']),
		reverse=True,
	)
	ordered_te = [row['te_type'] for row in overall_rows_sorted[: args.te_type_top]]

	if not args.no_plot:
		plot_group_proportions(f'{args.out_prefix}.group_stacked.png', group_rows_out)
		plot_te_type_top(f'{args.out_prefix}.te_type_top.png', overall_rows, args.te_type_top)
		if not args.no_group_plot:
			plot_te_type_grouped(
				f'{args.out_prefix}.te_type_by_group.png',
				group_te_counts,
				ordered_te,
			)


if __name__ == '__main__':
	main()
