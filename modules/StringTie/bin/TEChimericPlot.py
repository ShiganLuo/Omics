#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import csv
import os
from collections import OrderedDict, defaultdict
from typing import Dict, Iterable, List, Tuple, Optional, Sequence
import numpy as np
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt


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
	"""Summarize transcript counts with TE in 5'/3'/end/internal regions."""
	counts = defaultdict(int)
	with open(path, 'r') as fh:
		reader = csv.DictReader(fh, delimiter='\t')
		for row in reader:
			te_count = parse_int(row.get('TE_count', '0'))
			if require_te and te_count == 0:
				continue
			te_5p = parse_int(row.get('TE_5p_count', '0'))
			te_3p = parse_int(row.get('TE_3p_count', '0'))
			te_end = parse_int(row.get('TE_end_unknown_count', '0'))
			te_internal = parse_int(row.get('TE_internal_count', '0'))

			if te_5p > 0:
				counts['five_any'] += 1
			if te_3p > 0:
				counts['three_any'] += 1
			if te_end > 0:
				counts['end_any'] += 1
			if te_internal > 0:
				counts['internal_any'] += 1
			if te_5p == 0 and te_3p == 0 and te_end == 0 and te_internal == 0:
				counts['none'] += 1
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



from collections import defaultdict
from typing import Dict, List, Optional

import numpy as np
from matplotlib.patches import Polygon


def plot_group_proportions(
	path: str,
	group_te_counts: Dict[str, Dict[str, int]],
	ordered_te: Optional[List[str]] = None,
	group_colors: Optional[Dict[str, str]] = None,
	fig_title:str = "",
	x_label: str = "Group",
	y_label: str = "Transcript",
	legend_title: str = "TE insert position",
	normalize: bool = False,
	fill_between_groups: bool = True,
	fill_alpha: float = 0.18,
	bar_width: float = 1,
	component_height: float = 1,
	group_gap: float = 0.35,
	dpi: int = 300,
) -> None:
	"""
	Plot stacked horizontal bar chart for TE composition across groups.

	Parameters
	----------
	path
		Output image path.

	group_te_counts
		Structure:

		{
			group_name: {
				te_name: count
			}
		}

	ordered_te
		TE display order.
		If None, automatically inferred from all groups.

	group_colors
		Color mapping for each TE component.

	normalize
		If True:
			draw proportions (0~1)

		If False:
			draw raw counts

	fill_between_groups
		Whether to fill polygons between adjacent groups
		for the same TE component.

	fill_alpha
		Polygon transparency.

	component_height
		Component height.

	group_gap
		Vertical gap between groups.

	Notes
	-----
	Each group is one stacked horizontal bar.

	Polygon filling connects the same TE component
	between adjacent groups.
	"""

	# ------------------------------------------------------------------
	# infer ordered_te
	# ------------------------------------------------------------------
	if ordered_te is None:

		all_te = set()

		for te_map in group_te_counts.values():
			all_te.update(te_map.keys())

		ordered_te = sorted(all_te)

	if not ordered_te:
		return

	group_names = list(group_te_counts.keys())
	n_groups = len(group_names)
	n_legend = len(group_te_counts[group_names[0]])


	if group_colors is None:

		palette = [
			'#4E79A7',
			'#F28E2B',
			'#59A14F',
			'#E15759',
			'#B07AA1',
			'#76B7B2',
			'#EDC948',
			'#9C755F',
			'#BAB0AC',
		]

		group_colors = {
			te: palette[i % len(palette)]
			for i, te in enumerate(ordered_te)
		}


	plot_values = defaultdict(dict)

	for group in group_names:

		te_map = group_te_counts[group]

		numeric_te_map = {}

		for te, v in te_map.items():

			try:
				numeric_te_map[te] = float(v)
			except (TypeError, ValueError):
				numeric_te_map[te] = 0.0

		total = sum(numeric_te_map.values())

		for te in ordered_te:

			value = numeric_te_map.get(te, 0.0)

			if normalize:

				if total > 0:
					value = value / total
				else:
					value = 0.0

			plot_values[group][te] = value

	# ------------------------------------------------------------------
	# figure
	# ------------------------------------------------------------------
	fig_width = n_groups * (bar_width + group_gap) + 2

	fig_height = component_height * n_legend
	fig, ax = plt.subplots(figsize=(fig_width, fig_height))


	x_positions = {}

	current_x = 0.0

	for group in group_names:

		x_positions[group] = current_x

		current_x += bar_width + group_gap

	# ------------------------------------------------------------------
	# boundaries
	#
	# boundaries[group][te] = (y0, y1)
	# ------------------------------------------------------------------
	boundaries = defaultdict(dict)

	# ------------------------------------------------------------------
	# draw stacked bars
	# ------------------------------------------------------------------
	for group in group_names:

		x = x_positions[group]

		bottom = 0.0

		for te in ordered_te:

			value = plot_values[group][te]

			y0 = bottom
			y1 = bottom + value

			boundaries[group][te] = (y0, y1)

			if value > 0:

				ax.bar(
					x=x,
					height=value,
					bottom=bottom,
					width=bar_width,
					color=group_colors[te],
					edgecolor='none',
					label=te,
					zorder=3,
				)

			bottom = y1

	# ------------------------------------------------------------------
	# polygon filling between adjacent groups
	# ------------------------------------------------------------------
	if fill_between_groups and n_groups >= 2:

		for te in ordered_te:

			for i in range(n_groups - 1):

				group_a = group_names[i]
				group_b = group_names[i + 1]

				x_a = x_positions[group_a]
				x_b = x_positions[group_b]

				y0_a, y1_a = boundaries[group_a].get(
					te,
					(0, 0),
				)

				y0_b, y1_b = boundaries[group_b].get(
					te,
					(0, 0),
				)

				# skip empty
				if (
					(y1_a - y0_a) == 0
					and
					(y1_b - y0_b) == 0
				):
					continue

				polygon = Polygon(
					[
						(x_a + bar_width / 2, y0_a),
						(x_b - bar_width / 2, y0_b),
						(x_b - bar_width / 2, y1_b),
						(x_a + bar_width / 2, y1_a),
					],
					closed=True,
					facecolor=group_colors[te],
					alpha=fill_alpha,
					edgecolor='none',
					zorder=1,
				)

				ax.add_patch(polygon)


	ax.set_xticks(
		[
			x_positions[group]
			for group in group_names
		]
	)

	ax.set_xticklabels(
		group_names,
		rotation=45,
		ha='right',
	)

	# ------------------------------------------------------------------
	# aesthetics
	# ------------------------------------------------------------------
	if normalize:
		ax.set_ylim(0, 1)
		ax.set_ylabel(y_label + " (proportion)")
	else:
		ax.set_ylabel(y_label+ "(count)")

	ax.set_xlabel(x_label)

	ax.set_title(fig_title)

	ax.grid(
		axis='y',
		linestyle='--',
		alpha=0.3,
		zorder=0,
	)

	# ------------------------------------------------------------------
	# legend
	# ------------------------------------------------------------------
	handles, labels = ax.get_legend_handles_labels()

	unique = dict(zip(labels, handles))

	ax.legend(
		unique.values(),
		unique.keys(),
		frameon=False,
		loc=(1.02, 0.5),
		title=legend_title,
	)

	plt.tight_layout()

	fig.savefig(
		path,
		dpi=dpi,
		bbox_inches='tight',
	)

	plt.close(fig)

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

		total = (
			counts['five_any']
			+ counts['three_any']
			+ counts['end_any']
			+ counts['internal_any']
			+ counts['none']
		)
		total_chimeric = (
			counts['five_any']
			+ counts['three_any']
			+ counts['end_any']
			+ counts['internal_any']
		)
		row = {
			'sample': sample,
			'group': group,
			'total_tx': str(total),
			'chimeric_tx': str(total_chimeric),
			'five_any': str(counts['five_any']),
			'three_any': str(counts['three_any']),
			'end_any': str(counts['end_any']),
			'internal_any': str(counts['internal_any']),
			'none': str(counts['none']),
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
		'five_any', 'three_any', 'end_any', 'internal_any', 'none',
	]
	write_summary(f'{args.out_prefix}.sample_summary.tsv', sample_summaries, sample_header)

	group_rows_out = {}
	for group, counts in group_summary.items():
		total_chimeric = (
			counts['five_any']
			+ counts['three_any']
			+ counts['end_any']
			+ counts['internal_any']
		)
		denom = total_chimeric if total_chimeric > 0 else 1
		row = {
			'five_any': counts['five_any'],
			'three_any': counts['three_any'],
			'end_any': counts['end_any'],
			'internal_any': counts['internal_any'],
		}
		group_rows_out[group] = row
	group_header = [
		'group', 'chimeric_tx',
		'five_any', 'three_any', 'end_any', 'internal_any',
		'five_any_prop', 'three_any_prop', 'end_any_prop', 'internal_any_prop',
	]
	# write_summary(f'{args.out_prefix}.group_summary.tsv', group_rows_out, group_header)

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
		plot_group_proportions(
			f'{args.out_prefix}.group_stacked.png',
			group_te_counts=group_rows_out,			
		)
		plot_te_type_top(f'{args.out_prefix}.te_type_top.png', overall_rows, args.te_type_top)
		if not args.no_group_plot:
			plot_te_type_grouped(
				f'{args.out_prefix}.te_type_by_group.png',
				group_te_counts,
				ordered_te,
			)


if __name__ == '__main__':
	main()
