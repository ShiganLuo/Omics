from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Set


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


@dataclass
class GTFStats:
	total_lines: int = 0
	header_lines: int = 0
	malformed_lines: int = 0
	kept_feature_lines: int = 0
	kept_features: int = 0


def _open_text_auto(path: Path):
	"""Open plain/gzip text file for reading."""
	if str(path).endswith(".gz"):
		return gzip.open(path, "rt", encoding="utf-8", errors="replace")
	return open(path, "r", encoding="utf-8", errors="replace")


def _ensure_parent(path: Path) -> None:
	path.parent.mkdir(parents=True, exist_ok=True)


def _run_command(cmd: list[str], stdout_path: Optional[Path] = None) -> None:
	"""Run external command, optionally writing stdout to a file."""
	if stdout_path is None:
		subprocess.run(cmd, check=True)
		return

	_ensure_parent(stdout_path)
	with open(stdout_path, "wb") as out_fh:
		subprocess.run(cmd, check=True, stdout=out_fh)


def _filter_features(
	input_gtf: Path,
	output_gtf: Path,
	keep_features: Set[str],
	keep_header: bool,
) -> GTFStats:
	"""Filter GTF rows by feature type (column 3)."""
	stats = GTFStats(kept_features=len(keep_features))
	_ensure_parent(output_gtf)

	with _open_text_auto(input_gtf) as in_fh, open(output_gtf, "w", encoding="utf-8") as out_fh:
		for line in in_fh:
			stats.total_lines += 1

			if line.startswith("#"):
				stats.header_lines += 1
				if keep_header:
					out_fh.write(line)
				continue

			parts = line.rstrip("\n").split("\t")
			if len(parts) < 9:
				stats.malformed_lines += 1
				continue

			if parts[2] in keep_features:
				out_fh.write(line)
				stats.kept_feature_lines += 1

	return stats


def _sort_gtf(input_gtf: Path, output_gtf: Path, tmp_dir: Optional[Path], sort_memory: str) -> None:
	"""Sort GTF by contig and start coordinate using GNU sort."""
	if shutil.which("sort") is None:
		raise RuntimeError("'sort' command not found. Please install coreutils.")

	cmd = [
		"sort",
		"-k1,1",
		"-k4,4n",
		"-S",
		sort_memory,
		str(input_gtf),
	]
	if tmp_dir is not None:
		cmd.extend(["-T", str(tmp_dir)])
	_run_command(cmd, stdout_path=output_gtf)


def _compress_with_bgzip(input_gtf: Path, output_gz: Path, threads: int) -> None:
	"""Compress GTF with bgzip."""
	if shutil.which("bgzip") is None:
		raise RuntimeError("'bgzip' command not found. Install htslib or use --compress gzip.")
	cmd = ["bgzip", "-@", str(max(1, threads)), "-c", str(input_gtf)]
	_run_command(cmd, stdout_path=output_gz)


def _compress_with_gzip(input_gtf: Path, output_gz: Path) -> None:
	"""Compress GTF with Python gzip (no tabix random access)."""
	_ensure_parent(output_gz)
	with open(input_gtf, "rb") as in_fh, gzip.open(output_gz, "wb") as out_fh:
		shutil.copyfileobj(in_fh, out_fh)


def _build_tabix_index(input_gz: Path) -> Path:
	"""Build tabix index for compressed GTF/GFF."""
	if shutil.which("tabix") is None:
		raise RuntimeError("'tabix' command not found. Install htslib or use --no-tabix.")
	subprocess.run(["tabix", "-p", "gff", str(input_gz)], check=True)
	return input_gz.with_suffix(input_gz.suffix + ".tbi")


def _normalize_features(features: Iterable[str]) -> Set[str]:
	normalized = {item.strip() for item in features if item and item.strip()}
	if not normalized:
		raise ValueError("At least one feature must be provided via --keep-feature")
	return normalized


def preprocess_gtf(
	input_path: Path,
	output_path: Path,
	keep_features: Set[str],
	keep_header: bool,
	sort_records: bool,
	compress: str,
	tabix_index: bool,
	threads: int,
	sort_memory: str,
	tmp_dir: Optional[Path],
) -> None:
	"""Run GTF preprocessing pipeline."""
	if not input_path.exists():
		raise FileNotFoundError(f"Input GTF not found: {input_path}")

	_ensure_parent(output_path)

	with tempfile.TemporaryDirectory(prefix="gtf_preprocess_") as work_dir_str:
		work_dir = Path(work_dir_str)
		filtered_path = work_dir / "filtered.gtf"
		sorted_path = work_dir / "sorted.gtf"

		logging.info("Filtering features: %s", ", ".join(sorted(keep_features)))
		stats = _filter_features(
			input_gtf=input_path,
			output_gtf=filtered_path,
			keep_features=keep_features,
			keep_header=keep_header,
		)

		source_for_compression = filtered_path
		if sort_records:
			logging.info("Sorting filtered GTF by contig and start position")
			_sort_gtf(
				input_gtf=filtered_path,
				output_gtf=sorted_path,
				tmp_dir=tmp_dir,
				sort_memory=sort_memory,
			)
			source_for_compression = sorted_path

		if compress == "none":
			shutil.copy2(source_for_compression, output_path)
			logging.info("Wrote uncompressed output: %s", output_path)
		elif compress == "bgzip":
			_compress_with_bgzip(source_for_compression, output_path, threads=threads)
			logging.info("Wrote bgzip-compressed output: %s", output_path)
		elif compress == "gzip":
			_compress_with_gzip(source_for_compression, output_path)
			logging.info("Wrote gzip-compressed output: %s", output_path)
		else:
			raise ValueError(f"Unsupported compression mode: {compress}")

		if tabix_index:
			if compress != "bgzip":
				raise ValueError("tabix indexing requires --compress bgzip")
			index_path = _build_tabix_index(output_path)
			logging.info("Wrote tabix index: %s", index_path)

		logging.info(
			"Done. total=%d, headers=%d, malformed=%d, kept=%d, features=%d",
			stats.total_lines,
			stats.header_lines,
			stats.malformed_lines,
			stats.kept_feature_lines,
			stats.kept_features,
		)


def _build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description=(
			"Preprocess GTF by keeping selected feature types, sorting records, "
			"compressing output, and optionally building tabix index."
		)
	)
	parser.add_argument("-i", "--input", required=True, help="Input GTF(.gz) path")
	parser.add_argument("-o", "--output", required=True, help="Output path (.gtf, .gtf.gz)")
	parser.add_argument(
		"-k",
		"--keep-feature",
		required=True,
		nargs="+",
		help="Feature types to keep, e.g. gene exon transcript",
	)
	parser.add_argument("--drop-header", action="store_true", help="Drop comment/header lines")
	parser.add_argument("--no-sort", action="store_true", help="Skip sorting")
	parser.add_argument(
		"--compress",
		choices=["bgzip", "gzip", "none"],
		default="bgzip",
		help="Compression mode (default: bgzip)",
	)
	parser.add_argument("--threads", type=int, default=4, help="Threads for bgzip")
	parser.add_argument(
		"--sort-memory",
		default="2G",
		help="Memory for GNU sort, passed to sort -S (default: 2G)",
	)
	parser.add_argument("--tmp-dir", default="", help="Temporary directory for sorting")
	parser.add_argument("--no-tabix", action="store_true", help="Do not build tabix index")
	return parser


def main() -> None:
	parser = _build_parser()
	args = parser.parse_args()

	input_path = Path(args.input)
	output_path = Path(args.output)
	tmp_dir = Path(args.tmp_dir) if args.tmp_dir else None
	keep_features = _normalize_features(args.keep_feature)

	preprocess_gtf(
		input_path=input_path,
		output_path=output_path,
		keep_features=keep_features,
		keep_header=not args.drop_header,
		sort_records=not args.no_sort,
		compress=args.compress,
		tabix_index=not args.no_tabix,
		threads=args.threads,
		sort_memory=args.sort_memory,
		tmp_dir=tmp_dir,
	)


if __name__ == "__main__":
	main()
