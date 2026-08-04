#!/usr/bin/env python3
"""Create strict per-gene FASTQ, reference, and annotation inputs.

Reads are assigned once from a coordinate-sorted whole-genome BAM. Read names
that overlap more than one annotated gene are either excluded or duplicated,
according to ``--ambiguous``. Every emitted reference and GTF uses a gene-local
coordinate system oriented in the transcriptional direction.
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import IO, Iterable


def run(command: list[str], stdout: IO[str] | None = None) -> None:
    """Run a command and fail immediately when it returns a non-zero status."""
    subprocess.run(command, check=True, stdout=stdout)


def load_bed(path: Path) -> dict[str, list[str]]:
    """Load unique BED records keyed by gene ID."""
    records: dict[str, list[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if not line.strip() or line.startswith("#"):
                continue
            if len(fields) < 6:
                raise ValueError(f"BED line {line_number} has fewer than six fields")
            gene_id = fields[3]
            if gene_id in records and records[gene_id][:6] != fields[:6]:
                raise ValueError(f"Gene ID {gene_id!r} occurs at multiple BED loci")
            records[gene_id] = fields
    if not records:
        raise ValueError(f"No gene records found in {path}")
    return records


def normalize_read_name(name: str) -> str:
    """Normalize FASTQ/BAM mate suffixes without altering internal slashes."""
    if name.endswith("/1") or name.endswith("/2"):
        return name[:-2]
    return name


def read_overlap_assignments(
    overlap_path: Path,
    ambiguous: str,
) -> dict[str, list[str]]:
    """Return gene-to-read-name assignments from bedtools intersect output."""
    read_genes: dict[str, set[str]] = defaultdict(set)
    with overlap_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            # bedtools -bed emits six A fields followed by all B fields.
            if len(fields) < 12:
                raise ValueError(
                    f"Unexpected bedtools intersect record at line {line_number}: "
                    f"expected at least 12 fields, observed {len(fields)}"
                )
            read_name = normalize_read_name(fields[3])
            gene_id = fields[-4]
            read_genes[read_name].add(gene_id)

    assignments: dict[str, list[str]] = defaultdict(list)
    for read_name, genes in sorted(read_genes.items()):
        if ambiguous == "exclude" and len(genes) != 1:
            continue
        for gene_id in sorted(genes):
            assignments[gene_id].append(read_name)
    return assignments


def fastq_record_count(path: Path) -> int:
    """Count FASTQ records in a gzip-compressed file."""
    line_count = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_count, _ in enumerate(handle, start=1):
            pass
    if line_count % 4:
        raise ValueError(f"FASTQ has an incomplete record: {path}")
    return line_count // 4


def write_gene_reference(
    gene_id: str,
    fields: list[str],
    genome_fasta: Path,
    gene_dir: Path,
    flank: int,
    bedtools: str,
) -> tuple[Path, Path]:
    """Write a strand-oriented gene-local FASTA and matching GTF."""
    chrom, start_text, end_text, _, _, strand = fields[:6]
    start, end = int(start_text), int(end_text)
    local_start_genomic = max(0, start - flank)
    local_end_genomic = end + flank

    reference_bed = gene_dir / "reference.bed"
    reference_bed.write_text(
        "\t".join(
            [chrom, str(local_start_genomic), str(local_end_genomic), gene_id, ".", strand]
        )
        + "\n",
        encoding="utf-8",
    )

    raw_fasta = gene_dir / "reference.raw.fa"
    with raw_fasta.open("w", encoding="utf-8") as output_handle:
        run(
            [
                bedtools,
                "getfasta",
                "-fi",
                str(genome_fasta),
                "-bed",
                str(reference_bed),
                "-nameOnly",
                "-s",
            ],
            stdout=output_handle,
        )

    raw_lines = raw_fasta.read_text(encoding="utf-8").splitlines()
    if len(raw_lines) < 2:
        raise ValueError(f"Failed to extract reference sequence for {gene_id}")
    sequence = "".join(raw_lines[1:]).upper()
    reference_fasta = gene_dir / "reference.fa"
    reference_fasta.write_text(f">{gene_id}\n{sequence}\n", encoding="utf-8")
    raw_fasta.unlink()

    upstream_bases = start - local_start_genomic
    downstream_bases = local_end_genomic - end
    gene_length = end - start
    if strand == "+":
        local_gene_start = upstream_bases + 1
    elif strand == "-":
        local_gene_start = downstream_bases + 1
    else:
        raise ValueError(f"Unsupported strand {strand!r} for {gene_id}")
    local_gene_end = local_gene_start + gene_length - 1

    # bedtools getfasta -s already orients the sequence 5' to 3'. The local
    # annotation therefore uses '+' regardless of the genomic source strand.
    attributes = f'gene_id "{gene_id}"; gene_name "{gene_id}";'
    transcript_attributes = attributes + f' transcript_id "{gene_id}";'
    reference_gtf = gene_dir / "reference.gtf"
    reference_gtf.write_text(
        "\n".join(
            [
                f"{gene_id}\tOmics\tgene\t{local_gene_start}\t{local_gene_end}\t.\t+\t.\t{attributes}",
                f"{gene_id}\tOmics\ttranscript\t{local_gene_start}\t{local_gene_end}\t.\t+\t.\t{transcript_attributes}",
                f"{gene_id}\tOmics\texon\t{local_gene_start}\t{local_gene_end}\t.\t+\t.\t{transcript_attributes} exon_number \"1\";",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return reference_fasta, reference_gtf


def prepare_gene_fastq(
    read_names: Iterable[str],
    bam: Path,
    gene_dir: Path,
    samtools: str,
) -> tuple[Path, Path, str, int]:
    """Extract primary alignments and convert them to per-gene FASTQ files."""
    names_path = gene_dir / "read_names.txt"
    names_path.write_text("\n".join(read_names) + "\n", encoding="utf-8")
    assigned_bam = gene_dir / "assigned.primary.bam"
    name_sorted_bam = gene_dir / "assigned.primary.name_sorted.bam"
    fastq1 = gene_dir / "reads_1.fq.gz"
    fastq2 = gene_dir / "reads_2.fq.gz"

    run(
        [
            samtools,
            "view",
            "-N",
            str(names_path),
            "-F",
            "0x900",
            "-b",
            "-o",
            str(assigned_bam),
            str(bam),
        ]
    )
    run([samtools, "sort", "-n", "-o", str(name_sorted_bam), str(assigned_bam)])
    run(
        [
            samtools,
            "fastq",
            "-1",
            str(fastq1),
            "-2",
            str(fastq2),
            "-0",
            "/dev/null",
            "-s",
            "/dev/null",
            "-n",
            str(name_sorted_bam),
        ]
    )

    read1_count = fastq_record_count(fastq1)
    read2_count = fastq_record_count(fastq2)
    if read1_count == 0:
        raise ValueError(f"No FASTQ records were generated in {gene_dir}")
    if read2_count == 0:
        layout = "SE"
    elif read1_count == read2_count:
        layout = "PE"
    else:
        raise ValueError(
            f"Unbalanced paired FASTQ for {gene_dir.name}: R1={read1_count}, R2={read2_count}"
        )
    return fastq1, fastq2, layout, read1_count


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--genome-fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bedtools", default="bedtools")
    parser.add_argument("--ambiguous", choices=("exclude", "duplicate"), default="exclude")
    parser.add_argument("--flank", type=int, default=50)
    return parser.parse_args()


def main() -> None:
    """Create all per-gene inputs and an atomic manifest."""
    args = parse_args()
    if args.flank < 0:
        raise ValueError("--flank must be non-negative")
    if not args.bam.is_file() or not args.bed.is_file() or not args.genome_fasta.is_file():
        raise FileNotFoundError("BAM, BED, and genome FASTA inputs must exist")

    args.outdir.mkdir(parents=True, exist_ok=True)
    bed_records = load_bed(args.bed)
    overlap_path = args.outdir / "read_gene_overlaps.tsv"
    with overlap_path.open("w", encoding="utf-8") as output_handle:
        run(
            [
                args.bedtools,
                "intersect",
                "-abam",
                str(args.bam),
                "-b",
                str(args.bed),
                "-bed",
                "-wa",
                "-wb",
            ],
            stdout=output_handle,
        )

    assignments = read_overlap_assignments(overlap_path, args.ambiguous)
    manifest_tmp = args.outdir / "genes.tsv.tmp"
    manifest = args.outdir / "genes.tsv"
    emitted = 0
    with manifest_tmp.open("w", encoding="utf-8") as manifest_handle:
        manifest_handle.write(
            "gene_id\tfastq1\tfastq2\tfasta\tgtf\tlayout\tassigned_records\n"
        )
        for gene_id, read_names in sorted(assignments.items()):
            fields = bed_records.get(gene_id)
            if fields is None:
                raise ValueError(f"Assigned gene {gene_id!r} is absent from the BED records")
            gene_dir = args.outdir / gene_id
            if gene_dir.exists():
                shutil.rmtree(gene_dir)
            gene_dir.mkdir(parents=True)
            fastq1, fastq2, layout, record_count = prepare_gene_fastq(
                read_names, args.bam, gene_dir, args.samtools
            )
            reference_fasta, reference_gtf = write_gene_reference(
                gene_id,
                fields,
                args.genome_fasta,
                gene_dir,
                args.flank,
                args.bedtools,
            )
            manifest_handle.write(
                "\t".join(
                    [
                        gene_id,
                        str(fastq1),
                        str(fastq2),
                        str(reference_fasta),
                        str(reference_gtf),
                        layout,
                        str(record_count),
                    ]
                )
                + "\n"
            )
            emitted += 1
    if emitted == 0:
        raise ValueError("No unambiguous gene-specific read assignments were produced")
    os.replace(manifest_tmp, manifest)
    print(f"Prepared {emitted} gene-specific inputs in {args.outdir}")


if __name__ == "__main__":
    main()
