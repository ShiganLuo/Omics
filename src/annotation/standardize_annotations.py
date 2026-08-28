#!/usr/bin/env python3
"""Standardize gene GTF and TE BED files for consistent downstream analysis.

Handles:
1. Chromosome naming: normalize to Ensembl convention (no 'chr' prefix)
2. Mitochondrial genes: add 'MT-' prefix to gene_name on chrMT(only Human and Mouse is that standard)
3. TE BED: strip 'chr' prefix to match GTF convention

Usage:
    python standardize_annotations.py --gtf <gene.gtf> [--te <te.bed>]

Examples:
    python standardize_annotations.py --gtf Macaca_mulatta.Mmul_10.116.gtf --te rheMac10_rmsk_TE.bed
    python standardize_annotations.py --gtf Homo_sapiens.GRCh38.116.gtf
"""

import argparse
import gzip
import os
import sys
import tempfile
from typing import Dict, List, Set, Tuple

# Mitochondrial chromosome names to recognize
MT_CHROMS: Set[str] = {'MT', 'chrM', 'chrMT', 'M'}


def is_mitochondrial_chrom(chrom: str) -> bool:
    """Check if chromosome is mitochondrial."""
    return chrom in MT_CHROMS


def normalize_chrom(chrom: str) -> str:
    """Normalize chromosome name to Ensembl convention (no 'chr' prefix).

    Special handling for mitochondrial chromosome: always returns 'MT'.
    """
    # Strip 'chr' prefix
    if chrom.startswith('chr'):
        chrom = chrom[3:]

    # Normalize mitochondrial chromosome names
    if chrom in ('M', 'MT', 'chrM', 'chrMT'):
        return 'MT'

    return chrom


def standardize_gtf(
    gtf_path: str,
    dry_run: bool = False,
) -> Dict[str, int]:
    """Standardize a gene GTF file in-place.

    Operations:
    1. Normalize chromosome names (strip 'chr', M/chrM → MT)
    2. Add 'MT-' prefix to mitochondrial gene_names (if not already prefixed)

    Parameters
    ----------
    gtf_path : str
        Path to the GTF file.
    dry_run : bool
        If True, only report what would be changed.

    Returns
    -------
    Dict[str, int]
        Statistics: total_lines, chrom_fixed, mt_genes_prefixed.
    """
    stats = {'total_lines': 0, 'chrom_fixed': 0, 'mt_genes_prefixed': 0}

    if not os.path.isfile(gtf_path):
        print(f"  ERROR: GTF not found: {gtf_path}")
        return stats

    print(f"Standardizing GTF: {os.path.basename(gtf_path)}")
    is_gzipped = gtf_path.endswith('.gz')

    # Read all lines
    opener = gzip.open if is_gzipped else open
    with opener(gtf_path, 'rt', encoding='utf-8') as f:
        lines = f.readlines()

    new_lines: List[str] = []
    for line in lines:
        stats['total_lines'] += 1

        # Skip comments and empty lines
        if not line or line.startswith('#'):
            new_lines.append(line)
            continue

        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            new_lines.append(line)
            continue

        # Column 1: chromosome
        orig_chrom = parts[0]
        new_chrom = normalize_chrom(orig_chrom)
        if new_chrom != orig_chrom:
            stats['chrom_fixed'] += 1
            parts[0] = new_chrom

        # Column 9: attributes - fix mitochondrial gene_name
        if is_mitochondrial_chrom(new_chrom):
            attrs = parts[8]
            if 'gene_name' in attrs:
                # Extract gene_name value
                try:
                    before, rest = attrs.split('gene_name "', 1)
                    gene_name, after = rest.split('"', 1)

                    # Add MT- prefix if not already present
                    if not gene_name.startswith('MT-'):
                        new_gene_name = f'MT-{gene_name}'
                        attrs = f'{before}gene_name "{new_gene_name}"{after}'
                        parts[8] = attrs
                        stats['mt_genes_prefixed'] += 1
                except ValueError:
                    pass  # Malformed attribute, skip

        new_lines.append('\t'.join(parts) + '\n')

    if dry_run:
        print(f"  DRY RUN: would fix {stats['chrom_fixed']} chromosomes, "
              f"{stats['mt_genes_prefixed']} MT gene names")
        return stats

    # Write back
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(gtf_path), suffix='.tmp')
    try:
        if is_gzipped:
            with gzip.open(tmp, 'wt', encoding='utf-8') as f:
                f.writelines(new_lines)
        else:
            with os.fdopen(fd, 'w', encoding='utf-8') as f:
                f.writelines(new_lines)
        os.replace(tmp, gtf_path)
    except BaseException:
        os.unlink(tmp)
        raise

    print(f"  Chromosomes normalized: {stats['chrom_fixed']}")
    print(f"  MT genes prefixed: {stats['mt_genes_prefixed']}")
    return stats


def standardize_bed(
    bed_path: str,
    dry_run: bool = False,
) -> Dict[str, int]:
    """Standardize a TE BED file in-place.

    Operations:
    1. Normalize chromosome names (strip 'chr', M/chrM → MT)

    Parameters
    ----------
    bed_path : str
        Path to the BED file.
    dry_run : bool
        If True, only report what would be changed.

    Returns
    -------
    Dict[str, int]
        Statistics: total_lines, chrom_fixed.
    """
    stats = {'total_lines': 0, 'chrom_fixed': 0}

    if not os.path.isfile(bed_path):
        print(f"  ERROR: BED not found: {bed_path}")
        return stats

    print(f"Standardizing BED: {os.path.basename(bed_path)}")

    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(bed_path), suffix='.tmp')
    try:
        with open(bed_path, 'r', encoding='utf-8') as fin, \
             os.fdopen(fd, 'w', encoding='utf-8') as fout:
            for line in fin:
                stats['total_lines'] += 1

                if not line or line.startswith('#'):
                    fout.write(line)
                    continue

                parts = line.split('\t', 1)
                orig_chrom = parts[0]
                new_chrom = normalize_chrom(orig_chrom)

                if new_chrom != orig_chrom:
                    stats['chrom_fixed'] += 1
                    parts[0] = new_chrom

                fout.write('\t'.join(parts))

        if dry_run:
            os.unlink(tmp)
            print(f"  DRY RUN: would fix {stats['chrom_fixed']} chromosomes")
            return stats

        os.replace(tmp, bed_path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise

    print(f"  Chromosomes normalized: {stats['chrom_fixed']}")
    return stats


def find_annotation_files(base_dir: str) -> List[Tuple[str, str, str]]:
    """Find GTF and BED files in the reference directory structure.

    Expected layout:
        base_dir/<species>/<source>/<build>/
            *.gtf or *.gtf.gz
            *.bed

    Returns
    -------
    List[Tuple[str, str, str]]
        List of (species_dir, gtf_path, bed_path) tuples.
    """
    results = []
    for species_dir in sorted(os.listdir(base_dir)):
        species_path = os.path.join(base_dir, species_dir)
        if not os.path.isdir(species_path):
            continue
        for source_dir in os.listdir(species_path):
            source_path = os.path.join(species_path, source_dir)
            if not os.path.isdir(source_path):
                continue
            for build_dir in os.listdir(source_path):
                build_path = os.path.join(source_path, build_dir)
                if not os.path.isdir(build_path):
                    continue

                gtf_path = None
                bed_path = None

                for f in os.listdir(build_path):
                    if f.endswith('.gtf') or f.endswith('.gtf.gz'):
                        gtf_path = os.path.join(build_path, f)
                    elif f.endswith('.bed') and 'rmsk' in f.lower():
                        bed_path = os.path.join(build_path, f)

                if gtf_path:
                    results.append((species_dir, gtf_path, bed_path))

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Standardize gene GTF and TE BED annotation files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--gtf', type=str, nargs='+',
        help='Gene GTF file(s) to standardize.',
    )
    parser.add_argument(
        '--te', type=str, nargs='+',
        help='TE BED file(s) to standardize.',
    )
    parser.add_argument(
        '--scan', type=str, default=None,
        help='Scan a directory for annotation files and report status.',
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Only report what would be changed, do not modify files.',
    )

    args = parser.parse_args()

    if args.scan:
        # Scan mode: find and report annotation files
        print(f"Scanning: {args.scan}")
        print("=" * 60)
        files = find_annotation_files(args.scan)
        for species, gtf, bed in files:
            print(f"\n{species}:")
            print(f"  GTF: {gtf}")
            if bed:
                print(f"  BED: {bed}")
            else:
                print(f"  BED: (not found)")
        print(f"\nTotal: {len(files)} species")
        return

    if not args.gtf and not args.te:
        parser.print_help()
        sys.exit(1)

    total_stats = {'gtf': 0, 'bed': 0, 'chrom_fixed': 0, 'mt_prefixed': 0}

    # Process GTF files
    if args.gtf:
        for gtf_path in args.gtf:
            stats = standardize_gtf(gtf_path, args.dry_run)
            total_stats['gtf'] += 1
            total_stats['chrom_fixed'] += stats['chrom_fixed']
            total_stats['mt_prefixed'] += stats['mt_genes_prefixed']

    # Process BED files
    if args.te:
        for bed_path in args.te:
            stats = standardize_bed(bed_path, args.dry_run)
            total_stats['bed'] += 1
            total_stats['chrom_fixed'] += stats['chrom_fixed']

    print("\n" + "=" * 60)
    print(f"Summary:")
    print(f"  GTF files processed: {total_stats['gtf']}")
    print(f"  BED files processed: {total_stats['bed']}")
    print(f"  Chromosomes normalized: {total_stats['chrom_fixed']}")
    print(f"  MT genes prefixed: {total_stats['mt_prefixed']}")


if __name__ == '__main__':
    main()
