#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pyfaidx import Fasta

GENOME = None

def parse_fasta(path, keep_names=None):
    """
    Parse a FASTA file into a dictionary of sequences.

    Parameters
    ----------
    path : str
        Path to the FASTA file.
    keep_names : set or None, optional
        If provided, only sequences with names in this set are kept.

    Returns
    -------
    dict
        Mapping from sequence name to uppercase sequence string.
    """
    seqs = {}
    name = None
    chunks = []
    keep = keep_names if keep_names else None

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name and (keep is None or name in keep):
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if keep is None or name in keep:
                    chunks.append(line)
        if name and (keep is None or name in keep):
            seqs[name] = "".join(chunks).upper()
    return seqs

def revcomp(seq):
    """
    Compute the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    str
        Reverse complemented sequence.
    """
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def parse_gtf_attrs(attr_str):
    """
    Parse the attribute column of a GTF line into a dict.

    Parameters
    ----------
    attr_str : str
        Raw attribute string from GTF (9th column).

    Returns
    -------
    dict
        Parsed attributes as key-value pairs.
    """
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            k, v = part.split(" ", 1)
            v = v.strip().strip('"')
            attrs[k] = v
    return attrs

def is_rrna_feature(feature, attrs):
    """
    Determine whether a GTF record describes an rRNA feature.

    Parameters
    ----------
    feature : str
        Feature type (3rd column of GTF).
    attrs : dict
        Parsed attribute dictionary.

    Returns
    -------
    bool
        True if the record is considered rRNA, otherwise False.
    """
    if feature.lower() == "rrna":
        return True
    for key in ("gene_biotype", "gene_type", "transcript_biotype", "transcript_type", "biotype", "type"):
        if attrs.get(key, "").lower() == "rrna":
            return True
    for key in ("gene_name", "transcript_name", "gene_id", "transcript_id"):
        if "rrna" in attrs.get(key, "").lower():
            return True
    return False

def init_fasta(fasta_path):
    """
    Initialize a global FASTA index for worker processes.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    None
        Sets a global indexed FASTA handle.
    """
    global GENOME
    GENOME = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)

def extract_rrna_record(task):
    """
    Extract a single rRNA record sequence.

    Parameters
    ----------
    task : tuple
        Tuple of (rid, feats, meta).

    Returns
    -------
    tuple
        (rid, header, sequence, error_message). error_message is None on success.
    """
    rid, feats, meta = task
    exons = [x for x in feats if x[4].lower() == "exon"]
    use = exons if exons else feats

    chrom = use[0][0]
    strand = use[0][3]

    if chrom not in GENOME:
        return rid, None, None, f"Chromosome {chrom} not in genome FASTA; skip {rid}"

    use_sorted = sorted(use, key=lambda x: x[1])
    seq = "".join(GENOME[chrom][s-1:e] for _, s, e, _, _ in use_sorted)
    if strand == "-":
        seq = revcomp(seq)

    m = meta.get(rid, {})
    header = rid
    if m.get("gene_name"):
        header += f"|{m['gene_name']}"
    if m.get("transcript_name"):
        header += f"|{m['transcript_name']}"
    if m.get("gene_type"):
        header += f"|{m['gene_type']}"
    header += f"|{chrom}:{use_sorted[0][1]}-{use_sorted[-1][2]}({strand})"

    return rid, header, seq, None

def main():
    """
    Extract rRNA sequences from genome FASTA and GTF and write to FASTA.

    Returns
    -------
    None
        Writes output FASTA and exits with non-zero status on errors.
    """
    ap = argparse.ArgumentParser(description="Extract rRNA sequences from genome FASTA and GTF.")
    ap.add_argument("-f", "--fasta", required=True, help="Genome FASTA")
    ap.add_argument("-t", "--gtf", required=True, help="GTF annotation")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA for rRNA")
    ap.add_argument("-p", "--threads", type=int, default=max(1, os.cpu_count() or 1),
                    help="Number of parallel worker processes")
    args = ap.parse_args()

    intervals = defaultdict(list)
    meta = {}

    # First pass: parse GTF and collect rRNA intervals
    with open(args.gtf, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attr_str = parts

            if "rrna" not in feature.lower() and "rrna" not in attr_str.lower():
                continue

            attrs = parse_gtf_attrs(attr_str)
            if not is_rrna_feature(feature, attrs):
                continue

            rid = attrs.get("transcript_id") or attrs.get("gene_id")
            if not rid:
                continue

            start = int(start)
            end = int(end)
            intervals[rid].append((chrom, start, end, strand, feature))
            if rid not in meta:
                meta[rid] = {
                    "gene_id": attrs.get("gene_id", ""),
                    "transcript_id": attrs.get("transcript_id", ""),
                    "gene_name": attrs.get("gene_name", ""),
                    "transcript_name": attrs.get("transcript_name", ""),
                    "gene_type": attrs.get("gene_biotype", "") or attrs.get("gene_type", ""),
                }

    if not intervals:
        sys.stderr.write("No rRNA features found in GTF.\n")
        sys.exit(1)

    tasks = [(rid, feats, meta) for rid, feats in intervals.items()]

    with open(args.output, "w") as out, ProcessPoolExecutor(
        max_workers=args.threads,
        initializer=init_fasta,
        initargs=(args.fasta,)
    ) as ex:
        futures = [ex.submit(extract_rrna_record, t) for t in tasks]
        for fut in as_completed(futures):
            rid, header, seq, err = fut.result()
            if err:
                sys.stderr.write(err + "\n")
                continue
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

if __name__ == "__main__":
    main()