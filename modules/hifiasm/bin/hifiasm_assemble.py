#!/usr/bin/env python3
"""Assemble HiFi BAM into primary-contig FASTA with hifiasm.

Streams ``samtools fastq`` straight into hifiasm's stdin (no intermediate
FASTQ on disk), then converts the unitig GFA to the primary-contig FASTA that
downstream rules consume.
"""
from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
from typing import List

logger = logging.getLogger("hifiasm_assemble")


def setup_logging() -> None:
    """Configure stderr logging for the engine entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def assemble_from_bam(bam: str, samtools: str, hifiasm: str,
                      out_prefix: str, threads: int) -> None:
    """Pipe samtools fastq output directly into hifiasm.

    Args:
        bam: input PacBio HiFi BAM (with colocated .bam.pbi index).
        samtools: samtools executable.
        hifiasm: hifiasm executable.
        out_prefix: hifiasm -o prefix; writes <prefix>.bp.p_ctg.gfa.
        threads: threads for both tools.

    Raises:
        RuntimeError: either tool exited non-zero or no reads were streamed.
    """
    logger.info("Executing: %s fastq -@ %d %s | %s -o %s -t %d /dev/stdin",
                samtools, threads, bam, hifiasm, out_prefix, threads)
    sam = subprocess.Popen(
        [samtools, "fastq", "-@", str(threads), bam],
        stdout=subprocess.PIPE,
    )
    assert sam.stdout is not None
    hi = subprocess.Popen(
        [hifiasm, "-o", out_prefix, "-t", str(threads), "/dev/stdin"],
        stdin=sam.stdout,
    )
    sam.stdout.close()
    hi_rc = hi.wait()
    sam_rc = sam.wait()
    if sam_rc != 0:
        raise RuntimeError(
            f"samtools fastq failed (exit {sam_rc}) for {bam}; "
            f"check the rule log for tool output.")
    if hi_rc != 0:
        raise RuntimeError(
            f"hifiasm failed (exit {hi_rc}) for {bam}; "
            f"check the rule log for tool output.")


def gfa_to_fasta(gfa: str, fasta: str) -> int:
    """Convert hifiasm unitig GFA (S lines) to FASTA; return contig count.

    Args:
        gfa: hifiasm assembly GFA (e.g. <prefix>.bp.p_ctg.gfa).
        fasta: output FASTA path.

    Returns:
        Number of contigs written.

    Raises:
        RuntimeError: the GFA holds no S records.
    """
    n = 0
    with open(gfa) as fin, open(fasta, "w") as fout:
        for line in fin:
            if not line.startswith("S"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            fout.write(f">{fields[1]}\n{fields[2]}\n")
            n += 1
    if n == 0:
        raise RuntimeError(
            f"No contigs found in {gfa}; hifiasm assembly failed or the GFA "
            f"format changed (expected 'S' records with name and sequence).")
    return n


def verify_output(path: str) -> None:
    """Assert a declared output exists and is non-empty."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise RuntimeError(
            f"Expected output missing or empty: {path}; "
            f"check the rule log for the failing stage.")


def main() -> None:
    """Parse arguments and run the assembly stages in order."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True,
                        help="Input PacBio HiFi BAM (with colocated .bam.pbi)")
    parser.add_argument("--out-prefix", dest="out_prefix", required=True,
                        help="hifiasm -o prefix; writes <prefix>.bp.p_ctg.gfa")
    parser.add_argument("--out-fasta", dest="out_fasta", required=True,
                        help="Primary-contig FASTA to write "
                             "(must equal <out-prefix>.bp.p_ctg.fa)")
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--samtools", default="samtools",
                        help="samtools executable (default: samtools)")
    parser.add_argument("--hifiasm", default="hifiasm",
                        help="hifiasm executable (default: hifiasm)")
    args = parser.parse_args()
    setup_logging()

    expected_fasta = args.out_prefix + ".bp.p_ctg.fa"
    if os.path.abspath(args.out_fasta) != os.path.abspath(expected_fasta):
        raise RuntimeError(
            f"--out-fasta ({args.out_fasta}) must equal "
            f"--out-prefix + '.bp.p_ctg.fa' ({expected_fasta}); "
            f"the module wires downstream rules to that exact path.")

    os.makedirs(os.path.dirname(args.out_fasta) or ".", exist_ok=True)
    assemble_from_bam(args.bam, args.samtools, args.hifiasm,
                      args.out_prefix, args.threads)
    n = gfa_to_fasta(args.out_prefix + ".bp.p_ctg.gfa", args.out_fasta)
    verify_output(args.out_fasta)
    logger.info("Wrote %d primary contigs to %s", n, args.out_fasta)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("%s", exc)
        sys.exit(1)
