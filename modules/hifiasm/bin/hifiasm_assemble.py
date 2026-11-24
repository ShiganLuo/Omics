#!/usr/bin/env python3
"""Assemble HiFi BAM into primary-contig FASTA with hifiasm.

Streams ``samtools fastq`` straight into hifiasm's stdin (no intermediate
FASTQ on disk), then converts the unitig GFA to the primary-contig FASTA that
downstream rules consume. Emits stage-boundary log lines with elapsed times
so the rule log stays informative while hifiasm itself is quiet: hifiasm's
own [M::...] stderr is appended to the same log, but prints in bulk near the
end of the run.

The hifiasm GFA is kept and reused across reruns when a fingerprint stamp
matches the input BAM, so a failure in the cheap GFA->FASTA stage never
forces a reassembly. hifiasm 0.25 writes raw non-UTF8 bytes into GFA 'A'
(arc) records, so the GFA is parsed in binary mode.
"""
from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import subprocess
import sys
import time
from typing import List, NamedTuple

logger = logging.getLogger("hifiasm_assemble")


class FastaStats(NamedTuple):
    """Summary of the converted primary-contig FASTA."""
    n_contigs: int
    total_bp: int
    n50_bp: int


def setup_logging() -> None:
    """Configure stderr logging for the engine entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def tool_version(exe: str, args: List[str]) -> str:
    """Return the first output line of ``exe args`` for provenance logging.

    Args:
        exe: tool executable.
        args: version-query arguments (e.g. ['--version']).

    Returns:
        First stdout/stderr line, or 'unknown' when the probe fails.
    """
    try:
        out = subprocess.run([exe] + args, capture_output=True, text=True,
                             timeout=60)
        lines = (out.stdout or out.stderr).strip().splitlines()
        return lines[0] if lines else "unknown"
    except Exception:
        return "unknown"


def bam_fingerprint(bam: str) -> dict:
    """Return a path+size+mtime fingerprint of the input BAM."""
    st = os.stat(bam)
    return {"bam": os.path.abspath(bam), "size": st.st_size,
            "mtime": st.st_mtime}


def stamp_path(out_prefix: str) -> str:
    """Return the fingerprint stamp path for an hifiasm output prefix."""
    return out_prefix + ".hifiasm.stamp"


def can_reuse_gfa(gfa: str, out_prefix: str, bam: str) -> bool:
    """Decide whether an existing hifiasm GFA can be reused.

    Reuse is granted on a matching fingerprint stamp; an existing GFA
    without a stamp is reused with a warning (a crashed run leaves no
    stamp but its GFA is complete once hifiasm printed its final summary).

    Args:
        gfa: hifiasm primary-contig GFA path.
        out_prefix: hifiasm -o prefix (locates the stamp).
        bam: input BAM whose fingerprint must match.

    Returns:
        True when the GFA is complete and safe to skip reassembly for.
    """
    if not os.path.exists(gfa) or os.path.getsize(gfa) == 0:
        return False
    stamp = stamp_path(out_prefix)
    if os.path.exists(stamp):
        try:
            with open(stamp) as f:
                recorded = json.load(f)
        except (OSError, ValueError):
            recorded = None
        if recorded == bam_fingerprint(bam):
            logger.info("Stage 1/2 skipped: reusing existing GFA %s "
                        "(fingerprint stamp matches %s)", gfa, bam)
            return True
        logger.info("Existing GFA %s has a stale fingerprint stamp; "
                    "reassembling", gfa)
        return False
    logger.warning(
        "Reusing existing GFA %s without a fingerprint stamp (likely left "
        "by a crashed run). Delete %s (and %s) to force a clean reassembly.",
        gfa, gfa, out_prefix + ".bp.hap1.p_ctg.gfa")
    return True


def write_stamp(out_prefix: str, bam: str) -> None:
    """Record the input BAM fingerprint after a successful assembly."""
    with open(stamp_path(out_prefix), "w") as f:
        json.dump(bam_fingerprint(bam), f)


def assemble_from_bam(bam: str, samtools: str, hifiasm: str,
                      out_prefix: str, threads: int) -> float:
    """Pipe samtools fastq output directly into hifiasm.

    Args:
        bam: input PacBio HiFi BAM (with colocated .bam.pbi index).
        samtools: samtools executable.
        hifiasm: hifiasm executable.
        out_prefix: hifiasm -o prefix; writes <prefix>.bp.p_ctg.gfa.
        threads: threads for both tools.

    Returns:
        Elapsed seconds for the streaming assembly.

    Raises:
        RuntimeError: either tool exited non-zero.
    """
    logger.info(
        "Stage 1/2: streaming %s into hifiasm (hifiasm is quiet during read "
        "loading/correction; its own [M::...] output appears below in bulk "
        "near the end)", bam)
    start = time.monotonic()
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
    if hi_rc != 0:
        raise RuntimeError(
            f"hifiasm failed (exit {hi_rc}) for {bam}; "
            f"check the rule log for tool output.")
    if sam_rc != 0 and sam_rc not in (-13, 141):
        raise RuntimeError(
            f"samtools fastq failed (exit {sam_rc}) for {bam}; "
            f"check the rule log for tool output.")
    if sam_rc != 0:
        logger.warning(
            "samtools fastq terminated by SIGPIPE (hifiasm stopped reading "
            "early, e.g. when resuming from its cached overlaps); assembly "
            "completed successfully")
    elapsed = time.monotonic() - start
    logger.info("Stage 1/2 finished: hifiasm assembly completed in %.1f s",
                elapsed)
    return elapsed


def gfa_to_fasta(gfa: str, fasta: str) -> FastaStats:
    """Convert hifiasm unitig GFA (S lines) to FASTA and summarize it.

    The GFA is read in binary mode: hifiasm 0.25 emits raw non-UTF8 bytes
    in 'A' (arc) records, so only the ASCII name/sequence fields of 'S'
    records are decoded.

    Args:
        gfa: hifiasm assembly GFA (e.g. <prefix>.bp.p_ctg.gfa).
        fasta: output FASTA path.

    Returns:
        Contig count, total bp and N50 of the written FASTA.

    Raises:
        RuntimeError: the GFA holds no S records.
    """
    lengths: List[int] = []
    with open(gfa, "rb") as fin, open(fasta, "wb") as fout:
        for line in fin:
            if not line.startswith(b"S"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            seq = fields[2]
            fout.write(b">" + fields[1] + b"\n" + seq + b"\n")
            lengths.append(len(seq))
    if not lengths:
        raise RuntimeError(
            f"No contigs found in {gfa}; hifiasm assembly failed or the GFA "
            f"format changed (expected 'S' records with name and sequence).")
    total = sum(lengths)
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= total / 2:
            n50 = length
            break
    return FastaStats(len(lengths), total, n50)


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
    parser.add_argument("--force-reassemble", dest="force_reassemble",
                        action="store_true",
                        help="Ignore any existing GFA and rerun hifiasm")
    args = parser.parse_args()
    setup_logging()

    expected_fasta = args.out_prefix + ".bp.p_ctg.fa"
    if os.path.abspath(args.out_fasta) != os.path.abspath(expected_fasta):
        raise RuntimeError(
            f"--out-fasta ({args.out_fasta}) must equal "
            f"--out-prefix + '.bp.p_ctg.fa' ({expected_fasta}); "
            f"the module wires downstream rules to that exact path.")

    logger.info("Versions: hifiasm %s; samtools %s",
                tool_version(args.hifiasm, ["--version"]),
                tool_version(args.samtools, ["--version"]))
    os.makedirs(os.path.dirname(args.out_fasta) or ".", exist_ok=True)
    gfa = args.out_prefix + ".bp.p_ctg.gfa"
    if args.force_reassemble:
        stale = [p for p in glob.glob(args.out_prefix + ".*")
                 if os.path.isfile(p)]
        for path in stale:
            os.remove(path)
        logger.info("--force-reassemble: removed %d cached hifiasm artifact(s)",
                    len(stale))
    if args.force_reassemble or not can_reuse_gfa(gfa, args.out_prefix,
                                                  args.bam):
        assemble_from_bam(args.bam, args.samtools, args.hifiasm,
                          args.out_prefix, args.threads)
        write_stamp(args.out_prefix, args.bam)
    logger.info("Stage 2/2: converting GFA to FASTA")
    stats = gfa_to_fasta(gfa, args.out_fasta)
    verify_output(args.out_fasta)
    logger.info(
        "Stage 2/2 finished: wrote %d primary contigs to %s "
        "(total %d bp, N50 %d bp)",
        stats.n_contigs, args.out_fasta, stats.total_bp, stats.n50_bp)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("%s", exc)
        sys.exit(1)
