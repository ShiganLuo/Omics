#!/usr/bin/env python3
"""Canonical three-pass STAR alignment for small non-coding RNA reads.

Directly executes the alignment pipeline via subprocess:

  1. Whole-genome end-to-end alignment (pass 1)
  2. Small-RNA reference end-to-end alignment (pass 2)
  3. Genome local alignment of pass-2 mapped and unmapped reads (pass 3a/3b)
  4. Merge final BAM and index

A command log is written for reproducibility.
"""

from __future__ import annotations

import argparse
import gzip
import os
import shlex
import subprocess
import time
from dataclasses import dataclass, field


# ── Per-pass STAR option assembly ────────────────────────────────────────

_PASS_PARAMS = [
    ("out_filter_multimap_nmax", "--outFilterMultimapNmax", int),
    ("out_filter_multimap_score_range", "--outFilterMultimapScoreRange", int),
    ("out_filter_mismatch_nover_lmax", "--outFilterMismatchNoverLmax", float),
    ("align_intron_min", "--alignIntronMin", int),
    ("align_mates_gap_max", "--alignMatesGapMax", int),
    ("align_ends_type", "--alignEndsType", str),
]

_PASS2_PARAMS = [
    ("out_filter_mismatch_nover_read_lmax", "--outFilterMismatchNoverReadLmax", float),
    ("clip5p_nbases", "--clip5pNbases", str),
    ("clip3p_nbases", "--clip3pNbases", str),
]


def _star_pass_options(pass_name: str, paired: bool, args: argparse.Namespace,
                       force_end_to_end: bool, hard_clip_5p: int) -> list[str]:
    """Assemble STAR flag tokens for one pass from typed CLI arguments."""
    opts: dict[str, str] = {}

    for attr, flag, cast in _PASS_PARAMS:
        val = getattr(args, f"{pass_name}_{attr}")
        if val is not None:
            opts[flag] = str(cast(val))

    if pass_name == "pass2":
        for attr, flag, cast in _PASS2_PARAMS:
            val = getattr(args, f"pass2_{attr}")
            if val is not None:
                opts[flag] = str(cast(val))

    opts.pop("--outReadsUnmapped", None)

    # clip5p/clip3p defaults depend on paired/SE; apply when user did not
    # explicitly override (argparse default is None).
    if pass_name == "pass2":
        if opts.get("--clip5pNbases") is None:
            opts["--clip5pNbases"] = "20 0" if paired else "20"
        if opts.get("--clip3pNbases") is None:
            opts["--clip3pNbases"] = "0 20" if paired else "0"

    if force_end_to_end:
        opts["--alignEndsType"] = "EndToEnd"
    if hard_clip_5p:
        opts["--clip5pNbases"] = f"{hard_clip_5p} 0" if paired else str(hard_clip_5p)

    tokens: list[str] = []
    for flag, value in opts.items():
        tokens.append(flag)
        tokens.extend(value.split())
    return tokens


def _samtools_fastq_cmds(samtools: str, threads: int, source_bam: str,
                         read1: str, read2: str, paired: bool) -> list[list[str]]:
    """Return samtools fastq commands as a list of arg-lists."""
    if paired:
        return [[samtools, "fastq", "-@", str(threads), "-n",
                 "-1", read1, "-2", read2,
                 "-0", "/dev/null", "-s", "/dev/null", source_bam]]
    return [[samtools, "fastq", "-@", str(threads), "-n",
             "-0", read1, "-s", "/dev/null", source_bam]]


# ── Command plan ─────────────────────────────────────────────────────────

@dataclass
class Step:
    """One executable step in the pipeline.

    ``cmd`` is the primary command.  If ``pipe_to`` is set, stdout of ``cmd``
    is piped into ``pipe_to``.  If ``stdout`` is set, stdout of ``cmd`` (or
    ``pipe_to``) is written to that path instead of the log.
    """
    cmd: list[str]
    desc: str = ""
    pipe_to: list[str] | None = None
    stdout: str | None = None      # redirect final stdout to this file
    skip_log: bool = False         # if True, don't log stdout/stderr to log


def build_steps(args: argparse.Namespace, paired: bool) -> list[Step]:
    """Build the ordered list of execution steps."""
    star = args.star
    samtools = args.samtools
    bedtools = args.bedtools
    threads = args.threads

    pass1_prefix = os.path.dirname(args.pass1_bam) + "/star."
    pass2_prefix = os.path.dirname(args.pass2_bam) + "/star."
    pass3a_prefix = os.path.dirname(args.pass3a_bam) + "/star."
    pass3b_prefix = os.path.dirname(args.pass3b_bam) + "/star."
    pass1_name_bam = os.path.dirname(args.pass1_bam) + "/name_sorted.bam"
    pass2_name_bam = os.path.dirname(args.pass2_bam) + "/name_sorted.bam"
    pass2_mapped_fq1 = os.path.dirname(args.pass2_bam) + "/mapped_1.fq.gz"
    pass2_mapped_fq2 = os.path.dirname(args.pass2_bam) + "/mapped_2.fq.gz"
    pass3a_raw_bam = os.path.dirname(args.pass3a_bam) + "/raw.bam"

    steps: list[Step] = []

    # ── Pass 1: whole-genome end-to-end alignment ───────────────────
    steps.append(Step(
        cmd=[star, "--runThreadN", str(threads),
             "--genomeDir", args.genome_index,
             "--readFilesIn", *args.fastq,
             "--readFilesCommand", "zcat",
             *_star_pass_options("pass1", paired, args,
                                 args.force_end_to_end, args.hard_clip_5p),
             "--outSAMtype", "BAM", "SortedByCoordinate",
             "--outFileNamePrefix", pass1_prefix],
        desc="pass1: whole-genome end-to-end alignment",
    ))

    # Extract small-RNA reads from pass 1 (bedtools | samtools).
    steps.append(Step(
        cmd=[bedtools, "intersect", "-abam",
             pass1_prefix + "Aligned.sortedByCoord.out.bam",
             "-b", args.smallrna_bed, "-u"],
        pipe_to=[samtools, "sort", "-@", str(threads),
                 "-o", args.pass1_bam, "-"],
        desc="pass1: extract small-RNA reads",
    ))

    # Index, name-sort, and extract FASTQ from pass 1.
    steps.append(Step(
        cmd=[samtools, "index", "-@", str(threads), args.pass1_bam],
        desc="pass1: index BAM",
    ))
    steps.append(Step(
        cmd=[samtools, "sort", "-n", "-@", str(threads),
             "-o", pass1_name_bam, args.pass1_bam],
        desc="pass1: name-sort BAM",
    ))
    for fq_cmd in _samtools_fastq_cmds(samtools, threads, pass1_name_bam,
                                       args.pass1_fq1, args.pass1_fq2, paired):
        steps.append(Step(cmd=fq_cmd, desc="pass1: extract FASTQ"))

    # ── Pass 2: small-RNA reference end-to-end alignment ────────────
    pass1_reads = [args.pass1_fq1]
    if paired:
        pass1_reads.append(args.pass1_fq2)
    steps.append(Step(
        cmd=[star, "--runThreadN", str(threads),
             "--genomeDir", args.smallrna_index,
             "--readFilesIn", *pass1_reads,
             "--readFilesCommand", "zcat",
             *_star_pass_options("pass2", paired, args,
                                 args.force_end_to_end, args.hard_clip_5p),
             "--outReadsUnmapped", "Fastx",
             "--outSAMtype", "BAM", "SortedByCoordinate",
             "--outFileNamePrefix", pass2_prefix],
        desc="pass2: small-RNA reference end-to-end alignment",
    ))
    steps.append(Step(
        cmd=["mv", pass2_prefix + "Aligned.sortedByCoord.out.bam", args.pass2_bam],
        desc="pass2: move BAM",
    ))
    steps.append(Step(
        cmd=[samtools, "index", "-@", str(threads), args.pass2_bam],
        desc="pass2: index BAM",
    ))
    # gzip unmapped reads
    steps.append(Step(
        cmd=["gzip", "-c", pass2_prefix + "Unmapped.out.mate1"],
        stdout=args.pass2_unmapped1,
        desc="pass2: compress unmapped mate1",
    ))
    if paired:
        steps.append(Step(
            cmd=["gzip", "-c", pass2_prefix + "Unmapped.out.mate2"],
            stdout=args.pass2_unmapped2,
            desc="pass2: compress unmapped mate2",
        ))
    else:
        # SE: create empty placeholder for mate2
        steps.append(Step(
            cmd=["python", "-c",
                 f"import gzip; gzip.open({args.pass2_unmapped2!r}, \"wb\").close()"],
            desc="pass2: create empty mate2 placeholder",
            skip_log=True,
        ))

    # Name-sort pass 2 BAM and extract mapped reads as FASTQ.
    steps.append(Step(
        cmd=[samtools, "sort", "-n", "-@", str(threads),
             "-o", pass2_name_bam, args.pass2_bam],
        desc="pass2: name-sort BAM",
    ))
    for fq_cmd in _samtools_fastq_cmds(samtools, threads, pass2_name_bam,
                                       pass2_mapped_fq1, pass2_mapped_fq2, paired):
        steps.append(Step(cmd=fq_cmd, desc="pass2: extract mapped FASTQ"))

    # ── Pass 3a: local alignment of pass 2 mapped reads ─────────────
    pass2_mapped_reads = [pass2_mapped_fq1]
    if paired:
        pass2_mapped_reads.append(pass2_mapped_fq2)
    steps.append(Step(
        cmd=[star, "--runThreadN", str(threads),
             "--genomeDir", args.genome_index,
             "--readFilesIn", *pass2_mapped_reads,
             "--readFilesCommand", "zcat",
             *_star_pass_options("pass3", paired, args,
                                 args.force_end_to_end, args.hard_clip_5p),
             "--outSAMtype", "BAM", "SortedByCoordinate",
             "--outFileNamePrefix", pass3a_prefix],
        desc="pass3a: local alignment of mapped reads",
    ))
    steps.append(Step(
        cmd=["mv", pass3a_prefix + "Aligned.sortedByCoord.out.bam", pass3a_raw_bam],
        desc="pass3a: move raw BAM",
    ))
    steps.append(Step(
        cmd=[bedtools, "intersect", "-abam", pass3a_raw_bam,
             "-b", args.smallrna_bed, "-u"],
        pipe_to=[samtools, "sort", "-@", str(threads),
                 "-o", args.pass3a_bam, "-"],
        desc="pass3a: extract and sort small-RNA reads",
    ))
    steps.append(Step(
        cmd=[samtools, "index", "-@", str(threads), args.pass3a_bam],
        desc="pass3a: index BAM",
    ))

    # ── Pass 3b: local alignment of pass 2 unmapped reads ───────────
    pass2_unmapped_reads = [args.pass2_unmapped1]
    if paired:
        pass2_unmapped_reads.append(args.pass2_unmapped2)
    steps.append(Step(
        cmd=[star, "--runThreadN", str(threads),
             "--genomeDir", args.genome_index,
             "--readFilesIn", *pass2_unmapped_reads,
             "--readFilesCommand", "zcat",
             *_star_pass_options("pass3", paired, args,
                                 args.force_end_to_end, args.hard_clip_5p),
             "--outSAMtype", "BAM", "SortedByCoordinate",
             "--outFileNamePrefix", pass3b_prefix],
        desc="pass3b: local alignment of unmapped reads",
    ))
    steps.append(Step(
        cmd=["mv", pass3b_prefix + "Aligned.sortedByCoord.out.bam", args.pass3b_bam],
        desc="pass3b: move BAM",
    ))
    steps.append(Step(
        cmd=[samtools, "index", "-@", str(threads), args.pass3b_bam],
        desc="pass3b: index BAM",
    ))

    # ── Merge and index final BAM ───────────────────────────────────
    steps.append(Step(
        cmd=[samtools, "merge", "-f", "-@", str(threads), args.bam,
             args.pass3a_bam, args.pass3b_bam],
        desc="merge pass3a + pass3b BAMs",
    ))
    steps.append(Step(
        cmd=[samtools, "index", "-@", str(threads), args.bam],
        desc="index final BAM",
    ))

    return steps


# ── Execution ────────────────────────────────────────────────────────────

def _run_step(step: Step, log_handle) -> None:
    """Execute one Step, logging the command and output."""
    cmd_str = shlex.join(step.cmd)
    if step.pipe_to:
        cmd_str += " | " + shlex.join(step.pipe_to)
    if step.stdout:
        cmd_str += " > " + shlex.quote(step.stdout)
    log_handle.write(f"\n{'=' * 60}\n{step.desc}\n$ {cmd_str}\n")
    log_handle.flush()

    if step.pipe_to:
        # Pipeline: cmd | pipe_to
        p1 = subprocess.Popen(step.cmd, stdout=subprocess.PIPE)
        assert p1.stdout is not None
        p2 = subprocess.Popen(step.pipe_to, stdin=p1.stdout,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()
        stdout, stderr = p2.communicate()
        if p1.wait() != 0:
            raise subprocess.CalledProcessError(p1.returncode, step.cmd)
        if p2.returncode != 0:
            raise subprocess.CalledProcessError(p2.returncode, step.pipe_to, stdout, stderr)
        if stdout:
            log_handle.write(stdout.decode("utf-8", errors="replace"))
        if stderr:
            log_handle.write(stderr.decode("utf-8", errors="replace"))
    elif step.stdout:
        # Redirect stdout to file
        with open(step.stdout, "wb") as out:
            subprocess.run(step.cmd, stdout=out, stderr=log_handle, check=True)
    else:
        result = subprocess.run(step.cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=True)
        if result.stdout and not step.skip_log:
            log_handle.write(result.stdout.decode("utf-8", errors="replace"))
    log_handle.flush()


# ── CLI ──────────────────────────────────────────────────────────────────

def _add_pass_args(group, prefix: str) -> None:
    """Add STAR parameters common to all passes for one pass group."""
    group.add_argument(f"--{prefix}-out-filter-multimap-nmax", type=int,
                       default=None, metavar="N",
                       help="STAR --outFilterMultimapNmax")
    group.add_argument(f"--{prefix}-out-filter-multimap-score-range", type=int,
                       default=None, metavar="N",
                       help="STAR --outFilterMultimapScoreRange")
    group.add_argument(f"--{prefix}-out-filter-mismatch-nover-lmax", type=float,
                       default=None, metavar="F",
                       help="STAR --outFilterMismatchNoverLmax")
    group.add_argument(f"--{prefix}-align-intron-min", type=int,
                       default=None, metavar="N",
                       help="STAR --alignIntronMin")
    group.add_argument(f"--{prefix}-align-mates-gap-max", type=int,
                       default=None, metavar="N",
                       help="STAR --alignMatesGapMax")
    group.add_argument(f"--{prefix}-align-ends-type",
                       choices=("Local", "EndToEnd"), default=None,
                       help="STAR --alignEndsType")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Canonical three-pass STAR alignment for small ncRNA reads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # I/O
    parser.add_argument("--fastq", nargs="+", required=True,
                        help="Input FASTQ file(s)")
    parser.add_argument("--genome-index", required=True,
                        help="Whole-genome STAR index directory")
    parser.add_argument("--smallrna-index", required=True,
                        help="Small-RNA STAR index directory")
    parser.add_argument("--smallrna-bed", required=True,
                        help="BED file of small-RNA gene loci")
    parser.add_argument("--pass1-bam", required=True)
    parser.add_argument("--pass1-fq1", required=True)
    parser.add_argument("--pass1-fq2", required=True)
    parser.add_argument("--pass2-bam", required=True)
    parser.add_argument("--pass2-unmapped1", required=True)
    parser.add_argument("--pass2-unmapped2", required=True)
    parser.add_argument("--pass3a-bam", required=True)
    parser.add_argument("--pass3b-bam", required=True)
    parser.add_argument("--bam", required=True, help="Final merged BAM output")
    parser.add_argument("--log", required=True, help="Log file path")

    # Tool binaries
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--star", default="STAR")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bedtools", default="bedtools")

    # Global behaviour flags
    parser.add_argument("--force-end-to-end", action="store_true", default=False,
                        help="Override alignEndsType to EndToEnd for all passes")
    parser.add_argument("--hard-clip-5p", type=int, default=0, metavar="N",
                        help="If >0, override clip5pNbases to N (PE: 'N 0')")

    # ── Per-pass STAR parameters (literature defaults applied in code) ──
    g1 = parser.add_argument_group("pass 1 STAR parameters (whole-genome E2E)")
    _add_pass_args(g1, "pass1")
    g1.set_defaults(
        pass1_out_filter_multimap_nmax=1000,
        pass1_align_intron_min=9999999,
        pass1_out_filter_multimap_score_range=1,
        pass1_out_filter_mismatch_nover_lmax=0.2,
    )

    g2 = parser.add_argument_group("pass 2 STAR parameters (small-RNA E2E)")
    _add_pass_args(g2, "pass2")
    g2.add_argument("--pass2-out-filter-mismatch-nover-read-lmax", type=float,
                    default=0.05, metavar="F",
                    help="STAR --outFilterMismatchNoverReadLmax")
    g2.add_argument("--pass2-clip5p-nbases", default=None, metavar="STR",
                    help="STAR --clip5pNbases (default: '20 0' PE, '20' SE)")
    g2.add_argument("--pass2-clip3p-nbases", default=None, metavar="STR",
                    help="STAR --clip3pNbases (default: '0 20' PE, '0' SE)")
    g2.set_defaults(
        pass2_out_filter_multimap_nmax=1000,
        pass2_out_filter_multimap_score_range=0,
        pass2_out_filter_mismatch_nover_lmax=0.2,
        pass2_align_intron_min=9999999,
        pass2_align_mates_gap_max=500,
        pass2_align_ends_type="EndToEnd",
    )

    g3 = parser.add_argument_group("pass 3 STAR parameters (genome local)")
    _add_pass_args(g3, "pass3")
    g3.set_defaults(
        pass3_out_filter_multimap_nmax=1000,
        pass3_out_filter_multimap_score_range=0,
        pass3_out_filter_mismatch_nover_lmax=0.025,
        pass3_align_intron_min=9999999,
        pass3_align_mates_gap_max=500,
        pass3_align_ends_type="Local",
    )

    return parser.parse_args()


def main() -> None:
    """Execute the three-pass alignment pipeline."""
    args = parse_args()
    log_path = args.log
    paired = len(args.fastq) == 2

    for path in (
        args.pass1_bam, args.pass1_fq1, args.pass2_bam,
        args.pass2_unmapped1, args.pass3a_bam, args.pass3b_bam,
        args.bam, log_path,
    ):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    steps = build_steps(args, paired)

    # Write command log for reproducibility
    stamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    cmd_log = os.path.join(os.path.dirname(args.bam), f"three_pass_{stamp}.cmd.log")

    try:
        with open(log_path, "w", encoding="utf-8") as log_handle:
            log_handle.write(f"three_pass_align.py — {stamp}\n")
            log_handle.write(f"paired={paired}\n")

            # Write command log header
            with open(cmd_log, "w", encoding="utf-8") as cmd_handle:
                cmd_handle.write(f"# three_pass_align.py command log — {stamp}\n")
                cmd_handle.write(f"# paired={paired}\n")
                for step in steps:
                    line = shlex.join(step.cmd)
                    if step.pipe_to:
                        line += " | " + shlex.join(step.pipe_to)
                    if step.stdout:
                        line += " > " + shlex.quote(step.stdout)
                    cmd_handle.write(line + "\n")

            for step in steps:
                _run_step(step, log_handle)

            log_handle.write(f"\n{'=' * 60}\nAll steps completed successfully.\n")
    except Exception as exc:
        with open(log_path, "a", encoding="utf-8") as handle:
            handle.write(f"three-pass alignment failed: {exc}\n")
        raise


if __name__ == "__main__":
    main()
