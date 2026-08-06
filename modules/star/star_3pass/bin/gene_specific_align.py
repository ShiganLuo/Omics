#!/usr/bin/env python3
"""Strict gene-specific three-pass local alignment and Tailer after canonical star_3pass.

Directly executes the pipeline via subprocess:

  1. Creates per-gene FASTQ, reference, and annotation via prepare_gene_inputs.py
  2. For each gene:
     a. Builds STAR index (pass 1 parameters)
     b. Pass 1: local alignment of all gene reads
     c. Pass 2: local alignment of pass-1 unmapped reads
     d. Pass 3: local alignment of pass-2 unmapped reads
     e. Runs Tailer on the pass-3 BAM
  3. Merges pass-3 gene BAMs and concatenates tail CSVs

A command log is written for reproducibility.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import time
from dataclasses import dataclass


@dataclass
class Step:
    """One executable step.

    ``cmd`` is the primary command.  If ``pipe_to`` is set, stdout of ``cmd``
    is piped into ``pipe_to``.
    """
    cmd: list[str]
    desc: str = ""
    pipe_to: list[str] | None = None


def _run_step(step: Step, log_handle) -> None:
    """Execute one Step, logging the command and output."""
    cmd_str = shlex.join(step.cmd)
    if step.pipe_to:
        cmd_str += " | " + shlex.join(step.pipe_to)
    log_handle.write(f"\n{'=' * 60}\n{step.desc}\n$ {cmd_str}\n")
    log_handle.flush()

    if step.pipe_to:
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
    else:
        result = subprocess.run(step.cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=True)
        if result.stdout:
            log_handle.write(result.stdout.decode("utf-8", errors="replace"))
    log_handle.flush()


def _read_manifest(manifest_path: str) -> list[dict[str, str]]:
    """Read the per-gene manifest TSV produced by prepare_gene_inputs.py.

    Columns: gene_id, ensembl_id, fastq1, fastq2, fasta, gtf, layout, assigned_records
    """
    genes: list[dict[str, str]] = []
    with open(manifest_path, encoding="utf-8") as handle:
        header = None
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            genes.append(dict(zip(header, fields)))
    return genes


# ── Per-pass STAR option assembly ────────────────────────────────────────

_PASS_PARAMS = [
    ("out_filter_multimap_nmax", "--outFilterMultimapNmax", int),
    ("out_filter_multimap_score_range", "--outFilterMultimapScoreRange", int),
    ("out_filter_mismatch_nover_lmax", "--outFilterMismatchNoverLmax", float),
    ("align_intron_min", "--alignIntronMin", int),
    ("align_mates_gap_max", "--alignMatesGapMax", int),
    ("align_ends_type", "--alignEndsType", str),
]


def _star_pass_options(pass_name: str, args: argparse.Namespace) -> list[str]:
    """Assemble STAR flag tokens for one pass from typed CLI arguments."""
    tokens: list[str] = []
    for attr, flag, cast in _PASS_PARAMS:
        val = getattr(args, f"{pass_name}_{attr}")
        if val is not None:
            tokens.append(flag)
            tokens.append(str(cast(val)))
    return tokens


def _build_align_cmd(args: argparse.Namespace, pass_name: str,
                     index_dir: str, read_files: list[str],
                     out_prefix: str) -> list[str]:
    """Build a STAR alignment command for one pass."""
    cmd = [
        args.star, "--runThreadN", str(args.threads),
        "--genomeDir", index_dir,
        "--readFilesIn", *read_files,
        "--readFilesCommand", "zcat",
        *_star_pass_options(pass_name, args),
        "--outReadsUnmapped", "Fastx",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFileNamePrefix", out_prefix,
    ]
    return cmd


def main() -> None:
    """Execute the gene-specific three-pass alignment pipeline."""
    args = parse_args()
    log_path = args.log

    os.makedirs(os.path.dirname(args.output_bam), exist_ok=True)
    os.makedirs(args.input_dir, exist_ok=True)

    stamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    cmd_log = os.path.join(
        os.path.dirname(args.output_bam), f"gene_specific_{stamp}.cmd.log"
    )

    try:
        with open(log_path, "w", encoding="utf-8") as log_handle:
            log_handle.write(f"gene_specific_align.py - {stamp}\n")

            all_cmds: list[str] = []

            # ── Step 1: Prepare per-gene inputs ──────────────────────
            prep_cmd = [
                "python", args.helper,
                "--bam", args.input_bam,
                "--bed", args.input_bed,
                "--genome-fasta", args.input_genome_fasta,
                "--outdir", args.input_dir,
                "--samtools", args.samtools,
                "--bedtools", args.bedtools,
                "--ambiguous", args.ambiguous,
                "--flank", str(args.flank),
            ]
            all_cmds.append(shlex.join(prep_cmd))
            _run_step(Step(cmd=prep_cmd, desc="prepare per-gene inputs"), log_handle)

            # ── Step 2: Per-gene three-pass loop ─────────────────────
            manifest_path = os.path.join(args.input_dir, "genes.tsv")
            genes = _read_manifest(manifest_path)

            gene_bams: list[str] = []
            gene_tails: list[str] = []

            for gene in genes:
                gene_id = gene["gene_id"]
                fq1 = gene["fastq1"]
                fq2 = gene["fastq2"]
                fasta = gene["fasta"]
                gtf = gene["gtf"]
                gene_dir = os.path.dirname(fasta)
                index_dir = os.path.join(gene_dir, "index")
                os.makedirs(index_dir, exist_ok=True)

                # ── STAR genomeGenerate (uses pass1 params) ─────────
                gen_cmd = [
                    args.star, "--runMode", "genomeGenerate",
                    "--runThreadN", str(args.threads),
                    "--genomeDir", index_dir,
                    "--genomeFastaFiles", fasta,
                    "--genomeSAindexNbases", str(args.pass1_genome_sa_index_nbases),
                ]
                all_cmds.append(shlex.join(gen_cmd))
                _run_step(Step(cmd=gen_cmd, desc=f"[{gene_id}] genomeGenerate"), log_handle)

                # ── Determine initial read files ─────────────────────
                current_reads = [fq1]
                if os.path.getsize(fq2) > 0:
                    current_reads.append(fq2)

                # ── Three passes: each uses prev pass's unmapped reads ─
                final_bam = None
                for pass_name, pass_num in [("pass1", 1), ("pass2", 2), ("pass3", 3)]:
                    pass_prefix = os.path.join(gene_dir, f"pass{pass_num}.")
                    align_cmd = _build_align_cmd(
                        args, pass_name, index_dir, current_reads, pass_prefix
                    )
                    all_cmds.append(shlex.join(align_cmd))
                    _run_step(Step(
                        cmd=align_cmd,
                        desc=f"[{gene_id}] pass{pass_num} local alignment",
                    ), log_handle)

                    pass_bam = pass_prefix + "Aligned.sortedByCoord.out.bam"
                    if not os.path.exists(pass_bam) or os.path.getsize(pass_bam) == 0:
                        log_handle.write(f"[{gene_id}] pass{pass_num} no aligned BAM, skipping gene\n")
                        final_bam = None
                        break

                    final_bam = pass_bam

                    # Feed unmapped reads into next pass
                    if pass_num < 3:
                        unmapped_mate1 = pass_prefix + "Unmapped.out.mate1"
                        unmapped_mate2 = pass_prefix + "Unmapped.out.mate2"
                        if os.path.exists(unmapped_mate2) and os.path.getsize(unmapped_mate2) > 0:
                            current_reads = [unmapped_mate1, unmapped_mate2]
                        else:
                            current_reads = [unmapped_mate1]

                if final_bam is None:
                    continue

                # ── Use the last non-empty pass BAM for Tailer ──────
                # STAR always emits a BAM (even with 0 reads, the header is
                # non-zero), so we must check read count, not file size.
                # Walk backwards from pass3 to find the last pass with reads.
                tail_bam = final_bam
                for pass_num in (3, 2, 1):
                    candidate = os.path.join(gene_dir, f"pass{pass_num}.Aligned.sortedByCoord.out.bam")
                    count_cmd = [args.samtools, "view", "-c", candidate]
                    result = subprocess.run(count_cmd, stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE, check=True)
                    read_count = int(result.stdout.decode().strip())
                    if read_count > 0:
                        tail_bam = candidate
                        log_handle.write(f"[{gene_id}] using pass{pass_num} BAM ({read_count} reads) for Tailer\n")
                        break
                else:
                    log_handle.write(f"[{gene_id}] all passes have 0 reads, skipping Tailer\n")
                    continue

                # ── samtools index the BAM for Tailer ───────────────
                idx_cmd = [args.samtools, "index", tail_bam]
                all_cmds.append(shlex.join(idx_cmd))
                _run_step(Step(cmd=idx_cmd, desc=f"[{gene_id}] index BAM for Tailer"), log_handle)

                # ── Tailer on the final non-empty BAM ───────────────
                tail_cmd = [
                    args.tailer, "-a", gtf,
                    "-read", "1", "-t", str(args.tailer_threshold),
                ]
                if args.tailer_rev_comp:
                    tail_cmd.append("--rev_comp")
                tail_cmd.append(tail_bam)
                all_cmds.append(shlex.join(tail_cmd))
                _run_step(Step(cmd=tail_cmd, desc=f"[{gene_id}] Tailer"), log_handle)

                gene_tail = tail_bam.rsplit(".bam", 1)[0] + "_tail.csv"
                if not os.path.exists(gene_tail) or os.path.getsize(gene_tail) == 0:
                    raise RuntimeError(f"[{gene_id}] Tailer produced no output: {gene_tail}")
                gene_bams.append(tail_bam)
                gene_tails.append(gene_tail)

            # ── Step 3: Merge gene BAMs (last non-empty pass) ────────
            if not gene_bams:
                raise RuntimeError("No gene BAMs were produced")

            merge_cmd = [
                args.samtools, "merge", "-f", "-@", str(args.threads),
                args.output_bam, *gene_bams,
            ]
            all_cmds.append(shlex.join(merge_cmd))
            _run_step(Step(cmd=merge_cmd, desc="merge pass3 gene BAMs"), log_handle)

            index_cmd = [
                args.samtools, "index", "-@", str(args.threads), args.output_bam,
            ]
            all_cmds.append(shlex.join(index_cmd))
            _run_step(Step(cmd=index_cmd, desc="index merged BAM"), log_handle)

            # ── Step 4: Concatenate tail CSVs ────────────────────────
            with open(args.output_tail, "w", encoding="utf-8") as out:
                for i, csv_path in enumerate(gene_tails):
                    with open(csv_path, encoding="utf-8") as csv_handle:
                        if i == 0:
                            out.write(csv_handle.read())
                        else:
                            next(csv_handle, None)
                            out.write(csv_handle.read())
            all_cmds.append(f"# concatenated {len(gene_tails)} tail CSVs > {args.output_tail}")
            log_handle.write(f"\n{'=' * 60}\nconcatenated {len(gene_tails)} tail CSVs > {args.output_tail}\n")

            log_handle.write(f"\n{'=' * 60}\nAll steps completed successfully.\n")

            with open(cmd_log, "w", encoding="utf-8") as cmd_handle:
                cmd_handle.write(f"# gene_specific_align.py command log - {stamp}\n")
                cmd_handle.write(f"# genes processed: {len(gene_bams)}\n")
                for line in all_cmds:
                    cmd_handle.write(line + "\n")

    except Exception as exc:
        with open(log_path, "a", encoding="utf-8") as handle:
            handle.write(f"gene-specific alignment/Tailer failed: {exc}\n")
        raise


# ── CLI ──────────────────────────────────────────────────────────────────

def _add_pass_args(group, prefix: str) -> None:
    """Add STAR alignment parameters for one pass group."""
    group.add_argument(f"--{prefix}-out-filter-multimap-nmax", type=int,
                       default=None, metavar="N", help="STAR --outFilterMultimapNmax")
    group.add_argument(f"--{prefix}-out-filter-multimap-score-range", type=int,
                       default=None, metavar="N", help="STAR --outFilterMultimapScoreRange")
    group.add_argument(f"--{prefix}-out-filter-mismatch-nover-lmax", type=float,
                       default=None, metavar="F", help="STAR --outFilterMismatchNoverLmax")
    group.add_argument(f"--{prefix}-align-intron-min", type=int,
                       default=None, metavar="N", help="STAR --alignIntronMin")
    group.add_argument(f"--{prefix}-align-mates-gap-max", type=int,
                       default=None, metavar="N", help="STAR --alignMatesGapMax")
    group.add_argument(f"--{prefix}-align-ends-type",
                       choices=("Local", "EndToEnd"), default=None,
                       help="STAR --alignEndsType")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Gene-specific three-pass local alignment and Tailer after star_3pass.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # I/O
    parser.add_argument("--input-bam", required=True,
                        help="Final BAM from canonical star_3pass")
    parser.add_argument("--input-bed", required=True,
                        help="BED file of small-RNA gene loci")
    parser.add_argument("--input-genome-fasta", required=True,
                        help="Genome FASTA for reference extraction")
    parser.add_argument("--output-bam", required=True,
                        help="Merged gene-specific BAM output (pass3 only)")
    parser.add_argument("--output-tail", required=True,
                        help="Concatenated Tailer CSV output")
    parser.add_argument("--input-dir", required=True,
                        help="Directory for per-gene intermediate inputs")
    parser.add_argument("--helper", required=True,
                        help="Path to prepare_gene_inputs.py")
    parser.add_argument("--log", required=True)

    # Tool binaries
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--star", default="STAR")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bedtools", default="bedtools")
    parser.add_argument("--tailer", default="Tailer")

    # Gene assignment parameters
    parser.add_argument("--ambiguous", choices=("exclude", "duplicate"),
                        default="exclude",
                        help="How to handle reads mapping to multiple genes")
    parser.add_argument("--flank", type=int, default=50, metavar="N",
                        help="Flanking bases added to each gene reference")

    # ── Pass 1 STAR parameters (also includes genomeGenerate) ──────
    g1 = parser.add_argument_group("pass 1 STAR parameters")
    g1.add_argument("--pass1-genome-sa-index-nbases", type=int, default=3, metavar="N",
                    help="STAR --genomeSAindexNbases for per-gene index")
    _add_pass_args(g1, "pass1")
    g1.set_defaults(
        pass1_out_filter_multimap_nmax=1000,
        pass1_out_filter_multimap_score_range=0,
        pass1_out_filter_mismatch_nover_lmax=0.2,
        pass1_align_intron_min=9999999,
        pass1_align_mates_gap_max=500,
        pass1_align_ends_type="Local",
    )

    # ── Pass 2 STAR parameters ─────────────────────────────────────
    g2 = parser.add_argument_group("pass 2 STAR parameters")
    _add_pass_args(g2, "pass2")
    g2.set_defaults(
        pass2_out_filter_multimap_nmax=1000,
        pass2_out_filter_multimap_score_range=0,
        pass2_out_filter_mismatch_nover_lmax=0.2,
        pass2_align_intron_min=9999999,
        pass2_align_mates_gap_max=500,
        pass2_align_ends_type="Local",
    )

    # ── Pass 3 STAR parameters ─────────────────────────────────────
    g3 = parser.add_argument_group("pass 3 STAR parameters")
    _add_pass_args(g3, "pass3")
    g3.set_defaults(
        pass3_out_filter_multimap_nmax=1000,
        pass3_out_filter_multimap_score_range=0,
        pass3_out_filter_mismatch_nover_lmax=0.025,
        pass3_align_intron_min=9999999,
        pass3_align_mates_gap_max=500,
        pass3_align_ends_type="Local",
    )

    # Tailer parameters
    g_tail = parser.add_argument_group("Tailer parameters")
    g_tail.add_argument("--tailer-threshold", type=int, default=50, metavar="N",
                        help="Tailer -t threshold")
    g_tail.add_argument("--tailer-rev-comp", action="store_true", default=False,
                        help="Pass --rev_comp to Tailer")

    return parser.parse_args()


if __name__ == "__main__":
    main()
