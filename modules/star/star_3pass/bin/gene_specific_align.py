#!/usr/bin/env python3
"""Strict gene-specific three-pass alignment and Tailer after canonical star_3pass.

Directly executes the pipeline via subprocess:

  1. Creates per-gene FASTQ, reference, and annotation via prepare_gene_inputs.py
  2. For each gene:
     a. Builds STAR index (per-gene reference)
     b. Pass 1: per-gene local alignment of all gene reads
     c. Pass 2: per-gene end-to-end alignment with clip5p/clip3p of all gene reads
        -> split into mapped and unmapped reads
     d. Pass 3a: per-gene local alignment of pass-2 mapped reads
     e. Pass 3b: per-gene local alignment of pass-2 unmapped reads
     f. Merges pass-3a and pass-3b BAMs
     g. Runs Tailer on the merged gene BAM
  3. Merges all gene BAMs and concatenates tail CSVs

A command log is written for reproducibility.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass
from typing import List, Optional, Tuple


# =========================================================================
# 1. STAR 参数中心化配置字典 (数据驱动)
# 想要增删改 STAR 参数？只需在这里修改，命令行和运行逻辑会自动同步！
# 值为 None 代表该参数默认不传给 STAR，但允许在命令行中指定。
# =========================================================================

PASS1_STAR_DEFAULTS = {
    "--outFilterMultimapNmax": 1000,
    "--outFilterMultimapScoreRange": 1,
    "--outFilterMismatchNoverLmax": 0.2,
    "--alignIntronMin": 9999999,
    "--alignEndsType": "Local",
}

PASS2_STAR_DEFAULTS = {
    "--outFilterMultimapNmax": 1000,
    "--outFilterMultimapScoreRange": 0,
    "--outFilterMismatchNoverLmax": 0.2,
    "--outFilterMismatchNoverReadLmax": 0.05,
    "--alignIntronMin": 9999999,
    "--alignMatesGapMax": 500,
    "--alignEndsType": "EndToEnd",
    "--clip5pNbases": None, 
    "--clip3pNbases": None,
    "--outReadsUnmapped": "Fastx"
}

PASS3_STAR_DEFAULTS = {
    "--outFilterMultimapNmax": 1000,
    "--outFilterMultimapScoreRange": 0,
    "--outFilterMismatchNoverLmax": 0.025,
    "--alignIntronMin": 9999999,
    "--alignMatesGapMax": 500,
    "--alignEndsType": "Local",
}

STAR_PASS_CONFIGS = {
    "pass1": PASS1_STAR_DEFAULTS,
    "pass2": PASS2_STAR_DEFAULTS,
    "pass3": PASS3_STAR_DEFAULTS,
}


# =========================================================================
# 2. Pipeline 核心执行逻辑
# =========================================================================

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
                                stderr=subprocess.STDOUT, check=False)
        if result.stdout:
            log_handle.write(result.stdout.decode("utf-8", errors="replace"))
        log_handle.flush()
        if result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode, step.cmd, output=result.stdout)
    log_handle.flush()


def _read_manifest(manifest_path: str) -> list[dict[str, str]]:
    """Read the per-gene manifest TSV produced by prepare_gene_inputs.py."""
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


def _star_pass_options(pass_name: str, args: argparse.Namespace) -> list[str]:
    """从 argparse.Namespace 中提取并组装对应 Pass 的 STAR 运行参数"""
    tokens: list[str] = []
    defaults = STAR_PASS_CONFIGS[pass_name]
    
    # 1. 提取字典中已定义的参数
    for star_flag in defaults.keys():
        stripped_flag = star_flag.lstrip("-")
        # argparse 会将 "--pass1-alignIntronMin" 解析为 "pass1_alignIntronMin"
        attr_name = f"{pass_name}_{stripped_flag}"
        val = getattr(args, attr_name)
        
        if val is not None:
            tokens.append(star_flag)
            tokens.append(str(val))
            
    # 2. 解析并通过 shlex 注入额外的 raw 字符串参数 (万能通道)
    extra_args_str = getattr(args, f"{pass_name}_extra_args")
    if extra_args_str:
        tokens.extend(shlex.split(extra_args_str))
        
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
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFileNamePrefix", out_prefix,
    ]
    return cmd


def _bam_has_reads(samtools: str, bam_path: str) -> Tuple[bool, int]:
    """Check if a BAM file has any reads. Returns (has_reads, count)."""
    if not os.path.exists(bam_path) or os.path.getsize(bam_path) == 0:
        return False, 0
    result = subprocess.run([samtools, "view", "-c", bam_path],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    count = int(result.stdout.decode().strip())
    return count > 0, count


def _extract_fastq_from_bam(samtools: str, threads: int, bam_path: str,
                            fq1: str, fq2: str, paired: bool) -> list[list[str]]:
    """Return samtools commands to name-sort a BAM and extract FASTQ."""
    name_sorted = bam_path.rsplit(".bam", 1)[0] + ".name_sorted.bam"
    cmds = [
        [samtools, "sort", "-n", "-@", str(threads), "-o", name_sorted, bam_path],
    ]
    if paired:
        cmds.append([samtools, "fastq", "-@", str(threads), "-n",
                      "-1", fq1, "-2", fq2,
                      "-0", "/dev/null", "-s", "/dev/null", name_sorted])
    else:
        cmds.append([samtools, "fastq", "-@", str(threads), "-n",
                      "-0", fq1, "-s", "/dev/null", name_sorted])
    return cmds


# =========================================================================
# 3. Main 函数 
# =========================================================================

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

                paired = os.path.exists(fq2) and os.path.getsize(fq2) > 0
                initial_reads = [fq1]
                if paired:
                    initial_reads.append(fq2)

                # ── STAR genomeGenerate ───────────────────────────────
                gen_prefix = os.path.join(gene_dir, "gen.")
                gen_cmd = [
                    args.star, "--runMode", "genomeGenerate",
                    "--runThreadN", str(args.threads),
                    "--genomeDir", index_dir,
                    "--genomeFastaFiles", fasta,
                    "--genomeSAindexNbases", str(args.index_genome_sa_index_nbases),
                    "--outFileNamePrefix", gen_prefix,
                ]
                all_cmds.append(shlex.join(gen_cmd))
                _run_step(Step(cmd=gen_cmd, desc=f"[{gene_id}] genomeGenerate"), log_handle)

                # ── Pass 1: per-gene E2E alignment of all reads ──────
                pass1_prefix = os.path.join(gene_dir, "pass1.")
                pass1_cmd = _build_align_cmd(
                    args, "pass1", index_dir, initial_reads, pass1_prefix
                )
                all_cmds.append(shlex.join(pass1_cmd))
                _run_step(Step(cmd=pass1_cmd, desc=f"[{gene_id}] pass1 E2E alignment"), log_handle)

                pass1_bam = pass1_prefix + "Aligned.sortedByCoord.out.bam"
                has_reads, n_reads = _bam_has_reads(args.samtools, pass1_bam)
                if not has_reads:
                    log_handle.write(f"[{gene_id}] pass1 no aligned reads, skipping gene\n")
                    continue
                log_handle.write(f"[{gene_id}] pass1: {n_reads} reads\n")

                # ── Pass 2: per-gene E2E+clip alignment of all reads ─
                pass2_prefix = os.path.join(gene_dir, "pass2.")
                pass2_cmd = _build_align_cmd(
                    args, "pass2", index_dir, initial_reads, pass2_prefix
                )
                all_cmds.append(shlex.join(pass2_cmd))
                _run_step(Step(cmd=pass2_cmd, desc=f"[{gene_id}] pass2 E2E+clip alignment"), log_handle)

                pass2_bam = pass2_prefix + "Aligned.sortedByCoord.out.bam"
                pass2_unmapped1 = pass2_prefix + "Unmapped.out.mate1"
                pass2_unmapped2 = pass2_prefix + "Unmapped.out.mate2"

                pass2_mapped_fq1 = os.path.join(gene_dir, "pass2_mapped_1.fq.gz")
                pass2_mapped_fq2 = os.path.join(gene_dir, "pass2_mapped_2.fq.gz")
                for fq_cmd in _extract_fastq_from_bam(
                    args.samtools, args.threads, pass2_bam,
                    pass2_mapped_fq1, pass2_mapped_fq2, paired
                ):
                    all_cmds.append(shlex.join(fq_cmd))
                    _run_step(Step(cmd=fq_cmd, desc=f"[{gene_id}] pass2 extract mapped FASTQ"), log_handle)

                pass2_unmapped_gz1 = pass2_unmapped1
                pass2_unmapped_gz2 = pass2_unmapped2
                if os.path.exists(pass2_unmapped1) and os.path.getsize(pass2_unmapped1) > 0:
                    pass2_unmapped_gz1 = pass2_unmapped1 + ".gz"
                    gzip_cmd = ["gzip", "-f", pass2_unmapped1]
                    all_cmds.append(shlex.join(gzip_cmd))
                    _run_step(Step(cmd=gzip_cmd, desc=f"[{gene_id}] pass2 gzip unmapped mate1"), log_handle)
                if paired and os.path.exists(pass2_unmapped2) and os.path.getsize(pass2_unmapped2) > 0:
                    pass2_unmapped_gz2 = pass2_unmapped2 + ".gz"
                    gzip_cmd = ["gzip", "-f", pass2_unmapped2]
                    all_cmds.append(shlex.join(gzip_cmd))
                    _run_step(Step(cmd=gzip_cmd, desc=f"[{gene_id}] pass2 gzip unmapped mate2"), log_handle)

                # ── Pass 3a: per-gene Local alignment of mapped reads ─
                pass3a_prefix = os.path.join(gene_dir, "pass3a.")
                pass3a_reads = [pass2_mapped_fq1]
                if paired and os.path.exists(pass2_mapped_fq2) and os.path.getsize(pass2_mapped_fq2) > 0:
                    pass3a_reads.append(pass2_mapped_fq2)
                pass3a_cmd = _build_align_cmd(
                    args, "pass3", index_dir, pass3a_reads, pass3a_prefix
                )
                all_cmds.append(shlex.join(pass3a_cmd))
                _run_step(Step(cmd=pass3a_cmd, desc=f"[{gene_id}] pass3a Local alignment (mapped)"), log_handle)

                pass3a_bam = pass3a_prefix + "Aligned.sortedByCoord.out.bam"

                # ── Pass 3b: per-gene Local alignment of unmapped reads ─
                pass3b_prefix = os.path.join(gene_dir, "pass3b.")
                pass3b_reads = [pass2_unmapped_gz1]
                if paired and os.path.exists(pass2_unmapped_gz2) and os.path.getsize(pass2_unmapped_gz2) > 0:
                    pass3b_reads.append(pass2_unmapped_gz2)
                pass3b_cmd = _build_align_cmd(
                    args, "pass3", index_dir, pass3b_reads, pass3b_prefix
                )
                all_cmds.append(shlex.join(pass3b_cmd))
                _run_step(Step(cmd=pass3b_cmd, desc=f"[{gene_id}] pass3b Local alignment (unmapped)"), log_handle)

                pass3b_bam = pass3b_prefix + "Aligned.sortedByCoord.out.bam"

                # ── Merge pass3a + pass3b for this gene ───────────────
                gene_merged_bam = os.path.join(gene_dir, "merged.bam")
                merge_inputs = []
                if os.path.exists(pass3a_bam) and os.path.getsize(pass3a_bam) > 0:
                    merge_inputs.append(pass3a_bam)
                if os.path.exists(pass3b_bam) and os.path.getsize(pass3b_bam) > 0:
                    merge_inputs.append(pass3b_bam)
                if not merge_inputs:
                    log_handle.write(f"[{gene_id}] pass3a and pass3b both empty, skipping gene\n")
                    continue
                merge_cmd = [
                    args.samtools, "merge", "-f", "-@", str(args.threads),
                    gene_merged_bam, *merge_inputs,
                ]
                all_cmds.append(shlex.join(merge_cmd))
                _run_step(Step(cmd=merge_cmd, desc=f"[{gene_id}] merge pass3a+3b"), log_handle)

                has_reads, n_reads = _bam_has_reads(args.samtools, gene_merged_bam)
                if not has_reads:
                    log_handle.write(f"[{gene_id}] merged BAM has 0 reads, skipping Tailer\n")
                    continue
                log_handle.write(f"[{gene_id}] merged BAM: {n_reads} reads\n")

                # ── samtools index & Tailer ───────────────────────────
                idx_cmd = [args.samtools, "index", gene_merged_bam]
                all_cmds.append(shlex.join(idx_cmd))
                _run_step(Step(cmd=idx_cmd, desc=f"[{gene_id}] index BAM for Tailer"), log_handle)

                tail_cmd = [
                    args.tailer, "-a", gtf,
                    "-read", "1", "-t", str(args.tailer_threshold),
                ]
                if args.tailer_rev_comp:
                    tail_cmd.append("--rev_comp")
                tail_cmd.append(gene_merged_bam)
                all_cmds.append(shlex.join(tail_cmd))
                _run_step(Step(cmd=tail_cmd, desc=f"[{gene_id}] Tailer"), log_handle)

                gene_tail = gene_merged_bam.rsplit(".bam", 1)[0] + "_tail.csv"
                if not os.path.exists(gene_tail) or os.path.getsize(gene_tail) == 0:
                    raise RuntimeError(f"[{gene_id}] Tailer produced no output: {gene_tail}")
                gene_bams.append(gene_merged_bam)
                gene_tails.append(gene_tail)

            # ── Step 3 & 4: Merge & Concat ─────────────────────────────
            if not gene_bams:
                raise RuntimeError("No gene BAMs were produced")

            merge_cmd = [
                args.samtools, "merge", "-f", "-@", str(args.threads),
                args.output_bam, *gene_bams,
            ]
            all_cmds.append(shlex.join(merge_cmd))
            _run_step(Step(cmd=merge_cmd, desc="merge all gene BAMs"), log_handle)

            index_cmd = [
                args.samtools, "index", "-@", str(args.threads), args.output_bam,
            ]
            all_cmds.append(shlex.join(index_cmd))
            _run_step(Step(cmd=index_cmd, desc="index merged BAM"), log_handle)

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


# =========================================================================
# 4. CLI 参数解析与强校验
# =========================================================================

def valid_existing_file(path_str: str) -> str:
    """Fail-fast: 校验输入文件是否存在"""
    if not os.path.isfile(path_str):
        raise argparse.ArgumentTypeError(f"输入文件不存在: {path_str}")
    return path_str

def valid_executable(exe_str: str) -> str:
    """Fail-fast: 校验工具是否在 PATH 中或为有效路径"""
    if shutil.which(exe_str) is None:
        raise argparse.ArgumentTypeError(f"未在 PATH 中找到可执行文件，或路径无效: {exe_str}")
    return exe_str


def _add_dynamic_star_args(group: argparse._ArgumentGroup, prefix: str, defaults: dict) -> None:
    """遍历字典，动态将 STAR 参数注册到 argparse 中"""
    for star_flag, default_val in defaults.items():
        stripped_flag = star_flag.lstrip("-")
        cli_flag = f"--{prefix}-{stripped_flag}"
        
        # 根据字典中的默认值类型推断参数类型
        if isinstance(default_val, int):
            arg_type = int
        elif isinstance(default_val, float):
            arg_type = float
        else:
            arg_type = str
            
        group.add_argument(cli_flag, type=arg_type, default=default_val, 
                           metavar="VAL", help=f"STAR {star_flag} for {prefix}")

    # 万能逃逸通道，用于随时挂载字典里没有的冷门参数
    group.add_argument(f"--{prefix}-extra-args", type=str, default="",
                       metavar="STR", help=f"Raw extra STAR args for {prefix} (e.g. '--limitBAMsortRAM 10000')")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Gene-specific three-pass alignment and Tailer after star_3pass.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    g_io = parser.add_argument_group("Input/Output parameters")
    g_io.add_argument("--input-bam", type=valid_existing_file, required=True,
                      help="Final BAM from canonical star_3pass")
    g_io.add_argument("--input-bed", type=valid_existing_file, required=True,
                      help="BED file of small-RNA gene loci")
    g_io.add_argument("--input-genome-fasta", type=valid_existing_file, required=True,
                      help="Genome FASTA for reference extraction")
    g_io.add_argument("--helper", type=valid_existing_file, required=True,
                      help="Path to prepare_gene_inputs.py")
    
    g_io.add_argument("--output-bam", required=True, help="Merged gene-specific BAM output")
    g_io.add_argument("--output-tail", required=True, help="Concatenated Tailer CSV output")
    g_io.add_argument("--input-dir", required=True, help="Directory for per-gene intermediate inputs")
    g_io.add_argument("--log", required=True)

    # ── Tools ─────────────────────────────────────
    g_tools = parser.add_argument_group("Tool binaries & Compute")
    g_tools.add_argument("--threads", type=int, default=2)
    g_tools.add_argument("--star", type=valid_executable, default="STAR")
    g_tools.add_argument("--samtools", type=valid_executable, default="samtools")
    g_tools.add_argument("--bedtools", type=valid_executable, default="bedtools")
    g_tools.add_argument("--tailer", type=valid_executable, default="Tailer")

    # ── Gene Assignment ───────────────────────────
    g_gene = parser.add_argument_group("Gene assignment parameters")
    g_gene.add_argument("--ambiguous", choices=("exclude", "duplicate"), default="exclude",
                        help="How to handle reads mapping to multiple genes")
    g_gene.add_argument("--flank", type=int, default=50, metavar="N",
                        help="Flanking bases added to each gene reference")
    g_gene.add_argument("--index-genome-sa-index-nbases", type=int, default=3, metavar="N",
                        help="STAR --genomeSAindexNbases for per-gene index")

    # ── Dynamic STAR Params ───────────────────────
    for pass_name in ["pass1", "pass2", "pass3"]:
        group = parser.add_argument_group(f"{pass_name} STAR parameters")
        _add_dynamic_star_args(group, pass_name, STAR_PASS_CONFIGS[pass_name])

    # ── Tailer Params ─────────────────────────────
    g_tail = parser.add_argument_group("Tailer parameters")
    g_tail.add_argument("--tailer-threshold", type=int, default=50, metavar="N",
                        help="Tailer -t threshold")
    
    if hasattr(argparse, 'BooleanOptionalAction'):
        g_tail.add_argument("--tailer-rev-comp", action=argparse.BooleanOptionalAction, default=False,
                            help="Pass --rev_comp to Tailer")
    else:
        g_tail.add_argument("--tailer-rev-comp", action="store_true", default=False,
                            help="Pass --rev_comp to Tailer")

    return parser.parse_args()


if __name__ == "__main__":
    main()