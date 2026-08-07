#!/usr/bin/env python3
"""Canonical three-pass STAR alignment for small non-coding RNA reads.

Directly executes the alignment pipeline via subprocess:

  1. Whole-genome local alignment (pass 1)
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
import shutil
import subprocess
import time
from dataclasses import dataclass


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
    "--clip5pNbases": '20 0',  # 动态默认值：PE为 '20 0', SE为 '20'
    "--clip3pNbases": '0 20',  # 动态默认值：PE为 '0 20', SE为 '0'
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


def _star_pass_options(pass_name: str, paired: bool, args: argparse.Namespace) -> list[str]:
    """从 argparse.Namespace 中动态提取并组装对应 Pass 的 STAR 运行参数。"""
    tokens: list[str] = []
    defaults = STAR_PASS_CONFIGS[pass_name]
    
    for star_flag in defaults.keys():
        stripped_flag = star_flag.lstrip("-")
        attr_name = f"{pass_name}_{stripped_flag}"
        val = getattr(args, attr_name)                
        if val is not None:
            tokens.append(star_flag)
            # 兼容 "20 0" 这种带空格的值，打散为多个 token 传给 subprocess
            tokens.extend(str(val).split())
            
    # 解析万能逃逸通道 extra-args
    extra_args_str = getattr(args, f"{pass_name}_extra_args")
    if extra_args_str:
        tokens.extend(shlex.split(extra_args_str))
        
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


# =========================================================================
# 2. Pipeline 执行核心
# =========================================================================

@dataclass
class Step:
    """One executable step in the pipeline."""
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
             *_star_pass_options("pass1", paired, args),
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
             *_star_pass_options("pass2", paired, args),
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
             *_star_pass_options("pass3", paired, args),
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
             *_star_pass_options("pass3", paired, args),
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
        with open(step.stdout, "wb") as out:
            subprocess.run(step.cmd, stdout=out, stderr=log_handle, check=True)
    else:
        result = subprocess.run(step.cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=True)
        if result.stdout and not step.skip_log:
            log_handle.write(result.stdout.decode("utf-8", errors="replace"))
    log_handle.flush()


# =========================================================================
# 3. CLI 参数解析与强校验
# =========================================================================

def valid_existing_file(path_str: str) -> str:
    """Fail-fast: 校验输入文件是否存在"""
    if not os.path.isfile(path_str):
        raise argparse.ArgumentTypeError(f"文件不存在: {path_str}")
    return path_str

def valid_existing_dir(path_str: str) -> str:
    """Fail-fast: 校验输入目录是否存在"""
    if not os.path.isdir(path_str):
        raise argparse.ArgumentTypeError(f"目录不存在: {path_str}")
    return path_str

def valid_executable(exe_str: str) -> str:
    """Fail-fast: 校验工具是否在 PATH 中"""
    if shutil.which(exe_str) is None:
        raise argparse.ArgumentTypeError(f"未找到可执行文件: {exe_str}")
    return exe_str


def _add_dynamic_star_args(group: argparse._ArgumentGroup, prefix: str, defaults: dict) -> None:
    """遍历字典，动态将 STAR 参数注册到 argparse 中"""
    for star_flag, default_val in defaults.items():
        stripped_flag = star_flag.lstrip("-")
        cli_flag = f"--{prefix}-{stripped_flag}"
        
        if isinstance(default_val, int):
            arg_type = int
        elif isinstance(default_val, float):
            arg_type = float
        else:
            arg_type = str
            
        group.add_argument(cli_flag, type=arg_type, default=default_val, 
                           metavar="VAL", help=f"STAR {star_flag} for {prefix}")

    # 万能逃逸通道
    group.add_argument(f"--{prefix}-extra-args", type=str, default="",
                       metavar="STR", help=f"Raw extra STAR args for {prefix} (e.g. '--limitBAMsortRAM 10000')")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Canonical three-pass STAR alignment for small ncRNA reads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── I/O 参数 (包含防呆校验) ──
    g_io = parser.add_argument_group("Input/Output parameters")
    g_io.add_argument("--fastq", nargs="+", type=valid_existing_file, required=True,
                      help="Input FASTQ file(s)")
    g_io.add_argument("--genome-index", type=valid_existing_dir, required=True,
                      help="Whole-genome STAR index directory")
    g_io.add_argument("--smallrna-index", type=valid_existing_dir, required=True,
                      help="Small-RNA STAR index directory")
    g_io.add_argument("--smallrna-bed", type=valid_existing_file, required=True,
                      help="BED file of small-RNA gene loci")
    
    g_io.add_argument("--pass1-bam", required=True)
    g_io.add_argument("--pass1-fq1", required=True)
    g_io.add_argument("--pass1-fq2", required=True)
    g_io.add_argument("--pass2-bam", required=True)
    g_io.add_argument("--pass2-unmapped1", required=True)
    g_io.add_argument("--pass2-unmapped2", required=True)
    g_io.add_argument("--pass3a-bam", required=True)
    g_io.add_argument("--pass3b-bam", required=True)
    g_io.add_argument("--bam", required=True, help="Final merged BAM output")
    g_io.add_argument("--log", required=True, help="Log file path")

    # ── Tools ──
    g_tools = parser.add_argument_group("Tool binaries & Compute")
    g_tools.add_argument("--threads", type=int, default=12)
    g_tools.add_argument("--star", type=valid_executable, default="STAR")
    g_tools.add_argument("--samtools", type=valid_executable, default="samtools")
    g_tools.add_argument("--bedtools", type=valid_executable, default="bedtools")

    # ── 动态 STAR 参数 ──
    for pass_name in ["pass1", "pass2", "pass3"]:
        group = parser.add_argument_group(f"{pass_name} STAR parameters")
        _add_dynamic_star_args(group, pass_name, STAR_PASS_CONFIGS[pass_name])

    return parser.parse_args()


# =========================================================================
# 4. 主程序入口
# =========================================================================

def main() -> None:
    """Execute the three-pass alignment pipeline."""
    args = parse_args()
    log_path = args.log
    paired = len(args.fastq) == 2

    # 自动创建所有输出文件所在的目录
    for path in (
        args.pass1_bam, args.pass1_fq1, args.pass2_bam,
        args.pass2_unmapped1, args.pass3a_bam, args.pass3b_bam,
        args.bam, log_path,
    ):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

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