#!/usr/bin/env python3
"""Call FIRE peaks with `ft call-peaks`, building the shuffled null for FDR mode.

Single entry point for the ft_call_peaks rule so its receipt stays ONE
`python call_peaks.py ...` command (modules.md: 支撑逻辑收敛到 bin/ 脚本 —
规则只负责校验和调用). Phases, in order:

1. Guard: count aligned reads; a BAM with none gets header-only outputs for
   BOTH declared files (the rule registers them unconditionally).
2. fdr mode: stream fiber spans to a scratch BED (`samtools view -F 0x904`,
   SAM fields -> BED4), then shuffle them with shuffle_bed.py to build the
   null that `ft call-peaks --shuffled` scores against.
3. Run `ft call-peaks` (fdr or fire_frac flag set).
4. fire_frac mode: write a placeholder FDR table (the output is declared).
5. Remove the scratch dir on success (kept on failure for debugging).

Usage:
    python call_peaks.py --bam in.bam --out-peaks peaks.bed \
        --fdr-table fdr.tsv --mode fdr --fai genome.fa.fai --threads 8 \
        --tmp-dir /path/tmp --sample-id s1
"""
from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence, TextIO


def setup_logger(
    logger_name: str,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    stream: Optional[TextIO] = None,
) -> logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


logger = setup_logger("CallPeaks", stream=sys.stdout)

# Header-only placeholder contents for guard / fire_frac branches.
PEAKS_HEADER = "#chrom\tstart\tend\tname\tscore\tstrand\n"
FDR_TABLE_HEADER = "#sample\tsd_cov\tthreshold\tmax_fdr\tn_peaks\n"
FDR_PLACEHOLDER = "#FDR table not computed (fire_frac mode)\n"


def count_aligned_reads(bam: str, samtools: str) -> int:
    """Count aligned reads (`samtools view -c -F 4`).

    Raises RuntimeError when samtools is missing or exits nonzero, naming
    the BAM so the failing input is identifiable.
    """
    cmd = [samtools, "view", "-c", "-F", "4", bam]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        raise RuntimeError(f"samtools not found on PATH (tried '{samtools}')")
    if r.returncode != 0:
        raise RuntimeError(
            f"samtools view -c failed with exit code {r.returncode} on {bam}: "
            f"{r.stderr.strip()}"
        )
    return int(r.stdout.strip())


def write_placeholder(path: str, content: str) -> None:
    """Write a header-only placeholder so the declared output always exists."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(content, encoding="utf-8")


def extract_fiber_spans(bam: str, out_bed: str, samtools: str) -> int:
    """Stream fiber alignment spans from SAM to a BED4 scratch file.

    Equivalent to the previous `samtools view -F 0x904 | awk` pipeline:
    chrom=RNAME, start=POS-1, end=start+len(SEQ), name=QNAME. Lines are
    streamed (fiber BAMs run to millions of reads); reads without SEQ
    ('*') are skipped and counted. Returns the number of spans written.
    """
    n_written, n_no_seq = 0, 0
    out_path = Path(out_bed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.Popen(
        [samtools, "view", "-F", "0x904", bam],
        stdout=subprocess.PIPE,
        text=True,
    )
    with open(out_path, "w") as out:
        assert proc.stdout is not None  # PIPE above guarantees a pipe stream
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or f[2] == "*" or f[9] == "*":
                n_no_seq += 1
                continue
            start = int(f[3]) - 1
            out.write(f"{f[2]}\t{start}\t{start + len(f[9])}\t{f[0]}\n")
            n_written += 1
    if proc.wait() != 0:
        raise RuntimeError(f"samtools view failed on {bam} (exit {proc.returncode})")
    if n_no_seq:
        logger.info(f"Skipped {n_no_seq} reads without SEQ/RNAME")
    return n_written


def build_shuffled_null(
    fibers_bed: str, fai: str, shuffled_bed: str, shuffle_script: str, seed: int
) -> None:
    """Shuffle fiber spans genome-wide via the sibling shuffle_bed.py CLI.

    Uses sys.executable so the shuffle runs under the same interpreter/env
    as this script (inside the rule's container that is fibertools.yaml).
    """
    cmd = [
        sys.executable, shuffle_script,
        "-i", fibers_bed,
        "-g", fai,
        "-o", shuffled_bed,
        "--seed", str(seed),
    ]
    logger.info("Building shuffled null: " + " ".join(cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0 or not Path(shuffled_bed).exists():
        raise RuntimeError(
            f"shuffle_bed.py failed (exit {r.returncode}); check the fai "
            f"path '{fai}' matches the BAM's reference (set genome.fai in "
            f"the workflow config JSON)."
        )


def call_ft_peaks(
    ft: str,
    bam: str,
    out_peaks: str,
    fdr_table: str,
    threads: int,
    mode: str,
    sd_cov: float,
    max_fdr: float,
    min_fire_frac: float,
    haps: bool,
    shuffled_bed: Optional[str],
) -> None:
    """Run `ft call-peaks` and verify the peak BED landed.

    Flag set mirrors the upstream CLI: fdr mode scores real coverage
    against --shuffled with --max-fdr and writes --fdr-table-out; fire_frac
    mode thresholds on --min-fire-frac instead (the FDR table placeholder
    is written by the caller).
    """
    cmd = [
        ft, "call-peaks",
        "-t", str(threads),
        "-o", out_peaks,
        "--sd-cov", str(sd_cov),
    ]
    if mode == "fdr":
        cmd += [
            "--max-fdr", str(max_fdr),
            "--shuffled", str(shuffled_bed),
            "--fdr-table-out", fdr_table,
        ]
    else:
        cmd += ["--min-fire-frac", str(min_fire_frac)]
    if haps:
        cmd.append("--haps")
    cmd.append(bam)
    logger.info("Calling FIRE peaks: " + " ".join(cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise RuntimeError(
            f"ft call-peaks failed with exit code {r.returncode}; the tool "
            f"log above names the failing stage."
        )
    if not Path(out_peaks).exists():
        raise RuntimeError(
            f"ft call-peaks exited 0 but the peak BED is missing: {out_peaks}"
        )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """CLI: every knob of the rule's former inline shell, parameterized."""
    p = argparse.ArgumentParser(
        prog="call_peaks.py",
        description="FIRE peak calling with FDR null (wrapper for ft call-peaks).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam", required=True, help="aligned Fiber-seq BAM (m6a/nuc/FIRE called)")
    p.add_argument("--out-peaks", required=True, help="output FIRE peaks BED")
    p.add_argument("--fdr-table", required=True, help="output FDR threshold table")
    p.add_argument("--mode", choices=["fdr", "fire_frac"], default="fdr",
                   help="ft call-peaks scoring mode")
    p.add_argument("--fai", default=None,
                   help="chrom sizes (.fai) for the shuffled null (fdr mode)")
    p.add_argument("--shuffle-script", type=Path,
                   default=Path(__file__).resolve().parent / "shuffle_bed.py",
                   help="path to the sibling shuffle_bed.py")
    p.add_argument("--ft", default="ft", help="fibertools-rs executable")
    p.add_argument("--samtools", default="samtools", help="samtools executable")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--sd-cov", type=float, default=5.0,
                   help="coverage SD threshold above the chromosome mean")
    p.add_argument("--max-fdr", type=float, default=0.05)
    p.add_argument("--min-fire-frac", type=float, default=0.1)
    p.add_argument("--haps", action="store_true", help="call peaks per haplotype")
    p.add_argument("--seed", type=int, default=42, help="shuffle seed")
    p.add_argument("--tmp-dir", type=Path, default=Path("/tmp"),
                   help="workflow tmp dir for scratch (null intermediates)")
    p.add_argument("--sample-id", default="", help="sample tag for scratch naming/logs")
    p.add_argument("--log", type=Path, default=None, help="log file (default stdout only)")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Orchestrate the phases; every exit path leaves both outputs present."""
    global logger
    args = parse_args(argv)
    if args.log:
        logger = setup_logger("CallPeaks", log_file=str(args.log), stream=sys.stdout)
    logger.info(f"Command: {' '.join(sys.argv)}")
    try:
        if not Path(args.bam).exists():
            logger.error(f"input BAM not found: {args.bam}")
            return 1
        if args.mode == "fdr":
            for f, key in ((args.fai, "genome.fai"), (args.shuffle_script, "shuffle_script")):
                if not f or not Path(f).exists():
                    logger.error(
                        f"fdr mode requires an existing {key} "
                        f"(got: {f}); set it in the workflow config JSON."
                    )
                    return 1

        # scratch dir: pre-clean stale runs, recreate fresh (modules.md
        # 输出目录卫生 — never inside the results dir)
        scratch = Path(args.tmp_dir) / "call_peaks" / (args.sample_id or "sample")
        shutil.rmtree(scratch, ignore_errors=True)
        scratch.mkdir(parents=True, exist_ok=True)

        # Phase 1: aligned-read guard — both outputs as header-only files
        n_aligned = count_aligned_reads(args.bam, args.samtools)
        logger.info(f"Aligned reads: {n_aligned}")
        if n_aligned == 0:
            logger.info("BAM has no aligned reads — writing empty outputs")
            write_placeholder(args.out_peaks, PEAKS_HEADER)
            write_placeholder(args.fdr_table, FDR_TABLE_HEADER)
            shutil.rmtree(scratch, ignore_errors=True)
            return 0

        # Phase 2: shuffled null (fdr mode only)
        shuffled_bed: Optional[str] = None
        if args.mode == "fdr":
            fibers_bed = str(scratch / f"{args.sample_id or 'sample'}.fibers.bed")
            shuffled_bed = str(scratch / f"{args.sample_id or 'sample'}.fibers.shuffled.bed")
            n_spans = extract_fiber_spans(args.bam, fibers_bed, args.samtools)
            logger.info(f"Extracted {n_spans} fiber spans -> {fibers_bed}")
            build_shuffled_null(fibers_bed, args.fai, shuffled_bed,
                                str(args.shuffle_script), args.seed)

        # Phase 3: the tool itself
        call_ft_peaks(
            args.ft, args.bam, args.out_peaks, args.fdr_table,
            args.threads, args.mode, args.sd_cov, args.max_fdr,
            args.min_fire_frac, args.haps, shuffled_bed,
        )

        # Phase 4: fire_frac has no FDR table — keep the declared output
        if args.mode != "fdr":
            write_placeholder(args.fdr_table, FDR_PLACEHOLDER)

        shutil.rmtree(scratch, ignore_errors=True)
        logger.info(f"FIRE peak calling completed for sample {args.sample_id}")
        return 0
    except Exception as e:
        logger.error(f"call_peaks failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
