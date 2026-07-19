#!/usr/bin/env python3
import os
import time
import shutil
import subprocess
import logging
import argparse
from pathlib import Path
from typing import List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd



# ============================================================
# logging 初始化
# ============================================================

def setup_logger(log_file: Path) -> logging.Logger:
    logger = logging.getLogger("SRA_DOWNLOAD")
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    fh = logging.FileHandler(log_file)
    fh.setFormatter(formatter)

    sh = logging.StreamHandler()
    sh.setFormatter(formatter)

    if not logger.handlers:
        logger.addHandler(fh)
        logger.addHandler(sh)

    return logger


# ============================================================
# ENA 路径构建 (ascp)
# ============================================================

def build_ena_paths(srr_id: str) -> List[str]:
    n = len(srr_id)

    if n == 11:
        x6 = srr_id[:6]
        x2 = f"0{srr_id[-2:]}"
        base = f"/vol1/fastq/{x6}/{x2}/{srr_id}"
    elif n == 10:
        x6 = srr_id[:6]
        x2 = f"00{srr_id[-1]}"
        base = f"/vol1/fastq/{x6}/{x2}/{srr_id}"
    elif n == 9:
        x6 = srr_id[:6]
        base = f"/vol1/fastq/{x6}/{srr_id}"
    else:
        raise ValueError(f"非法 SRR ID: {srr_id}")

    prefix = "era-fasp@fasp.sra.ebi.ac.uk:"
    return [
        f"{prefix}{base}/{srr_id}_1.fastq.gz",
        f"{prefix}{base}/{srr_id}_2.fastq.gz",
        f"{prefix}{base}/{srr_id}.fastq.gz",
    ]


# ============================================================
# 公共工具
# ============================================================

def gzip_test(path: Path) -> bool:
    if not path.exists():
        return False
    return subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    ).returncode == 0


def _log_subprocess_error(logger: Optional[logging.Logger], label: str, result: subprocess.CompletedProcess):
    msg = result.stderr.strip() or f"exit code {result.returncode}"
    if logger:
        logger.warning(f"[{label}] {msg}")


# ============================================================
# ascp 下载
# ============================================================

def ascp_download(remote: str, dest: Path, key: Path, logger: Optional[logging.Logger] = None) -> bool:
    cmd = [
        "ascp", "-k", "1", "-T", "-l", "200m",
        "-P", "33001",
        "--file-checksum=md5",
        "--overwrite=always",
        "-i", str(key),
        remote, str(dest)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        _log_subprocess_error(logger, "ascp", result)
        return False
    return True


def ena_download_single_srr(
    srr_id: str,
    library_type: str,
    dest: Path,
    key: Path,
    logger: logging.Logger,
):
    """ascp 单次尝试（不含自旋），返回成功/失败。"""
    paths = build_ena_paths(srr_id)
    local = [dest / Path(p).name for p in paths]
    dest.mkdir(parents=True, exist_ok=True)

    if library_type == "PAIRED":
        if (ascp_download(paths[0], dest, key, logger) and
                ascp_download(paths[1], dest, key, logger)):
            return gzip_test(local[0]) and gzip_test(local[1])
        if ascp_download(paths[2], dest, key, logger):
            return gzip_test(local[2])
        return False
    else:  # SINGLE
        for i in (2, 0, 1):
            if ascp_download(paths[i], dest, key, logger):
                if gzip_test(local[i]):
                    return True
        return False


# ============================================================
# SRA Toolkit 下载 (prefetch + fasterq-dump)
# ============================================================

def sra_download_single_srr(
    srr_id: str,
    library_type: str,
    dest: Path,
    logger: logging.Logger,
) -> bool:
    """prefetch 下载 .sra → fasterq-dump 转 fastq → gzip 压缩。"""
    dest.mkdir(parents=True, exist_ok=True)

    # --- 1. prefetch ---
    sra_dir = dest / srr_id
    sra_file = sra_dir / f"{srr_id}.sra"
    prefetch_cmd = ["prefetch", srr_id, "-O", str(dest)]
    result = subprocess.run(prefetch_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        _log_subprocess_error(logger, "prefetch", result)
        shutil.rmtree(sra_dir, ignore_errors=True)
        return False

    if not sra_file.exists():
        # prefetch 有时直接放到 dest/ 下而非 dest/SRR* 下
        alt = dest / f"{srr_id}.sra"
        if alt.exists():
            sra_file = alt
        else:
            if logger:
                logger.warning(f"[prefetch] {srr_id} 完成但 .sra 文件不存在: {sra_file}")
            return False

    # --- 2. fasterq-dump ---
    fq_cmd = [
        "fasterq-dump", str(sra_file),
        "--outdir", str(dest),
        "--split-files",
        "--threads", "4",
    ]
    result = subprocess.run(fq_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        _log_subprocess_error(logger, "fasterq-dump", result)
        sra_file.unlink(missing_ok=True)
        shutil.rmtree(sra_dir, ignore_errors=True)
        return False

    # --- 3. gzip ---
    generated = list(dest.glob(f"{srr_id}*.fastq"))
    if not generated:
        if logger:
            logger.warning(f"[fasterq-dump] {srr_id} 未生成 .fastq 文件")
        sra_file.unlink(missing_ok=True)
        shutil.rmtree(sra_dir, ignore_errors=True)
        return False

    for fq in generated:
        gz_cmd = ["gzip", "-f", str(fq)]
        result = subprocess.run(gz_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            _log_subprocess_error(logger, "gzip", result)
            sra_file.unlink(missing_ok=True)
            shutil.rmtree(sra_dir, ignore_errors=True)
            return False

    # --- 4. 清理 .sra ---
    sra_file.unlink(missing_ok=True)
    if sra_dir.exists() and not any(sra_dir.iterdir()):
        sra_dir.rmdir()

    # --- 5. 验证 ---
    ok = True
    for fq_gz in dest.glob(f"{srr_id}*.fastq.gz"):
        if not gzip_test(fq_gz):
            if logger:
                logger.warning(f"[gzip] {fq_gz.name} 校验失败")
            ok = False
    return ok


# ============================================================
# 自旋
# ============================================================

def spin_until_success(try_func, desc, logger, sleep_base, sleep_max):
    attempt = 1
    sleep_time = sleep_base

    while True:
        logger.info(f"[Attempt {attempt}] {desc}")

        if try_func():
            logger.info(f"[SUCCESS] {desc}")
            return

        logger.warning(f"[FAIL] {desc}，{sleep_time}s 后重试")
        time.sleep(sleep_time)
        sleep_time = min(sleep_time * 2, sleep_max)
        attempt += 1


# ============================================================
# 单 SRR 下载（内部自旋）
# ============================================================

def download_spin(
    srr_id: str,
    library_type: str,
    dest: Path,
    method: str,
    logger: logging.Logger,
    sleep_base: int,
    sleep_max: int,
    key: Optional[Path] = None,
):
    logger.info(f"{srr_id} {library_type} 开始下载 (method={method})")

    if method == "ascp":
        if not key:
            raise ValueError("ascp 方法需要 --key 参数")
        try_func = lambda: ena_download_single_srr(srr_id, library_type, dest, key, logger)
    else:  # sra
        try_func = lambda: sra_download_single_srr(srr_id, library_type, dest, logger)

    spin_until_success(try_func, f"{srr_id} {library_type}", logger, sleep_base, sleep_max)


# ============================================================
# SRR 解析
# ============================================================

def load_tasks(args) -> List[Tuple[str, str]]:
    """
    Load (SRR, library_type) tasks.

    Meta mode (recommended)
    -----------------------
    - Meta file MUST contain a header
    - Supports CSV / TSV automatically
    - Column access is strictly based on column names
    - Default columns:
        SRR     : SRR accession
        Layout  : PAIRED | SINGLE
    """
    tasks: List[Tuple[str, str]] = []

    if args.meta:
        try:
            df = pd.read_csv(
                args.meta,
                sep=None,
                engine="python",
                comment="#"
            )
        except Exception as e:
            raise ValueError(f"Failed to read meta file: {args.meta}") from e

        srr_col = args.srr_col_name
        lib_col = args.lib_col_name

        missing = {c for c in (srr_col, lib_col) if c not in df.columns}
        if missing:
            raise ValueError(
                f"Missing required columns in meta file: {missing}. "
                f"Available columns: {list(df.columns)}"
            )

        df = df[[srr_col, lib_col]].dropna()
        df[lib_col] = df[lib_col].str.upper()
        invalid = df[~df[lib_col].isin({"PAIRED", "SINGLE"})]
        if not invalid.empty:
            raise ValueError(
                f"Invalid library layout values found:\n{invalid}"
            )

        tasks = list(df.itertuples(index=False, name=None))

    elif args.srr_list:
        for line in args.srr_list.open():
            if line.strip():
                tasks.append((line.strip(), args.library_type))

    else:
        tasks.append((args.srr_id, args.library_type))

    return tasks



# ============================================================
# argparse
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        "SRA fastq downloader",
        description="Download SRA data via ascp (ENA) or prefetch+fasterq-dump (NCBI)"
    )

    # ---------- Input modes ----------
    p.add_argument("--srr-id",
                   help="Single SRR accession")

    p.add_argument("--srr-list", type=Path,
                   help="File with one SRR accession per line")

    p.add_argument("--meta", type=Path,
                   help="Meta table with header (recommended)")

    # ---------- Meta column names ----------
    p.add_argument("--srr-col-name",
                   default="SRR",
                   help="SRR column name in meta file (default: SRR)")

    p.add_argument("--lib-col-name",
                   default="Layout",
                   help="Library layout column name (default: Layout)")

    # ---------- Library type ----------
    p.add_argument("-t", "--library-type",
                   choices=["PAIRED", "SINGLE"],
                   help="Library type (used with --srr-id or --srr-list)")

    # ---------- Download method ----------
    p.add_argument("-m", "--method",
                   choices=["ascp", "sra"],
                   default="sra",
                   help="Download method: ascp (ENA Aspera) or sra (prefetch+fasterq-dump). Default: sra")

    # ---------- Output / logging ----------
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("-k", "--key", type=Path,
                   help="Aspera key file (required for --method ascp)")
    p.add_argument("-l", "--log", type=Path, required=True)

    # ---------- Parallel / retry ----------
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--sleep-base", type=int, default=10)
    p.add_argument("--sleep-max", type=int, default=300)

    args = p.parse_args()

    if args.method == "ascp" and not args.key:
        p.error("--method ascp 需要指定 --key 参数")

    return args



# ============================================================
# main
# ============================================================

def main():
    args = parse_args()
    logger = setup_logger(args.log)

    tasks = load_tasks(args)
    logger.info(f"共 {len(tasks)} 个 SRR，method={args.method}，jobs={args.jobs}")

    download_fn = lambda srr, lib: download_spin(
        srr, lib, args.outdir, args.method,
        logger, args.sleep_base, args.sleep_max,
        key=args.key,
    )

    if args.jobs == 1:
        for srr, lib in tasks:
            download_fn(srr, lib)
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futures = [
                ex.submit(download_fn, srr, lib)
                for srr, lib in tasks
            ]
            for _ in as_completed(futures):
                pass


if __name__ == "__main__":
    main()
