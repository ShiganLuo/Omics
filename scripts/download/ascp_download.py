#!/usr/bin/env python3
import time
import subprocess
import logging
import argparse
from pathlib import Path
from typing import List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd



# ============================================================
# logging 初始化
# ============================================================

def setup_logger(log_file: Path) -> logging.Logger:
    logger = logging.getLogger("ENA_DOWNLOAD")
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
# ENA 路径构建
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
# ascp + gzip
# ============================================================

def ascp_download(remote: str, dest: Path, key: Path) -> bool:
    cmd = [
        "ascp", "-k", "1", "-T", "-l", "200m",
        "-P", "33001",
        "--file-checksum=md5",
        "--overwrite=always",
        "-i", str(key),
        remote, str(dest)
    ]
    return subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    ).returncode == 0


def gzip_test(path: Path) -> bool:
    if not path.exists():
        return False
    return subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    ).returncode == 0


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

def ena_download_spin(
    srr_id: str,
    library_type: str,
    dest: Path,
    key: Path,
    logger: logging.Logger,
    sleep_base: int,
    sleep_max: int
):
    paths = build_ena_paths(srr_id)
    local = [dest / Path(p).name for p in paths]
    dest.mkdir(parents=True, exist_ok=True)

    logger.info(f"{srr_id} {library_type} 开始下载")

    if library_type == "PAIRED":

        def paired_try():
            if (ascp_download(paths[0], dest, key) and
                    ascp_download(paths[1], dest, key)):
                return gzip_test(local[0]) and gzip_test(local[1])

            if ascp_download(paths[2], dest, key):
                return gzip_test(local[2])

            return False

        spin_until_success(
            paired_try, f"{srr_id} PAIRED", logger,
            sleep_base, sleep_max
        )

    else:  # SINGLE

        def single_try():
            for i in (2, 0, 1):
                if ascp_download(paths[i], dest, key):
                    if gzip_test(local[i]):
                        return True
            return False

        spin_until_success(
            single_try, f"{srr_id} SINGLE", logger,
            sleep_base, sleep_max
        )


# ============================================================
# SRR 解析（关键：灵活）
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
        # ---------- Read meta table ----------
        try:
            df = pd.read_csv(
                args.meta,
                sep=None,          # auto-detect CSV / TSV
                engine="python",
                comment="#"
            )
        except Exception as e:
            raise ValueError(f"Failed to read meta file: {args.meta}") from e

        # ---------- Validate columns ----------
        srr_col = args.srr_col_name
        lib_col = args.lib_col_name

        missing = {c for c in (srr_col, lib_col) if c not in df.columns}
        if missing:
            raise ValueError(
                f"Missing required columns in meta file: {missing}. "
                f"Available columns: {list(df.columns)}"
            )

        # ---------- Normalize & validate ----------
        df = df[[srr_col, lib_col]].dropna()

        df[lib_col] = df[lib_col].str.upper()
        invalid = df[~df[lib_col].isin({"PAIRED", "SINGLE"})]
        if not invalid.empty:
            raise ValueError(
                f"Invalid library layout values found:\n{invalid}"
            )

        # ---------- Build tasks ----------
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
    """
    Parse command-line arguments for ENA downloader.

    Meta file rules
    ---------------
    - Meta file MUST contain a header
    - Column access is strictly based on column names
    - Default column names:
        SRR     : SRR accession
        Layout  : library layout (PAIRED | SINGLE)
    """
    p = argparse.ArgumentParser("ENA downloader (spin + parallel)")

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

    # ---------- Output / logging ----------
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("-k", "--key", type=Path, required=True)
    p.add_argument("-l", "--log", type=Path, required=True)

    # ---------- Parallel / retry ----------
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--sleep-base", type=int, default=10)
    p.add_argument("--sleep-max", type=int, default=300)

    return p.parse_args()



# ============================================================
# main
# ============================================================

def main():
    args = parse_args()
    logger = setup_logger(args.log)

    tasks = load_tasks(args)
    logger.info(f"共 {len(tasks)} 个 SRR，jobs={args.jobs}")

    if args.jobs == 1:
        for srr, lib in tasks:
            ena_download_spin(
                srr, lib, args.outdir,
                args.key, logger,
                args.sleep_base, args.sleep_max
            )
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futures = [
                ex.submit(
                    ena_download_spin,
                    srr, lib,
                    args.outdir, args.key,
                    logger,
                    args.sleep_base, args.sleep_max
                )
                for srr, lib in tasks
            ]
            for _ in as_completed(futures):
                pass


if __name__ == "__main__":
    main()
