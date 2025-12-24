#!/usr/bin/env python3
import time
import subprocess
import logging
import argparse
from pathlib import Path
from typing import List


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
        raise ValueError(
            f"非法 SRR ID: {srr_id}（长度应为 9/10/11）"
        )

    prefix = "era-fasp@fasp.sra.ebi.ac.uk:"
    return [
        f"{prefix}{base}/{srr_id}_1.fastq.gz",
        f"{prefix}{base}/{srr_id}_2.fastq.gz",
        f"{prefix}{base}/{srr_id}.fastq.gz",
    ]


# ============================================================
# ascp 下载
# ============================================================

def ascp_download(
    remote: str,
    dest: Path,
    key: Path,
    logger: logging.Logger
) -> bool:
    cmd = [
        "ascp",
        "-k", "1",
        "-T",
        "-l", "200m",
        "-P", "33001",
        "--file-checksum=md5",
        "--overwrite=always",
        "-i", str(key),
        remote,
        str(dest)
    ]

    logger.debug("CMD: " + " ".join(cmd))

    res = subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    return res.returncode == 0


# ============================================================
# 自旋重试
# ============================================================

def spin_until_success(
    try_func,
    desc: str,
    logger: logging.Logger,
    sleep_base: int,
    sleep_max: int
):
    attempt = 1
    sleep_time = sleep_base

    while True:
        logger.info(f"[Attempt {attempt}] {desc}")

        if try_func():
            logger.info(f"[SUCCESS] {desc}")
            return

        logger.warning(
            f"[FAIL] {desc}，{sleep_time}s 后重试"
        )

        time.sleep(sleep_time)
        sleep_time = min(sleep_time * 2, sleep_max)
        attempt += 1


# ============================================================
# 主下载逻辑
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
    logger.info(f"{srr_id} {library_type} 开始自旋下载")

    paths = build_ena_paths(srr_id)
    dest.mkdir(parents=True, exist_ok=True)

    if library_type == "PAIRED":

        def paired_try():
            if (ascp_download(paths[0], dest, key, logger) and
                    ascp_download(paths[1], dest, key, logger)):
                return True
            return ascp_download(paths[2], dest, key, logger)

        spin_until_success(
            paired_try,
            f"{srr_id} PAIRED 下载",
            logger,
            sleep_base,
            sleep_max
        )

    elif library_type == "SINGLE":

        def single_try():
            if ascp_download(paths[2], dest, key, logger):
                return True
            if ascp_download(paths[0], dest, key, logger):
                return True
            return ascp_download(paths[1], dest, key, logger)

        spin_until_success(
            single_try,
            f"{srr_id} SINGLE 下载",
            logger,
            sleep_base,
            sleep_max
        )

    else:
        raise ValueError("Library type 必须是 PAIRED 或 SINGLE")


# ============================================================
# argparse CLI
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="ENA fastq downloader with infinite retry (spin)"
    )

    parser.add_argument(
        "-i", "--srr-id",
        required=True,
        help="SRR accession, e.g. SRR16119550"
    )

    parser.add_argument(
        "-t", "--library-type",
        required=True,
        choices=["PAIRED", "SINGLE"],
        help="Library layout: PAIRED or SINGLE"
    )

    parser.add_argument(
        "-o", "--outdir",
        required=True,
        type=Path,
        help="Output directory for fastq files"
    )

    parser.add_argument(
        "-k", "--key",
        required=True,
        type=Path,
        help="Aspera private key (asperaweb_id_dsa.openssh)"
    )

    parser.add_argument(
        "-l", "--log",
        required=True,
        type=Path,
        help="Log file path"
    )

    parser.add_argument(
        "--sleep-base",
        type=int,
        default=10,
        help="Initial sleep seconds between retries (default: 10)"
    )

    parser.add_argument(
        "--sleep-max",
        type=int,
        default=300,
        help="Maximum sleep seconds between retries (default: 300)"
    )

    return parser.parse_args()


# ============================================================
# main
# ============================================================

def main():
    args = parse_args()

    logger = setup_logger(args.log)

    ena_download_spin(
        srr_id=args.srr_id,
        library_type=args.library_type,
        dest=args.outdir,
        key=args.key,
        logger=logger,
        sleep_base=args.sleep_base,
        sleep_max=args.sleep_max
    )


if __name__ == "__main__":
    main()
