#!/usr/bin/env python3
import os
import time
import shutil
import subprocess
import logging
import argparse
import sys
import urllib.parse
from pathlib import Path
from typing import List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

_SRC_DIR = Path(__file__).resolve().parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

try:
    from common.util.SepUtil import detect_delimiter
except ImportError:
    # 为了保证代码独立运行不报错，提供一个 fallback 函数
    def detect_delimiter(file_path):
        return "\t" if str(file_path).endswith((".tsv", ".txt")) else ","

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
# ENA 路径构建 (ascp & globus)
# ============================================================

def build_ena_paths(srr_id: str) -> List[str]:
    """为 ascp 构建远端路径"""
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


def build_ena_globus_dir(srr_id: str) -> str:
    """为 globus 构建 ENA 的根目录路径"""
    n = len(srr_id)
    if n == 11:
        x6 = srr_id[:6]
        x2 = f"0{srr_id[-2:]}"
        return f"/vol1/fastq/{x6}/{x2}/{srr_id}"
    elif n == 10:
        x6 = srr_id[:6]
        x2 = f"00{srr_id[-1]}"
        return f"/vol1/fastq/{x6}/{x2}/{srr_id}"
    elif n == 9:
        x6 = srr_id[:6]
        return f"/vol1/fastq/{x6}/{srr_id}"
    else:
        raise ValueError(f"非法 SRR ID: {srr_id}")


def build_ena_http_base(srr_id: str) -> str:
    """构建 ENA HTTP 下载的 base URL（不含文件名）"""
    n = len(srr_id)
    if n == 11:
        x6 = srr_id[:6]
        x2 = f"0{srr_id[-2:]}"
        return f"http://ftp.sra.ebi.ac.uk/vol1/fastq/{x6}/{x2}/{srr_id}"
    elif n == 10:
        x6 = srr_id[:6]
        x2 = f"00{srr_id[-1]}"
        return f"http://ftp.sra.ebi.ac.uk/vol1/fastq/{x6}/{x2}/{srr_id}"
    elif n == 9:
        x6 = srr_id[:6]
        return f"http://ftp.sra.ebi.ac.uk/vol1/fastq/{x6}/{srr_id}"
    else:
        raise ValueError(f"非法 SRR ID: {srr_id}")


def build_ena_http_urls(srr_id: str) -> List[str]:
    """构建 ENA HTTP 下载 URL 列表，顺序: _1, _2, .fastq.gz"""
    base = build_ena_http_base(srr_id)
    return [
        f"{base}/{srr_id}_1.fastq.gz",
        f"{base}/{srr_id}_2.fastq.gz",
        f"{base}/{srr_id}.fastq.gz",
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
# Globus SDK 下载
# ============================================================

def get_globus_client(client_id: str, token: Optional[str]) -> "globus_sdk.TransferClient":
    try:
        import globus_sdk
    except ImportError:
        raise ImportError("请先安装 globus-sdk: pip install globus-sdk")

    if token:
        authorizer = globus_sdk.AccessTokenAuthorizer(token)
    else:
        # Native App Auth 流程 (如果没有 token，走网页授权)
        client = globus_sdk.NativeAppAuthClient(client_id)
        
        # 【修改处】明确请求 Globus Transfer API 的权限范围
        client.oauth2_start_flow(
            requested_scopes=[globus_sdk.TransferClient.scopes.all]
        )
        
        authorize_url = client.oauth2_get_authorize_url()
        print(f"\n================ Globus 认证 ================\n"
              f"请在浏览器中打开此URL并登录:\n\n{authorize_url}\n")
        auth_code = input("请输入获取到的 Authorization Code: ").strip()
        token_response = client.oauth2_exchange_code_for_tokens(auth_code)
        
        # 提取 Transfer API 的 access_token
        transfer_token = token_response.by_resource_server["transfer.api.globus.org"]["access_token"]
        authorizer = globus_sdk.AccessTokenAuthorizer(transfer_token)
        print("Globus 鉴权成功！\n===========================================\n")

    return globus_sdk.TransferClient(authorizer=authorizer)


def globus_download_single_srr(
    tc: "globus_sdk.TransferClient",
    source_ep: str,
    dest_ep: str,
    srr_id: str,
    library_type: str,
    dest: Path,
    logger: logging.Logger,
) -> bool:
    import globus_sdk
    
    base_dir = build_ena_globus_dir(srr_id)
    dest.mkdir(parents=True, exist_ok=True)
    
    # 1. 尝试获取该 SRR 目录下的文件列表
    try:
        ls_res = tc.operation_ls(source_ep, path=base_dir)
        files_present = [
            item["name"] for item in ls_res 
            if item["type"] == "file" and item["name"].startswith(srr_id)
        ]
    except globus_sdk.TransferAPIError as e:
        logger.error(f"[Globus] 无法列出目录 {base_dir}: {e}")
        return False

    if not files_present:
        logger.warning(f"[Globus] 没有在远端找到 {srr_id} 相关的 fastq 文件")
        return False

    # 2. 根据 Layout 判断需要下载哪些文件
    target_files = []
    if library_type == "PAIRED":
        if f"{srr_id}_1.fastq.gz" in files_present and f"{srr_id}_2.fastq.gz" in files_present:
            target_files = [f"{srr_id}_1.fastq.gz", f"{srr_id}_2.fastq.gz"]
        elif f"{srr_id}.fastq.gz" in files_present:
            target_files = [f"{srr_id}.fastq.gz"]
        else:
            return False
    else:  # SINGLE
        if f"{srr_id}.fastq.gz" in files_present:
            target_files = [f"{srr_id}.fastq.gz"]
        elif f"{srr_id}_1.fastq.gz" in files_present:
            # 有时 SINGLE 也会被命名为 _1
            target_files = [f"{srr_id}_1.fastq.gz"]
        else:
            return False

    # 3. 创建 Transfer Data
    tdata = globus_sdk.TransferData(
        source_ep, dest_ep,
        label=f"SRA Download {srr_id}",
        sync_level="checksum" # 开启校验
    )
    
    for f in target_files:
        source_path = f"{base_dir}/{f}"
        dest_path = f"{str(dest.resolve())}/{f}"
        tdata.add_item(source_path, dest_path)

    # 4. 提交并等待 Task 完成
    try:
        res = tc.submit_transfer(tdata)
        task_id = res["task_id"]
        logger.info(f"[Globus] 已提交 Transfer Task (ID: {task_id}) 正在等待完成...")
        
        while not tc.task_wait(task_id, timeout=60):
            task_status = tc.get_task(task_id)["status"]
            logger.info(f"[Globus] {srr_id} task 状态: {task_status}...")

        task = tc.get_task(task_id)
        if task["status"] == "SUCCEEDED":
            # 二次本地校验
            for f in target_files:
                if not gzip_test(dest / f):
                    logger.warning(f"[Globus] {f} 下载成功，但 gzip 校验失败")
                    return False
            return True
        else:
            logger.error(f"[Globus] {srr_id} 下载失败，状态: {task['status']}")
            return False
            
    except globus_sdk.TransferAPIError as e:
        logger.error(f"[Globus] 传输请求错误: {e}")
        return False


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
# aria2c 下载 (ENA HTTP 多线程)
# ============================================================

def aria2c_download_single(
    url: str,
    dest_dir: Path,
    logger: Optional[logging.Logger] = None,
    connections: int = 8,
    split: int = 8,
    min_split_size: str = "1M",
    timeout: int = 600,
) -> bool:
    """使用 aria2c 多线程下载单个文件，实时解析并输出下载速率"""
    filename = url.rsplit("/", 1)[-1]
    cmd = [
        "aria2c",
        "-x", str(connections),
        "-s", str(split),
        "-k", min_split_size,
        "--timeout", str(timeout),
        "--retry-wait", "5",
        "--max-tries", "5",
        "--continue=true",
        "--auto-file-renaming=false",
        "--allow-overwrite=true",
        "--summary-interval=5",
        "--console-log-level=warn",
        "-d", str(dest_dir),
        url,
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert proc.stdout is not None
    last_log_time = 0.0
    last_progress = ""
    for line in proc.stdout:
        line = line.strip()
        if not line:
            continue
        # aria2c 进度行格式: [#gid MiB/TotalMiB(%) CN:xx DL:xxMiB ETA:xxs]
        if line.startswith("[#") and "DL:" in line:
            last_progress = line
            now = time.time()
            if now - last_log_time >= 5.0:
                if logger:
                    logger.info(f"[aria2c] {filename}  {line}")
                last_log_time = now
        elif "Download Results" in line or "Status Legend" in line:
            # 输出最后一条进度（距上次日志 >= 2s 才输出，避免重复）
            if last_progress and logger and time.time() - last_log_time >= 2.0:
                logger.info(f"[aria2c] {filename}  {last_progress}")
            break
    proc.wait()
    rc = proc.returncode
    if rc != 0:
        if logger:
            logger.warning(f"[aria2c] {filename} 退出码 {rc}")
        return False
    return True


def aria2c_download_single_srr(
    srr_id: str,
    library_type: str,
    dest: Path,
    logger: logging.Logger,
) -> bool:
    """从 ENA HTTP 使用 aria2c 多线程下载单个 SRR 的 fastq.gz"""
    urls = build_ena_http_urls(srr_id)
    local = [dest / Path(urllib.parse.urlparse(u).path).name for u in urls]
    dest.mkdir(parents=True, exist_ok=True)

    if library_type == "PAIRED":
        # 优先尝试 paired (_1 + _2)
        ok1 = aria2c_download_single(urls[0], dest, logger)
        if ok1 and gzip_test(local[0]):
            ok2 = aria2c_download_single(urls[1], dest, logger)
            if ok2 and gzip_test(local[1]):
                return True
            logger.warning(f"[aria2c] {srr_id}_2 下载失败，尝试 single-end fallback")
            local[1].unlink(missing_ok=True)
        else:
            logger.warning(f"[aria2c] {srr_id}_1 下载失败，尝试 single-end fallback")
            local[0].unlink(missing_ok=True)

        # fallback: single-end (.fastq.gz)
        if aria2c_download_single(urls[2], dest, logger):
            if gzip_test(local[2]):
                # 清理可能残留的 _1/_2
                local[0].unlink(missing_ok=True)
                local[1].unlink(missing_ok=True)
                return True
        return False
    else:  # SINGLE
        for i in (2, 0, 1):
            if aria2c_download_single(urls[i], dest, logger):
                if gzip_test(local[i]):
                    # 清理不需要的文件
                    for j in range(3):
                        if j != i:
                            local[j].unlink(missing_ok=True)
                    return True
                else:
                    logger.warning(f"[aria2c] {local[i].name} gzip 校验失败")
                    local[i].unlink(missing_ok=True)
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
# 自旋重试逻辑
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
# 单 SRR 下载（主调度）
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
    globus_tc: Optional["globus_sdk.TransferClient"] = None,
    globus_src_ep: Optional[str] = None,
    globus_dest_ep: Optional[str] = None,
):
    logger.info(f"{srr_id} {library_type} 开始下载 (method={method})")

    if method == "ascp":
        if not key:
            raise ValueError("ascp 方法需要 --key 参数")
        try_func = lambda: ena_download_single_srr(srr_id, library_type, dest, key, logger)
    elif method == "aria2c":
        try_func = lambda: aria2c_download_single_srr(srr_id, library_type, dest, logger)
    elif method == "globus":
        if not globus_tc or not globus_dest_ep:
            raise ValueError("globus 方法缺失 transfer_client 或 destination endpoint")
        try_func = lambda: globus_download_single_srr(globus_tc, globus_src_ep, globus_dest_ep, srr_id, library_type, dest, logger)
    else:  # sra
        try_func = lambda: sra_download_single_srr(srr_id, library_type, dest, logger)

    spin_until_success(try_func, f"{srr_id} {library_type}", logger, sleep_base, sleep_max)


# ============================================================
# SRR 解析
# ============================================================

def load_tasks(args) -> List[Tuple[str, str]]:
    tasks: List[Tuple[str, str]] = []

    if args.meta:
        sep = detect_delimiter(str(args.meta))
        read_sep = r"\s+" if sep == "whitespace" else sep
        try:
            df = pd.read_csv(args.meta, sep=read_sep, comment="#")
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
            raise ValueError(f"Invalid library layout values found:\n{invalid}")

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
        description="Download SRA data via ascp (ENA), aria2c (ENA HTTP), prefetch+fasterq-dump (NCBI), or globus (ENA Endpoint)"
    )

    # ---------- Input modes ----------
    p.add_argument("--srr-id", help="Single SRR accession")
    p.add_argument("--srr-list", type=Path, help="File with one SRR accession per line")
    p.add_argument("--meta", type=Path, help="Meta table with header (recommended)")

    # ---------- Meta column names ----------
    p.add_argument("--srr-col-name", default="SRR", help="SRR column name in meta file (default: SRR)")
    p.add_argument("--lib-col-name", default="Layout", help="Library layout column name (default: Layout)")

    # ---------- Library type ----------
    p.add_argument("-t", "--library-type", choices=["PAIRED", "SINGLE"], help="Library type")

    # ---------- Download method ----------
    p.add_argument("-m", "--method", choices=["ascp", "aria2c", "sra", "globus"], default="sra", 
                   help="Download method: ascp (ENA fasp), aria2c (ENA HTTP multi-thread), sra (NCBI prefetch, default), globus")

    # ---------- Output / logging ----------
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("-l", "--log", type=Path, required=True)

    # ---------- ascp param ----------
    p.add_argument("-k", "--key", type=Path, help="Aspera key file (required for --method ascp)")

    # ---------- Globus SDK params ----------
    # 默认 UUID: 61338d24-54d5-408f-a10d-66c06b59f6d2 为 Globus 官方 Tutorial App
    p.add_argument("--globus-client-id", default="61338d24-54d5-408f-a10d-66c06b59f6d2", 
                   help="Globus Client ID (Optional)")
    p.add_argument("--globus-token", help="Globus Access Token (可选，不提供则走命令行交互网页授权)")
    # 默认 UUID: 1d547d2a-e85d-11e8-963d-0a1d4c5c824a 为 ENA Public
    p.add_argument("--globus-source-ep", default="47772002-3e5b-4fd3-b97c-18cee38d6df2", 
                   help="Globus 源端 Endpoint UUID (默认: ENA public)")
    p.add_argument("--globus-dest-ep", default="76fbc1da-94c8-11f1-9c2b-02ce27bde401",
                   help="Globus 目标端 Endpoint UUID (需要提供，即当前机器的 Endpoint UUID)")

    # ---------- Parallel / retry ----------
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--sleep-base", type=int, default=10)
    p.add_argument("--sleep-max", type=int, default=300)

    args = p.parse_args()

    if args.method == "ascp" and not args.key:
        p.error("--method ascp 需要指定 --key 参数")
        
    if args.method == "globus" and not args.globus_dest_ep:
        p.error("--method globus 需要提供 --globus-dest-ep (你的本地或集群 Globus Endpoint UUID)")

    return args


# ============================================================
# main
# ============================================================

def main():
    args = parse_args()
    logger = setup_logger(args.log)

    tasks = load_tasks(args)
    logger.info(f"共 {len(tasks)} 个 SRR，method={args.method}，jobs={args.jobs}")

    # 若使用 Globus 则提前初始化客户端避免重复登录
    globus_tc = None
    if args.method == "globus":
        globus_tc = get_globus_client(args.globus_client_id, args.globus_token)

    download_fn = lambda srr, lib: download_spin(
        srr, lib, args.outdir, args.method,
        logger, args.sleep_base, args.sleep_max,
        key=args.key,
        globus_tc=globus_tc,
        globus_src_ep=args.globus_source_ep,
        globus_dest_ep=args.globus_dest_ep
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