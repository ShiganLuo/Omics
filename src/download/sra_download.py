#!/usr/bin/env python3
import os
import time
import hashlib
import shutil
import subprocess
import logging
import argparse
import shlex
import sys
import threading
import urllib.parse
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.SepUtil import detect_delimiter
from src.common.util.LogUtil import setup_logger

# ============================================================
# 全局 logger：由 main() 调用 LogUtil.setup_logger 初始化，
# 各函数直接使用模块级 logger，不通过参数传递
# ============================================================
logger: logging.Logger = logging.getLogger("SRA_DOWNLOAD")

# 网络劫持中止标志：检测到下载内容被认证门户劫持时置位，剩余下载立即跳过
_ABORT_EVENT = threading.Event()


class NetworkHijackedError(RuntimeError):
    """Downloaded payload is a captive-portal login page, not real data.

    Raised when a "successful" download returns an HTML page or an empty
    body, which indicates the network gateway (e.g. Dr.COM campus portal)
    intercepted the request. Retrying cannot help until the network is
    authenticated, so the run is aborted instead of retrying.
    """


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
        return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{x6}/{x2}/{srr_id}"
    elif n == 10:
        x6 = srr_id[:6]
        x2 = f"00{srr_id[-1]}"
        return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{x6}/{x2}/{srr_id}"
    elif n == 9:
        x6 = srr_id[:6]
        return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{x6}/{srr_id}"
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
    """Return True when the file exists and passes `gzip -t`.

    Pure validation, never raises. Use :func:`verify_fastq_gzip` at download
    verification sites where hijack detection is required.

    Parameters
    ----------
    path : Path
        gzip file to test.

    Returns
    -------
    bool
        True if the file exists and is a valid gzip stream.
    """
    if not path.exists():
        return False
    return subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    ).returncode == 0


def is_html_payload(path: Path) -> bool:
    """Return True when the file content is an HTML/XML page, not binary data.

    Parameters
    ----------
    path : Path
        Downloaded file to sniff.

    Returns
    -------
    bool
        True if the first non-whitespace bytes look like an HTML document
        (e.g. a captive-portal login page injected by the network gateway).
    """
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open("rb") as fh:
        head = fh.read(1024).lstrip().lower()
    return head.startswith((b"<!doctype", b"<html", b"<?xml", b"<head", b"<script"))


def verify_fastq_gzip(path: Path) -> bool:
    """Validate a downloaded fastq.gz file, detecting captive-portal hijacks.

    Parameters
    ----------
    path : Path
        Downloaded fastq.gz file.

    Returns
    -------
    bool
        True when the file passes ``gzip -t``.

    Raises
    ------
    NetworkHijackedError
        If the downloaded content is an HTML page or a zero-byte file, i.e.
        the request was intercepted by a network authentication portal and
        retrying cannot help.
    """
    if is_html_payload(path):
        raise NetworkHijackedError(
            f"{path.name} 下载内容是 HTML 页面而非 fastq 数据，"
            f"疑似被网络认证门户(如 Dr.COM 校园网登录页)劫持: {path}"
        )
    if path.exists() and path.stat().st_size == 0:
        raise NetworkHijackedError(
            f"{path.name} 下载结果为空文件(0 字节)，疑似被网络网关拦截/重置: {path}"
        )
    if not gzip_test(path):
        size = path.stat().st_size if path.exists() else 0
        logger.warning(f"{path.name} gzip 校验失败 (size={size}B)")
        return False
    return True


def _log_subprocess_error(label: str, result: subprocess.CompletedProcess):
    msg = result.stderr.strip() or f"exit code {result.returncode}"
    logger.warning(f"[{label}] {msg}")


def md5_of_file(path: Path, chunk_size: int = 1048576) -> str:
    """Return the md5 hex digest of a file, computed in chunks.

    Parameters
    ----------
    path : Path
        File to hash.
    chunk_size : int, default=1MiB
        Read buffer size in bytes.

    Returns
    -------
    str
        Lowercase hex md5 digest.
    """
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_ena_fastq_md5(srr_id: str, connect_timeout: int = 60) -> Dict[str, str]:
    """Fetch expected md5 checksums for a run's fastq files from the ENA API.

    Parameters
    ----------
    srr_id : str
        SRA run accession.
    connect_timeout : int, default=60
        curl --connect-timeout in seconds.

    Returns
    -------
    Dict[str, str]
        Mapping from fastq file name (e.g. "SRR14664595.fastq.gz") to the
        expected md5 hex digest. Empty dict when the lookup fails.
    """
    api = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={srr_id}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
    )
    result = subprocess.run(
        ["curl", "-sL", "--connect-timeout", str(connect_timeout), api],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        logger.warning(f"[md5] ENA md5 查询失败: {srr_id}")
        return {}
    lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
    if len(lines) < 2:
        logger.warning(f"[md5] ENA 未返回 {srr_id} 的 fastq_md5")
        return {}
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    try:
        ftp_idx = header.index("fastq_ftp")
        md5_idx = header.index("fastq_md5")
    except ValueError:
        logger.warning(f"[md5] ENA 返回缺少 fastq_ftp/fastq_md5 列: {header}")
        return {}
    if ftp_idx >= len(row) or md5_idx >= len(row):
        logger.warning(f"[md5] ENA 返回行与表头不一致: {srr_id}")
        return {}
    files = row[ftp_idx].split(";")
    md5s = row[md5_idx].split(";")
    if len(files) != len(md5s):
        logger.warning(f"[md5] ENA 返回的 fastq_ftp/fastq_md5 数量不一致: {srr_id}")
        return {}
    return {
        Path(name.strip()).name: md5.strip()
        for name, md5 in zip(files, md5s)
    }


def verify_fastq_md5(path: Path, expected_md5: Optional[str]) -> bool:
    """Verify a downloaded file against the expected md5 digest.

    Parameters
    ----------
    path : Path
        Downloaded file.
    expected_md5 : Optional[str]
        Expected md5 hex digest; None skips the check with a warning.

    Returns
    -------
    bool
        True when the digest matches (or no expectation is available).
    """
    if not expected_md5:
        logger.warning(f"[md5] {path.name} 无预期 md5（ENA 查询失败或字段为空），跳过 md5 校验")
        return True
    actual = md5_of_file(path)
    if actual != expected_md5:
        logger.warning(f"[md5] {path.name} md5 不匹配: 实际 {actual} != 预期 {expected_md5}")
        return False
    logger.info(f"[md5] {path.name} md5 校验通过 ({actual})")
    return True


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

def ascp_download(remote: str, dest: Path, key: Path) -> bool:
    cmd = [
        "ascp", "-k", "1", "-T", "-l", "200m",
        "-P", "33001",
        "--file-checksum=md5",
        "--overwrite=always",
        "-i", str(key),
        remote, str(dest)
    ]
    logger.info(f"[ascp] 开始下载 {remote} -> {dest}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        _log_subprocess_error("ascp", result)
        return False
    return True


def ena_download_single_srr(
    srr_id: str,
    library_type: str,
    dest: Path,
    key: Path,
):
    paths = build_ena_paths(srr_id)
    local = [dest / Path(p).name for p in paths]
    dest.mkdir(parents=True, exist_ok=True)
    logger.info(f"[ascp] {srr_id} 候选远端路径:\n    " + "\n    ".join(paths))

    if library_type == "PAIRED":
        if (ascp_download(paths[0], dest, key) and
                ascp_download(paths[1], dest, key)):
            return verify_fastq_gzip(local[0]) and verify_fastq_gzip(local[1])
        if ascp_download(paths[2], dest, key):
            return verify_fastq_gzip(local[2])
        return False
    else:  # SINGLE
        for i in (2, 0, 1):
            if ascp_download(paths[i], dest, key):
                if verify_fastq_gzip(local[i]):
                    return True
        return False


# ============================================================
# aria2c 下载 (ENA HTTP 多线程)
# ============================================================

def aria2c_download_single(
    url: str,
    dest_dir: Path,
    connections: int = 8,
    split: int = 8,
    min_split_size: str = "1M",
    timeout: int = 600,
) -> bool:
    """使用 aria2c 多线程下载单个文件，实时解析并输出下载速率

    Parameters
    ----------
    url : str
        Complete download URL, logged verbatim before downloading.
    dest_dir : Path
        Destination directory.
    connections : int, default=8
        aria2c -x: max connections per server.
    split : int, default=8
        aria2c -s: split into N chunks.
    min_split_size : str, default="1M"
        aria2c -k: min split size.
    timeout : int, default=600
        aria2c --timeout in seconds.

    Returns
    -------
    bool
        True when aria2c exits with code 0.
    """
    filename = url.rsplit("/", 1)[-1]
    logger.info(f"[aria2c] 开始下载 {url} -> {dest_dir / filename}")
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
    detail_lines: List[str] = []
    for line in proc.stdout:
        line = line.strip()
        if not line:
            continue
        # aria2c 进度行格式: [#gid MiB/TotalMiB(%) CN:xx DL:xxMiB ETA:xxs]
        if line.startswith("[#") and "DL:" in line:
            last_progress = line
            now = time.time()
            if now - last_log_time >= 5.0:
                logger.info(f"[aria2c] {filename}  {line}")
                last_log_time = now
        elif "Download Results" in line or "Status Legend" in line:
            # 输出最后一条进度（距上次日志 >= 2s 才输出，避免重复）
            if last_progress and time.time() - last_log_time >= 2.0:
                logger.info(f"[aria2c] {filename}  {last_progress}")
            break
        else:
            # 错误/通知行（如 404、重定向），失败时输出
            detail_lines.append(line)
    proc.wait()
    rc = proc.returncode
    if rc != 0:
        for dl in detail_lines[-5:]:
            logger.warning(f"[aria2c] {dl}")
        logger.warning(f"[aria2c] {url} 退出码 {rc}")
        return False
    return True


def aria2c_download_single_srr(srr_id: str, library_type: str, dest: Path) -> bool:
    """从 ENA HTTP 使用 aria2c 多线程下载单个 SRR 的 fastq.gz

    Parameters
    ----------
    srr_id : str
        SRA run accession.
    library_type : str
        "PAIRED" or "SINGLE".
    dest : Path
        Destination directory.

    Returns
    -------
    bool
        True when a valid fastq.gz file is downloaded and verified.
    """
    urls = build_ena_http_urls(srr_id)
    local = [dest / Path(urllib.parse.urlparse(u).path).name for u in urls]
    dest.mkdir(parents=True, exist_ok=True)
    logger.info(f"[aria2c] {srr_id} 候选 URL:\n    " + "\n    ".join(urls))

    if library_type == "PAIRED":
        # 优先尝试 paired (_1 + _2)
        ok1 = aria2c_download_single(urls[0], dest)
        if ok1 and verify_fastq_gzip(local[0]):
            ok2 = aria2c_download_single(urls[1], dest)
            if ok2 and verify_fastq_gzip(local[1]):
                return True
            logger.warning(f"[aria2c] {srr_id}_2 下载失败，尝试 single-end fallback")
            local[1].unlink(missing_ok=True)
        else:
            logger.warning(f"[aria2c] {srr_id}_1 下载失败，尝试 single-end fallback")
            local[0].unlink(missing_ok=True)

        # fallback: single-end (.fastq.gz)
        if aria2c_download_single(urls[2], dest):
            if verify_fastq_gzip(local[2]):
                # 清理可能残留的 _1/_2
                local[0].unlink(missing_ok=True)
                local[1].unlink(missing_ok=True)
                return True
        return False
    else:  # SINGLE
        for i in (2, 0, 1):
            if aria2c_download_single(urls[i], dest):
                if verify_fastq_gzip(local[i]):
                    # 清理不需要的文件
                    for j in range(3):
                        if j != i:
                            local[j].unlink(missing_ok=True)
                    return True
                else:
                    logger.warning(f"[aria2c] {local[i].name} gzip 校验失败，删除后尝试下一个候选")
                    local[i].unlink(missing_ok=True)
        return False


# ============================================================
# curl 下载 (ENA HTTPS，支持 socks 等环境代理)
# ============================================================

def _socks_proxy_in_env() -> Optional[str]:
    """Return the first socks proxy URL found in proxy environment variables.

    Returns
    -------
    Optional[str]
        The proxy URL (e.g. "socks5://127.0.0.1:1080") or None when no
        socks proxy is configured.
    """
    for var in ("all_proxy", "ALL_PROXY", "https_proxy", "HTTPS_PROXY", "http_proxy", "HTTP_PROXY"):
        val = os.environ.get(var, "").strip()
        if val.lower().startswith(("socks4", "socks5")):
            return val
    return None


def _check_size(dest_file: Path, expected_total: Optional[int]) -> bool:
    """Check the downloaded file size against the expected total bytes.

    Parameters
    ----------
    dest_file : Path
        Downloaded file.
    expected_total : Optional[int]
        Expected size in bytes from the HTTP HEAD probe; None skips the check.

    Returns
    -------
    bool
        True when the size matches (or no expectation is available).
    """
    if expected_total is None:
        return True
    actual = dest_file.stat().st_size if dest_file.exists() else 0
    if actual != expected_total:
        logger.warning(f"[curl] {dest_file.name} 大小不符: 实际 {actual} != 预期 {expected_total} 字节")
        return False
    return True


def _probe_remote(url: str, connect_timeout: int = 60) -> Tuple[Optional[int], bool]:
    """Probe the remote file with an HTTP HEAD request.

    Parameters
    ----------
    url : str
        File URL to probe.
    connect_timeout : int, default=60
        curl --connect-timeout in seconds.

    Returns
    -------
    Tuple[Optional[int], bool]
        (content_length, accept_ranges): content length in bytes (None if
        unknown) and whether the server supports HTTP range requests.
    """
    result = subprocess.run(
        ["curl", "-sIL", "--connect-timeout", str(connect_timeout), url],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return None, False
    length: Optional[int] = None
    accept_ranges = False
    for line in result.stdout.splitlines():
        low = line.lower()
        if low.startswith("content-length:"):
            try:
                length = int(low.split(":", 1)[1].strip())
            except ValueError:
                pass
        elif low.startswith("accept-ranges:") and "bytes" in low:
            accept_ranges = True
    return length, accept_ranges


def _curl_download_range(
    url: str,
    dest_file: Path,
    byte_range: Optional[Tuple[int, int]],
    connect_timeout: int = 30,
    retries: int = 5,
) -> bool:
    """Download one file (or one byte range) with curl into dest_file.

    Parameters
    ----------
    url : str
        Complete download URL.
    dest_file : Path
        Output file path.
    byte_range : Optional[Tuple[int, int]]
        Inclusive (start, end) byte range, or None for the whole file.
    connect_timeout : int, default=30
        curl --connect-timeout in seconds.
    retries : int, default=5
        curl --retry: number of retries on transient errors.

    Returns
    -------
    bool
        True when curl exits with code 0.
    """
    cmd = [
        "curl",
        "-L", "--fail", "-sS",
        "--connect-timeout", str(connect_timeout),
        "--retry", str(retries),
        "--retry-delay", "5",
    ]
    if byte_range is not None:
        cmd += ["-r", f"{byte_range[0]}-{byte_range[1]}"]
    else:
        cmd += ["-C", "-"]
    cmd += ["-o", str(dest_file), url]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        for line in (result.stderr or "").strip().splitlines()[-5:]:
            logger.warning(f"[curl] {line}")
        logger.warning(f"[curl] {url} 退出码 {result.returncode}")
        return False
    return True


def curl_download_single(
    url: str,
    dest_dir: Path,
    connections: int = 8,
    connect_timeout: int = 30,
    retries: int = 5,
    min_segment_size: int = 64 * 1048576,
) -> bool:
    """使用 curl 分段并发下载单个文件，走环境代理(含 socks5)

    大文件按 HTTP Range 拆成多段并行下载（网络按连接限速时显著提速），
    小文件或服务器不支持 Range 时退化为单连接下载。

    Parameters
    ----------
    url : str
        Complete download URL, logged verbatim before downloading.
    dest_dir : Path
        Destination directory.
    connections : int, default=8
        Max parallel segments (also the max parallel connections).
    connect_timeout : int, default=30
        curl --connect-timeout in seconds.
    retries : int, default=5
        curl --retry: number of retries on transient errors.
    min_segment_size : int, default=64MiB
        Minimum bytes per segment; smaller segments are merged.

    Returns
    -------
    bool
        True when the file is fully downloaded.

    Notes
    -----
    Unlike aria2c, curl natively honors socks proxy environment variables,
    which is required on networks where direct connections are intercepted
    by a TLS MITM gateway.
    """
    filename = url.rsplit("/", 1)[-1]
    dest_file = dest_dir / filename
    total, ranges_ok = _probe_remote(url, connect_timeout)

    n_segments = 1
    if ranges_ok and total:
        n_segments = min(connections, max(1, total // min_segment_size))
    if n_segments <= 1:
        logger.info(f"[curl] 开始下载 {url} -> {dest_file} (单连接)")
        # 先写 .part 临时文件，完整后才改名为最终文件，保证最终名=完整文件
        tmp_part = dest_dir / f"{filename}.singlepart"
        if not _curl_download_range(url, tmp_part, None, connect_timeout, retries):
            return False
        if not _check_size(tmp_part, total):
            tmp_part.unlink(missing_ok=True)
            return False
        tmp_part.replace(dest_file)
        return True
    assert total is not None

    seg_size = (total + n_segments - 1) // n_segments
    segments: List[Tuple[Path, int, int]] = []
    for i in range(n_segments):
        start = i * seg_size
        end = min(total - 1, (i + 1) * seg_size - 1)
        segments.append((dest_dir / f"{filename}.part{i}", start, end))

    logger.info(
        f"[curl] 开始下载 {url} -> {dest_file} "
        f"(分 {n_segments} 段并发, 共 {total / 1048576:.1f} MiB)"
    )

    # 断点续传：已完成的分段（大小正确）直接跳过
    pending = []
    for part, start, end in segments:
        expected = end - start + 1
        if part.exists() and part.stat().st_size == expected:
            logger.info(f"[curl] {part.name} 已完成 ({expected / 1048576:.1f} MiB)，跳过")
            continue
        part.unlink(missing_ok=True)
        pending.append((part, start, end))

    procs = []
    for part, start, end in pending:
        cmd = [
            "curl",
            "-L", "--fail", "-sS",
            "--connect-timeout", str(connect_timeout),
            "--retry", str(retries),
            "--retry-delay", "5",
            "-r", f"{start}-{end}",
            "-o", str(part),
            url,
        ]
        procs.append((part, subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)))

    # 轮询分段进度
    last_log_time = 0.0
    while any(p.poll() is None for _, p in procs):
        time.sleep(5)
        now = time.time()
        if now - last_log_time >= 5.0:
            done = sum(part.stat().st_size for part, _, _ in segments if part.exists())
            logger.info(f"[curl] {filename}  分段进度 {done / 1048576:.1f}/{total / 1048576:.1f} MiB")
            last_log_time = now

    failed_parts = []
    for part, p in procs:
        if p.returncode != 0:
            err = (p.stderr.read() or "").strip() if p.stderr is not None else ""
            for line in err.splitlines()[-3:]:
                logger.warning(f"[curl] {line}")
            logger.warning(f"[curl] {part.name} 下载失败 (退出码 {p.returncode})")
            failed_parts.append(part)
    if failed_parts:
        return False

    # 分段合并
    tmp_file = dest_dir / f"{filename}.merging"
    with tmp_file.open("wb") as out:
        for part, _, _ in segments:
            with part.open("rb") as src:
                shutil.copyfileobj(src, out)
    tmp_file.replace(dest_file)
    for part, _, _ in segments:
        part.unlink(missing_ok=True)
    return _check_size(dest_file, total)


def curl_download_single_srr(srr_id: str, library_type: str, dest: Path) -> bool:
    """从 ENA HTTPS 使用 curl 下载单个 SRR 的 fastq.gz

    Parameters
    ----------
    srr_id : str
        SRA run accession.
    library_type : str
        "PAIRED" or "SINGLE".
    dest : Path
        Destination directory.

    Returns
    -------
    bool
        True when a valid fastq.gz file is downloaded and verified.
    """
    urls = build_ena_http_urls(srr_id)
    local = [dest / Path(urllib.parse.urlparse(u).path).name for u in urls]
    dest.mkdir(parents=True, exist_ok=True)
    logger.info(f"[curl] {srr_id} 候选 URL:\n    " + "\n    ".join(urls))
    md5_map = fetch_ena_fastq_md5(srr_id)

    def _verify(path: Path) -> bool:
        """gzip 完整性 + 字节级 md5 双重校验"""
        return verify_fastq_gzip(path) and verify_fastq_md5(path, md5_map.get(path.name))

    if library_type == "PAIRED":
        # 优先尝试 paired (_1 + _2)
        ok1 = curl_download_single(urls[0], dest)
        if ok1 and _verify(local[0]):
            ok2 = curl_download_single(urls[1], dest)
            if ok2 and _verify(local[1]):
                return True
            logger.warning(f"[curl] {srr_id}_2 下载失败，尝试 single-end fallback")
            local[1].unlink(missing_ok=True)
        else:
            logger.warning(f"[curl] {srr_id}_1 下载失败，尝试 single-end fallback")
            local[0].unlink(missing_ok=True)

        # fallback: single-end (.fastq.gz)
        if curl_download_single(urls[2], dest):
            if _verify(local[2]):
                # 清理可能残留的 _1/_2
                local[0].unlink(missing_ok=True)
                local[1].unlink(missing_ok=True)
                return True
        return False
    else:  # SINGLE
        for i in (2, 0, 1):
            if curl_download_single(urls[i], dest):
                if _verify(local[i]):
                    # 清理不需要的文件
                    for j in range(3):
                        if j != i:
                            local[j].unlink(missing_ok=True)
                    return True
                else:
                    logger.warning(f"[curl] {local[i].name} 校验失败，删除后尝试下一个候选")
                    local[i].unlink(missing_ok=True)
        return False


# ============================================================
# SRA Toolkit 下载 (prefetch + fasterq-dump)
# ============================================================

def sra_download_single_srr(srr_id: str, library_type: str, dest: Path) -> bool:
    """通过 NCBI SRA Toolkit (prefetch + fasterq-dump) 下载并转换单个 SRR

    Parameters
    ----------
    srr_id : str
        SRA run accession.
    library_type : str
        "PAIRED" or "SINGLE".
    dest : Path
        Destination directory.

    Returns
    -------
    bool
        True when fastq.gz files are generated and pass gzip validation.
    """
    dest.mkdir(parents=True, exist_ok=True)

    # --- 1. prefetch ---
    sra_dir = dest / srr_id
    sra_file = sra_dir / f"{srr_id}.sra"
    # prefetch 2.9.6 拒绝 --max-size 0（"Maximum requested file size is zero"），
    # 其默认上限是 20G，这里用 1T 作为事实上的"不限制"
    prefetch_cmd = ["prefetch", srr_id, "-O", str(dest), "--max-size", "1T"]
    logger.info(f"[prefetch] 执行: {shlex.join(prefetch_cmd)}")
    result = subprocess.run(prefetch_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        _log_subprocess_error("prefetch", result)
        shutil.rmtree(sra_dir, ignore_errors=True)
        return False

    if not sra_file.exists():
        alt = dest / f"{srr_id}.sra"
        if alt.exists():
            sra_file = alt
        else:
            logger.warning(f"[prefetch] {srr_id} 完成但 .sra 文件不存在: {sra_file}")
            return False

    # --- 2. fasterq-dump ---
    fq_cmd = [
        "fasterq-dump", str(sra_file),
        "--outdir", str(dest),
        "--split-files",
        "--threads", "4",
    ]
    logger.info(f"[fasterq-dump] 执行: {shlex.join(fq_cmd)}")
    result = subprocess.run(fq_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        _log_subprocess_error("fasterq-dump", result)
        sra_file.unlink(missing_ok=True)
        shutil.rmtree(sra_dir, ignore_errors=True)
        return False

    # --- 3. gzip ---
    generated = list(dest.glob(f"{srr_id}*.fastq"))
    if not generated:
        logger.warning(f"[fasterq-dump] {srr_id} 未生成 .fastq 文件")
        sra_file.unlink(missing_ok=True)
        shutil.rmtree(sra_dir, ignore_errors=True)
        return False

    for fq in generated:
        gz_cmd = ["gzip", "-f", str(fq)]
        logger.info(f"[gzip] 执行: {shlex.join(gz_cmd)}")
        result = subprocess.run(gz_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            _log_subprocess_error("gzip", result)
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
            logger.warning(f"[gzip] {fq_gz.name} 校验失败")
            ok = False
    return ok


# ============================================================
# 自旋重试逻辑
# ============================================================

def spin_until_success(try_func, desc: str, sleep_base: int, sleep_max: int, max_retries: int = 10) -> bool:
    """Spin with exponential backoff until try_func succeeds or retries run out.

    Parameters
    ----------
    try_func : Callable[[], bool]
        Zero-arg callable performing one download attempt.
    desc : str
        Description used in log messages.
    sleep_base : int
        Initial retry interval in seconds.
    sleep_max : int
        Upper bound of the retry interval in seconds.
    max_retries : int, default=10
        Maximum number of attempts.

    Returns
    -------
    bool
        True on success.

    Notes
    -----
    A :class:`NetworkHijackedError` from try_func aborts immediately (no
    retry) and sets the module-level abort event, because retrying cannot
    help until the network is authenticated.
    """
    attempt = 1
    sleep_time = sleep_base

    while attempt <= max_retries:
        logger.info(f"[Attempt {attempt}/{max_retries}] {desc}")

        try:
            ok = try_func()
        except NetworkHijackedError as e:
            logger.error(f"[ABORT] {desc}: {e}")
            logger.error("[ABORT] 网络被认证门户劫持，重试无效，中止剩余下载；请先完成网络认证(如校园网/Dr.COM 登录)后重新运行")
            _ABORT_EVENT.set()
            return False

        if ok:
            logger.info(f"[SUCCESS] {desc}")
            return True

        logger.warning(f"[FAIL] {desc}，{sleep_time}s 后重试")
        time.sleep(sleep_time)
        sleep_time = min(sleep_time * 2, sleep_max)
        attempt += 1

    logger.error(f"[GIVE UP] {desc}，已重试 {max_retries} 次")
    return False


# ============================================================
# 单 SRR 下载（主调度）
# ============================================================

def _fastq_exists(dest: Path, srr_id: str, library_type: str) -> bool:
    """检查目标 fastq 文件是否已存在且 gzip 完整"""
    if library_type == "PAIRED":
        f1 = dest / f"{srr_id}_1.fastq.gz"
        f2 = dest / f"{srr_id}_2.fastq.gz"
        if f1.exists() and f2.exists() and gzip_test(f1) and gzip_test(f2):
            return True
    fq = dest / f"{srr_id}.fastq.gz"
    if fq.exists() and gzip_test(fq):
        return True
    return False


def download_spin(
    srr_id: str,
    library_type: str,
    dest: Path,
    method: str,
    sleep_base: int,
    sleep_max: int,
    max_retries: int = 10,
    key: Optional[Path] = None,
    globus_tc: Optional["globus_sdk.TransferClient"] = None,
    globus_src_ep: Optional[str] = None,
    globus_dest_ep: Optional[str] = None,
) -> bool:
    if _ABORT_EVENT.is_set():
        logger.warning(f"[SKIP] {srr_id} {library_type}：已检测到网络劫持，跳过剩余下载")
        return False

    # 缓存：已存在且完整则跳过
    if _fastq_exists(dest, srr_id, library_type):
        logger.info(f"[SKIP] {srr_id} {library_type} 已存在，跳过")
        return True

    logger.info(f"{srr_id} {library_type} 开始下载 (method={method})")

    # aria2c 不支持 socks 代理(unrecognized proxy format)，直连又会被网关 TLS 劫持，
    # 检测到 socks 代理时自动切换 curl 引擎（curl 原生支持 socks 环境变量）
    use_curl = method == "curl"
    if method == "aria2c":
        socks = _socks_proxy_in_env()
        if socks:
            logger.warning(f"[aria2c] 检测到 socks 代理 {socks}，aria2c 无法使用，自动切换 curl 引擎")
            use_curl = True

    if method == "ascp":
        if not key:
            raise ValueError("ascp 方法需要 --key 参数")
        try_func = lambda: ena_download_single_srr(srr_id, library_type, dest, key)
    elif use_curl:
        try_func = lambda: curl_download_single_srr(srr_id, library_type, dest)
    elif method == "aria2c":
        try_func = lambda: aria2c_download_single_srr(srr_id, library_type, dest)
    elif method == "globus":
        if not globus_tc or not globus_dest_ep:
            raise ValueError("globus 方法缺失 transfer_client 或 destination endpoint")
        try_func = lambda: globus_download_single_srr(globus_tc, globus_src_ep, globus_dest_ep, srr_id, library_type, dest)
    else:  # sra
        try_func = lambda: sra_download_single_srr(srr_id, library_type, dest)

    return spin_until_success(try_func, f"{srr_id} {library_type}", sleep_base, sleep_max, max_retries)


# ============================================================
# SRR 解析
# ============================================================

def load_tasks(args) -> List[Tuple[str, str]]:
    """Parse input mode (meta table / srr list / single srr) into (srr, layout) tasks.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.

    Returns
    -------
    List[Tuple[str, str]]
        (srr_id, library_layout) pairs.

    Raises
    ------
    ValueError
        If the meta file cannot be read, the SRR column is missing, or
        library layout values are invalid.
    """
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

        if srr_col not in df.columns:
            raise ValueError(
                f"Missing required column '{srr_col}' in meta file. "
                f"Available columns: {list(df.columns)}"
            )

        # 只要求 srr 非空，library_layout 缺失时用 default_layout 兜底
        df = df.dropna(subset=[srr_col])

        if lib_col in df.columns:
            df[lib_col] = df[lib_col].fillna(args.default_layout).str.upper()
            invalid = df[~df[lib_col].isin({"PAIRED", "SINGLE"})]
            if not invalid.empty:
                raise ValueError(f"Invalid library layout values found:\n{invalid}")
        else:
            logger.warning(f"列 '{lib_col}' 不存在，全部使用 {args.default_layout}")
            df[lib_col] = args.default_layout

        tasks = list(df[[srr_col, lib_col]].itertuples(index=False, name=None))
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
    p.add_argument("-m", "--method", choices=["ascp", "aria2c", "curl", "sra", "globus"], default="sra",
                   help="Download method: ascp (ENA fasp), aria2c (ENA HTTP multi-thread), curl (ENA HTTPS, supports socks proxy), sra (NCBI prefetch, default), globus")

    # ---------- Output / logging ----------
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("-l", "--log", type=Path, default=None,
                   help="日志文件路径 (可选；缺省时日志输出到标准输出)")

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
    p.add_argument("--max-retries", type=int, default=10, help="单个 SRR 最大重试次数 (默认: 10)")
    p.add_argument("--default-layout", default="SINGLE", choices=["PAIRED", "SINGLE"],
                   help="library_layout 缺失时的默认值 (默认: SINGLE)")

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
    # 初始化全局 logger（LogUtil.setup_logger 按名字配置，与模块级 logger 是同一对象）
    # 未提供 --log 时日志仅输出到标准输出
    log_file = str(args.log) if args.log else None
    setup_logger("SRA_DOWNLOAD", log_file=log_file,
                 stream=None if log_file else sys.stdout)

    tasks = load_tasks(args)
    logger.info(f"共 {len(tasks)} 个 SRR，method={args.method}，jobs={args.jobs}")

    # 若使用 Globus 则提前初始化客户端避免重复登录
    globus_tc = None
    if args.method == "globus":
        globus_tc = get_globus_client(args.globus_client_id, args.globus_token)

    def download_fn(srr, lib):
        return download_spin(
            srr, lib, args.outdir, args.method,
            args.sleep_base, args.sleep_max,
            max_retries=args.max_retries,
            key=args.key,
            globus_tc=globus_tc,
            globus_src_ep=args.globus_source_ep,
            globus_dest_ep=args.globus_dest_ep
        )

    failed: List[str] = []
    if args.jobs == 1:
        for srr, lib in tasks:
            if not download_fn(srr, lib):
                failed.append(srr)
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            future_map = {
                ex.submit(download_fn, srr, lib): srr
                for srr, lib in tasks
            }
            for future in as_completed(future_map):
                srr = future_map[future]
                try:
                    if not future.result():
                        failed.append(srr)
                except Exception as e:
                    logger.error(f"[ERROR] {srr}: {e}")
                    failed.append(srr)

    if failed:
        if _ABORT_EVENT.is_set():
            logger.error("===== 运行中止：下载内容被网络认证门户(如 Dr.COM)劫持，先完成网络认证再重跑 =====")
        logger.error(f"===== {len(failed)}/{len(tasks)} 个 SRR 下载失败 =====")
        for srr in sorted(failed):
            logger.error(f"  FAILED: {srr}")
        # 写入失败列表文件
        failed_file = args.outdir / "download_failed.txt"
        args.outdir.mkdir(parents=True, exist_ok=True)
        with open(failed_file, "w") as f:
            for srr in sorted(failed):
                f.write(f"{srr}\n")
        logger.error(f"失败列表已写入: {failed_file}")
        sys.exit(1)
    else:
        logger.info(f"===== 全部 {len(tasks)} 个 SRR 下载成功 =====")


if __name__ == "__main__":
    main()