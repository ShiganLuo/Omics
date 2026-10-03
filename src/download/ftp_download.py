#!/usr/bin/env python3
"""Parameterized FTP batch downloader for omics raw data.

Collects FTP URLs from one of three input modes, then downloads them in
parallel with retries, atomic writes and skip-existing support.

Input modes
-----------
urls      Plain text file, one record per line; the URL is taken from a
          configurable column (whitespace delimited).
table     Delimited table (e.g. CNGB run.tsv); FTP columns are detected
          from the cell content (ftp:// URLs and HTTPS FTP-archive links
          such as https://ftp.pride.ebi.ac.uk/...) by default, or forced
          by name with --url-pattern.
metadata  CNGBdb experiment metadata + FTP reference list; records are
          filtered by a condition column and matched against run accessions.

Usage:
    # 1. plain URL list
    python ftp_download.py urls -f links.txt -o fastq/ -l log/ftp_download.log

    # 2. table: FTP columns auto-detected from cell content
    python ftp_download.py table -t run.tsv -o fastq/ -l log/ftp_download.log

    # or force columns by name
    python ftp_download.py table -t run.tsv -p DownLoad -o fastq/ -l log/ftp_download.log

    # 3. metadata-driven (condition filter + run accession match)
    python ftp_download.py metadata -m metadata_experiment.tsv \
        -R data_download_links_ftp.txt \
        -C experiment_title \
        -V "B2 CoD2 Human p65KO" \
        -o fastq/ -l log/ftp_download.log

Common options: -j/--jobs, -r/--retries, -T/--timeout, -w/--overwrite,
-n/--dry-run, -u/--username, -P/--password. Logs go to stdout unless a
log file is given via -l/--log.

Logging follows src/common/util/LogUtil.py:
    "<time> | <level> | FTPDownloader | <message>"
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys
import time
import ftplib
import http.client
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from urllib.parse import urlparse

import pandas as pd

_SRC_DIR = Path(__file__).resolve().parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger
from src.common.util.SepUtil import detect_delimiter
from src.common.util.MatchUtil import run_accession_match

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_TIMEOUT = 60        # FTP socket timeout in seconds
DEFAULT_JOBS = 4            # parallel download workers
DEFAULT_RETRIES = 3         # per-file retry attempts
BACKOFF_BASE = 10           # first retry sleep in seconds (doubles each try)
BACKOFF_MAX = 300           # retry sleep upper bound in seconds
PART_SUFFIX = ".part"       # temporary suffix for in-progress downloads

def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    """Configure the module logger.

    Logs go to stdout; when *log_file* is given they are written to that file
    as well (parent directories are created as needed).
    """
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("FTPDownloader", log_file=str(log_file))
    else:
        lg = setup_logger("FTPDownloader")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:  # console only, not FileHandler
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


# ---------------------------------------------------------------------------
# URL helpers
# ---------------------------------------------------------------------------
def parse_ftp_url(url: str) -> Tuple[str, str]:
    """Split an ftp:// URL into (server, absolute remote path).

    Parameters
    ----------
    url : str
        Full URL, e.g. "ftp://ftp.example.com/path/to/file.fastq.gz".

    Returns
    -------
    (server, remote_path) : Tuple[str, str]
        Host name and absolute remote path ("//"-free, starts with "/").

    Raises
    ------
    ValueError
        If the URL is not an ftp:// URL or has no host part.
    """
    if not url.startswith("ftp://"):
        raise ValueError(f"Not an ftp:// URL: {url}")
    parsed = urlparse(url)
    if not parsed.hostname:
        raise ValueError(f"Missing host in URL: {url}")
    remote_path = parsed.path or "/"
    return parsed.hostname, remote_path


def strip_server_prefix(remote_file_path: str, server: str) -> str:
    """Normalize a remote path that may still carry an "ftp://server" prefix.

    Accepts both "ftp://server/a/b.fastq" and "/a/b.fastq" and returns the
    bare absolute remote path.
    """
    if remote_file_path.startswith("ftp://"):
        parsed = urlparse(remote_file_path)
        remote_file_path = parsed.path or "/"
    return remote_file_path if remote_file_path.startswith("/") else "/" + remote_file_path


# ---------------------------------------------------------------------------
# FTP connection
# ---------------------------------------------------------------------------
def open_ftp(
    server: str,
    username: Optional[str] = None,
    password: Optional[str] = None,
    timeout: int = DEFAULT_TIMEOUT,
) -> ftplib.FTP:
    """Open and authenticate an FTP connection.

    Uses credential login when both username and password are given,
    anonymous login otherwise. The caller is responsible for closing.
    """
    ftp = ftplib.FTP(server, timeout=timeout)
    if username and password:
        ftp.login(user=username, passwd=password)
    else:
        ftp.login()
    return ftp


# ---------------------------------------------------------------------------
# Single-file download
# ---------------------------------------------------------------------------
def download_file_from_ftp(
    server: str,
    remote_file_path: str,
    local_file_path: str,
    username: Optional[str] = None,
    password: Optional[str] = None,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    overwrite: bool = False,
) -> bool:
    """Download a single file over FTP with retries and atomic rename.

    The file is first written to "<local_file_path>.part" and renamed on
    success, so an interrupted transfer never leaves a truncated file at the
    final path. An existing final file is skipped unless overwrite=True.

    Parameters
    ----------
    server : str
        FTP server host.
    remote_file_path : str
        Remote path, optionally prefixed with "ftp://server".
    local_file_path : str
        Destination file path; parent directories are created as needed.
    username, password : Optional[str]
        Credentials; anonymous login when omitted.
    timeout : int
        Socket timeout in seconds.
    retries : int
        Maximum download attempts (exponential backoff between attempts).
    overwrite : bool
        Re-download even if the destination file already exists.

    Returns
    -------
    bool
        True on success, False if all attempts failed.
    """
    remote_file_path = strip_server_prefix(remote_file_path, server)
    local_path = Path(local_file_path)

    if local_path.exists() and not overwrite:
        logger.info(f"[SKIP] already exists: {local_path}")
        return True

    local_path.parent.mkdir(parents=True, exist_ok=True)
    part_path = local_path.with_name(local_path.name + PART_SUFFIX)

    sleep_time = BACKOFF_BASE
    for attempt in range(1, retries + 1):
        ftp: Optional[ftplib.FTP] = None
        try:
            logger.info(
                f"[Attempt {attempt}/{retries}] RETR {remote_file_path} -> {local_path}"
            )
            ftp = open_ftp(server, username, password, timeout)
            with open(part_path, "wb") as fh:
                ftp.retrbinary(f"RETR {remote_file_path}", fh.write)
            ftp.quit()
            ftp = None
            os.replace(part_path, local_path)
            logger.info(f"[OK] {local_path}")
            return True
        except ftplib.all_errors as e:
            logger.warning(f"[FAIL] {remote_file_path}: {e}")
        except OSError as e:
            logger.warning(f"[FAIL] local I/O for {local_path}: {e}")
        finally:
            if ftp is not None:
                try:
                    ftp.close()
                except Exception:
                    pass
        if attempt < retries:
            logger.info(f"retry in {sleep_time}s")
            time.sleep(sleep_time)
            sleep_time = min(sleep_time * 2, BACKOFF_MAX)

    part_path.unlink(missing_ok=True)
    logger.error(f"[GIVE UP] {remote_file_path} after {retries} attempts")
    return False


def download_file_from_http(
    url: str,
    local_file_path: str,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    overwrite: bool = False,
) -> bool:
    """Download a single file over HTTP(S) with retries and atomic rename.

    Same semantics as download_file_from_ftp: the file is first written to
    "<local_file_path>.part" and renamed on success; an existing final file
    is skipped unless overwrite=True. Used for FTP archives served over
    HTTPS (e.g. https://ftp.pride.ebi.ac.uk/...).

    Parameters
    ----------
    url : str
        Full http:// or https:// URL.
    local_file_path : str
        Destination file path; parent directories are created as needed.
    timeout : int
        Socket timeout in seconds.
    retries : int
        Maximum download attempts (exponential backoff between attempts).
    overwrite : bool
        Re-download even if the destination file already exists.

    Returns
    -------
    bool
        True on success, False if all attempts failed.
    """
    local_path = Path(local_file_path)

    if local_path.exists() and not overwrite:
        logger.info(f"[SKIP] already exists: {local_path}")
        return True

    local_path.parent.mkdir(parents=True, exist_ok=True)
    part_path = local_path.with_name(local_path.name + PART_SUFFIX)

    sleep_time = BACKOFF_BASE
    for attempt in range(1, retries + 1):
        try:
            logger.info(f"[Attempt {attempt}/{retries}] GET {url} -> {local_path}")
            req = urllib.request.Request(url, headers={"User-Agent": "ftp_download"})
            with urllib.request.urlopen(req, timeout=timeout) as resp, \
                    open(part_path, "wb") as fh:
                shutil.copyfileobj(resp, fh)
            os.replace(part_path, local_path)
            logger.info(f"[OK] {local_path}")
            return True
        except (urllib.error.URLError, OSError, http.client.HTTPException) as e:
            logger.warning(f"[FAIL] {url}: {e}")
        if attempt < retries:
            logger.info(f"retry in {sleep_time}s")
            time.sleep(sleep_time)
            sleep_time = min(sleep_time * 2, BACKOFF_MAX)

    part_path.unlink(missing_ok=True)
    logger.error(f"[GIVE UP] {url} after {retries} attempts")
    return False


# ---------------------------------------------------------------------------
# Directory download
# ---------------------------------------------------------------------------
def download_files_from_ftp(
    server: str,
    remote_dir: str = "/",
    local_dir: str = "./downloads",
    username: Optional[str] = None,
    password: Optional[str] = None,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    overwrite: bool = False,
) -> List[str]:
    """Download every file of a remote directory (non-recursive).

    Parameters
    ----------
    server : str
        FTP server host.
    remote_dir : str
        Remote directory, optionally prefixed with "ftp://server".
    local_dir : str
        Local directory to save the downloaded files.
    username, password : Optional[str]
        Credentials; anonymous login when omitted.
    timeout : int
        Socket timeout in seconds.
    retries : int
        Per-file maximum download attempts.
    overwrite : bool
        Re-download even if a destination file already exists.

    Returns
    -------
    List[str]
        Remote paths that failed to download (empty when all succeeded).
    """
    remote_dir = strip_server_prefix(remote_dir, server)
    Path(local_dir).mkdir(parents=True, exist_ok=True)

    logger.info(f"Connecting to FTP server: {server}, remote dir: {remote_dir}")
    ftp = open_ftp(server, username, password, timeout)
    try:
        ftp.cwd(remote_dir)
        # mlsd gives file type info; fall back to a bare name list
        try:
            entries = list(ftp.mlsd())
            names = [name for name, facts in entries if facts.get("type") == "file"]
            skipped = [name for name, facts in entries if facts.get("type") != "file"]
        except ftplib.all_errors:
            names = ftp.nlst()
            skipped = []
        logger.info(f"Files found in {remote_dir}: {names}")
        if skipped:
            logger.info(f"Skipped non-file entries: {skipped}")
    finally:
        try:
            ftp.quit()
        except ftplib.all_errors:
            ftp.close()

    failed: List[str] = []
    for name in names:
        ok = download_file_from_ftp(
            server,
            f"{remote_dir.rstrip('/')}/{name}",
            os.path.join(local_dir, name),
            username=username,
            password=password,
            timeout=timeout,
            retries=retries,
            overwrite=overwrite,
        )
        if not ok:
            failed.append(name)
    return failed


# ---------------------------------------------------------------------------
# URL collection: plain list file
# ---------------------------------------------------------------------------
def collect_urls_from_file(
    url_file: Path,
    url_column: int = 1,
) -> List[str]:
    """Read whitespace-delimited lines and take the URL from one column.

    Blank lines and lines starting with "#" are ignored.

    Parameters
    ----------
    url_file : Path
        Input text file, one record per line.
    url_column : int
        0-based column index holding the URL (default: 1, i.e. second column).

    Returns
    -------
    List[str]
        URLs in file order.
    """
    urls: List[str] = []
    with open(url_file, "r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) <= url_column:
                logger.warning(f"{url_file}:{lineno} has no column {url_column}, skipped")
                continue
            urls.append(fields[url_column])
    logger.info(f"Collected {len(urls)} URLs from {url_file} (column {url_column})")
    return urls


# ---------------------------------------------------------------------------
# URL collection: table columns
# ---------------------------------------------------------------------------
def _is_ftp_value(value: str) -> bool:
    """True when a stripped cell value looks like a downloadable archive URL.

    Recognizes "ftp://" URLs and HTTP(S) URLs of FTP archives served over
    HTTPS (e.g. "https://ftp.pride.ebi.ac.uk/..." for PRIDE, whose host
    starts with "ftp.").
    """
    v = value.strip()
    low = v.lower()
    if low.startswith("ftp://"):
        return True
    if low.startswith(("http://", "https://")):
        host = urlparse(v).hostname or ""
        return host.lower().startswith("ftp.")
    return False


def collect_urls_from_table(
    table_file: Path,
    url_columns: Optional[Sequence[str]] = None,
) -> List[str]:
    """Collect FTP URLs from a delimited table.

    Column selection:
    - *url_columns* given (forced): every column whose name contains one of
      the given names is used whatever its content; if no column matches, a
      ValueError is raised.
    - otherwise (default): columns are detected from the cell content -- every
      column holding at least one download URL ("ftp://" or an HTTPS
      FTP-archive link such as "https://ftp.pride.ebi.ac.uk/...") is treated
      as an FTP column.

    Parameters
    ----------
    table_file : Path
        Delimited table with a header row (delimiter auto-detected).
    url_columns : Optional[Sequence[str]]
        Column-name substrings forcing the selection, e.g. ("DownLoad",).

    Returns
    -------
    List[str]
        Collected URLs. Content mode keeps only download-URL values; forced
        mode keeps all non-empty values of the named columns.

    Raises
    ------
    ValueError
        If no column matches the given names (forced mode), or no download
        URL is found in the table at all (content mode).
    """
    sep = detect_delimiter(str(table_file))
    read_sep = r"\s+" if sep == "whitespace" else sep
    df = pd.read_csv(table_file, sep=read_sep)

    urls: List[str] = []
    if url_columns:
        url_cols = [col for col in df.columns if any(p in col for p in url_columns)]
        if not url_cols:
            raise ValueError(
                f"No column of {table_file} matches {list(url_columns)}. "
                f"Available columns: {list(df.columns)}"
            )
        for col in url_cols:
            urls.extend(df[col].dropna().astype(str).tolist())
        logger.info(
            f"Collected {len(urls)} URLs from forced columns {url_cols} of {table_file}"
        )
    else:
        url_cols = []
        for col in df.columns:
            vals = [v.strip() for v in df[col].dropna().astype(str) if _is_ftp_value(v)]
            if vals:
                url_cols.append(col)
                urls.extend(vals)
        if not urls:
            raise ValueError(f"No download URLs (ftp:// or https://ftp.*) found in {table_file}.")
        logger.info(
            f"Collected {len(urls)} URLs from content-detected FTP columns "
            f"{url_cols} of {table_file}"
        )
    return urls


# ---------------------------------------------------------------------------
# URL collection: metadata + run accession match
# ---------------------------------------------------------------------------
def collect_urls_from_metadata(
    meta_file: Path,
    ftp_ref_file: Path,
    output_dir: Path,
    condition_col: str,
    condition_values: Sequence[str],
    run_accession_col: str = "run_accession",
    url_column: int = 1,
    exact_match: bool = False,
) -> List[str]:
    """Filter experiment metadata by condition and match run accessions to URLs.

    Workflow (CNGBdb metadata):
    1. keep metadata rows whose `condition_col` matches any condition value;
    2. take the run accessions of those rows;
    3. match them against the URLs of the FTP reference file;
    4. write "matched_urls.txt" and "unmatched_accessions.txt" to output_dir.

    Parameters
    ----------
    meta_file : Path
        Experiment metadata table (delimiter auto-detected).
    ftp_ref_file : Path
        Whitespace-delimited reference file holding FTP URLs in `url_column`.
    output_dir : Path
        Directory for the matched/unmatched reports.
    condition_col : str
        Column name to filter on.
    condition_values : Sequence[str]
        Accepted values of `condition_col`; a row matches when it matches any.
    run_accession_col : str
        Column holding the run accessions.
    url_column : int
        0-based column index of the URL inside `ftp_ref_file`.
    exact_match : bool
        Match the whole cell value instead of a substring; substring matching
        is literal (no regex interpretation).

    Returns
    -------
    List[str]
        Matched FTP URLs (may contain duplicates across condition values).

    Raises
    ------
    ValueError
        If a required column is missing from the metadata table.
    """
    sep = detect_delimiter(str(meta_file))
    read_sep = r"\s+" if sep == "whitespace" else sep
    df = pd.read_csv(meta_file, sep=read_sep)

    for col in (condition_col, run_accession_col):
        if col not in df.columns:
            raise ValueError(
                f"Missing required column '{col}' in {meta_file}. "
                f"Available columns: {list(df.columns)}"
            )

    mask = pd.Series(False, index=df.index)
    for value in condition_values:
        if exact_match:
            mask |= df[condition_col].astype(str).eq(value)
        else:
            mask |= df[condition_col].astype(str).str.contains(value, na=False, regex=False)
    df_filtered = df[mask]
    logger.info(
        f"Filtered metadata: {len(df_filtered)} records match {list(condition_values)} "
        f"in column '{condition_col}'"
    )

    ftp_urls = collect_urls_from_file(ftp_ref_file, url_column=url_column)
    run_accessions = df_filtered[run_accession_col].dropna().astype(str).tolist()
    matched_urls, unmatched_accessions = run_accession_match(ftp_urls, run_accessions)

    output_dir.mkdir(parents=True, exist_ok=True)
    unmatched_file = output_dir / "unmatched_accessions.txt"
    with open(unmatched_file, "w", encoding="utf-8") as fh:
        for acc in unmatched_accessions:
            fh.write(f"{acc}\n")
    matched_file = output_dir / "matched_urls.txt"
    with open(matched_file, "w", encoding="utf-8") as fh:
        for url in matched_urls:
            fh.write(f"{url}\n")
    logger.info(
        f"Matched {len(matched_urls)} URLs, {len(unmatched_accessions)} accessions unmatched "
        f"(reports: {matched_file}, {unmatched_file})"
    )
    return matched_urls


# ---------------------------------------------------------------------------
# Batch download
# ---------------------------------------------------------------------------
def plan_downloads(
    urls: Sequence[str],
    output_dir: Path,
) -> List[Tuple[str, str, str, str]]:
    """Resolve URLs to download tasks: (url, server, remote_path, local_path).

    Accepts ftp:// and http(s):// URLs; anything else is dropped with a
    warning. HTTP(S) URLs are transferred over HTTP(S), so FTP archives
    served over HTTPS (e.g. https://ftp.pride.ebi.ac.uk/...) work as-is.
    """
    tasks: List[Tuple[str, str, str, str]] = []
    for url in urls:
        if url.lower().startswith("ftp://"):
            try:
                server, remote_path = parse_ftp_url(url)
            except ValueError as e:
                logger.warning(f"[SKIP] {e}")
                continue
        elif url.lower().startswith(("http://", "https://")):
            parsed = urlparse(url)
            server = parsed.hostname or ""
            remote_path = parsed.path or "/"
            if not server:
                logger.warning(f"[SKIP] Missing host in URL: {url}")
                continue
        else:
            logger.warning(f"[SKIP] Unsupported URL scheme: {url}")
            continue
        local_path = output_dir / Path(remote_path).name
        tasks.append((url, server, remote_path, str(local_path)))
    return tasks


def download_urls(
    urls: Sequence[str],
    output_dir: Path,
    jobs: int = DEFAULT_JOBS,
    username: Optional[str] = None,
    password: Optional[str] = None,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    overwrite: bool = False,
) -> List[str]:
    """Download all ftp:// URLs into output_dir in parallel.

    Returns
    -------
    List[str]
        URLs that failed; also written to "<output_dir>/download_failed.txt".
    """
    tasks = plan_downloads(urls, output_dir)
    logger.info(f"Downloading {len(tasks)} files into {output_dir} (jobs={jobs})")

    def _run(task: Tuple[str, str, str, str]) -> Tuple[str, bool]:
        url, server, remote_path, local_path = task
        logger.info(f"[PLAN] {server} {remote_path} -> {local_path}")
        if url.lower().startswith("ftp://"):
            ok = download_file_from_ftp(
                server, remote_path, local_path,
                username=username, password=password,
                timeout=timeout, retries=retries, overwrite=overwrite,
            )
        else:
            ok = download_file_from_http(
                url, local_path,
                timeout=timeout, retries=retries, overwrite=overwrite,
            )
        return url, ok

    failed: List[str] = []
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
        future_map = {executor.submit(_run, task): task[0] for task in tasks}
        for future in as_completed(future_map):
            url = future_map[future]
            try:
                _, ok = future.result()
                if not ok:
                    failed.append(url)
            except Exception as e:
                logger.error(f"[ERROR] {url}: {e}")
                failed.append(url)

    if failed:
        failed_file = output_dir / "download_failed.txt"
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(failed_file, "w", encoding="utf-8") as fh:
            for url in sorted(failed):
                fh.write(f"{url}\n")
        logger.error(f"{len(failed)}/{len(tasks)} downloads failed; list: {failed_file}")
    else:
        logger.info(f"All {len(tasks)} downloads completed.")
    return failed


# ---------------------------------------------------------------------------
# argparse
# ---------------------------------------------------------------------------
def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="ftp_download.py",
        description="Parameterized FTP batch downloader (URL list / table / metadata input)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest="mode", required=True)

    def add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("-o", "--outdir", type=Path, required=True,
                        help="directory to save downloaded files")
        sp.add_argument("-l", "--log", type=Path, default=None,
                        help="log file path (default: log to stdout)")
        sp.add_argument("-j", "--jobs", type=int, default=DEFAULT_JOBS,
                        help="parallel download workers")
        sp.add_argument("-r", "--retries", type=int, default=DEFAULT_RETRIES,
                        help="max attempts per file")
        sp.add_argument("-T", "--timeout", type=int, default=DEFAULT_TIMEOUT,
                        help="FTP socket timeout in seconds")
        sp.add_argument("-u", "--username", help="FTP username (default: anonymous)")
        sp.add_argument("-P", "--password", help="FTP password (default: anonymous)")
        sp.add_argument("-w", "--overwrite", action="store_true",
                        help="re-download files that already exist locally")
        sp.add_argument("-n", "--dry-run", action="store_true",
                        help="only log planned downloads, do not transfer")

    # ---- urls mode ----
    sp_urls = sub.add_parser("urls", help="download URLs listed in a text file",
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp_urls.add_argument("-f", "--url-file", type=Path, required=True,
                         help="text file, one record per line")
    sp_urls.add_argument("-c", "--url-column", type=int, default=1,
                         help="0-based column index holding the URL")
    add_common(sp_urls)

    # ---- table mode ----
    sp_table = sub.add_parser("table", help="download URLs from matching table columns",
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp_table.add_argument("-t", "--table", type=Path, required=True,
                          help="delimited table with header (delimiter auto-detected)")
    sp_table.add_argument("-p", "--url-pattern", action="append", default=None,
                          help="force use of columns whose name contains this "
                               "string (repeatable; default: detect FTP columns "
                               "by cell content)")
    add_common(sp_table)

    # ---- metadata mode ----
    sp_meta = sub.add_parser(
        "metadata",
        help="filter metadata by condition, match run accessions, download URLs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sp_meta.add_argument("-m", "--meta", type=Path, required=True,
                         help="experiment metadata table (delimiter auto-detected)")
    sp_meta.add_argument("-R", "--ftp-ref", type=Path, required=True,
                         help="reference file holding FTP URLs")
    sp_meta.add_argument("-C", "--condition-col", default="experiment_title",
                         help="column to filter on")
    sp_meta.add_argument("-V", "--condition-value", action="append", required=True,
                         help="accepted value of the condition column (repeatable)")
    sp_meta.add_argument("-A", "--run-accession-col", default="run_accession",
                         help="column holding run accessions")
    sp_meta.add_argument("-c", "--url-column", type=int, default=1,
                         help="0-based column index of the URL in --ftp-ref")
    sp_meta.add_argument("-e", "--exact-match", action="store_true",
                         help="match the whole condition value instead of a substring")
    add_common(sp_meta)

    args = p.parse_args(argv)
    return args


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")

    # ---- collect URLs according to the input mode ----
    if args.mode == "urls":
        urls = collect_urls_from_file(args.url_file, url_column=args.url_column)
    elif args.mode == "table":
        urls = collect_urls_from_table(args.table, args.url_pattern)
    else:  # metadata
        urls = collect_urls_from_metadata(
            meta_file=args.meta,
            ftp_ref_file=args.ftp_ref,
            output_dir=args.outdir,
            condition_col=args.condition_col,
            condition_values=args.condition_value,
            run_accession_col=args.run_accession_col,
            url_column=args.url_column,
            exact_match=args.exact_match,
        )
    logger.info(f"Total FTP URLs to download: {len(urls)}")

    # ---- download ----
    tasks = plan_downloads(urls, args.outdir)
    if args.dry_run:
        for url, server, remote_path, local_path in tasks:
            logger.info(f"[DRY-RUN] {server} {remote_path} -> {local_path}")
        logger.info(f"[DRY-RUN] {len(tasks)} files planned, nothing downloaded")
        return 0

    failed = download_urls(
        urls,
        output_dir=args.outdir,
        jobs=args.jobs,
        username=args.username,
        password=args.password,
        timeout=args.timeout,
        retries=args.retries,
        overwrite=args.overwrite,
    )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
