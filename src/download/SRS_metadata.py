#!/usr/bin/env python3
"""
SRS Metadata Downloader & Extractor
====================================

Fetch sample- and run-level metadata for SRA sample accessions (SRS...)
from the ENA API and consolidate it into TSV tables.

Data Sources
------------
  1. ENA Browser API (sample XML):
       https://www.ebi.ac.uk/ena/browser/api/xml/<SRS>
     -> BioSample accession, sample alias/title, organism, and every
        SAMPLE_ATTRIBUTE (sex, tissue, strain, ...) as extra columns.
  2. ENA Portal API (filereport, run table):
       https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<SRS>&result=read_run
     -> one row per run: study/experiment/run accessions, library info,
        read counts, fastq_ftp + fastq_md5 + fastq_bytes.

Input Sources
-------------
  1. -s/--srs SRS1 SRS2 ...     Explicit accession list
  2. -i/--input file            File with an ID column (--id-column)
                                or one accession per line
  3. --group LABEL=ACC1,...     Accessions named in group specs are
                                collected too (see Group Labels)
  4. stdin (pipe)               Piped accession list

Group Labels
------------
  -g/--group LABEL=ACC1,ACC2,...   Repeatable, e.g. -g experiment=SRS1,SRS2
  --group-file FILE                Table with accession + group columns
                                   (--group-id-col / --group-name-col)
  A `group` column is always written into both output tables; samples
  without an assignment get an empty cell and a warning.

Outputs (in --outdir)
---------------------
  sample_metadata.tsv   one row per SRS: srs, group, biosample, sample_alias,
                        sample_title, organism, tax_id, + sample attributes
  run_metadata.tsv      one row per run: srs, group, run_accession, ...,
                        library_*, read_count, base_count, fastq_*
  cache/<SRS>.xml       raw ENA sample XML (reused on rerun unless -w)
  cache/<SRS>_runs.tsv  raw ENA run table (reused on rerun unless -w)
  srs_failed.txt        failed accessions (written only on failure,
                        exit code 1)

Examples
--------
  # Fetch a batch with group labels
  python SRS_metadata.py \\
    -g experiment=SRS18518057,SRS18518052,SRS18518053 \\
    -g control=SRS18527637,SRS18527634 \\
    -o /data/16T/Data/RNAseq_Ovary_CD38_PRJNA1002222 \\
    --jobs 4 --log /data/16T/Data/RNAseq_Ovary_CD38_PRJNA1002222/srs_metadata.log

  # Accessions from a file, groups from a second table
  python SRS_metadata.py -i srs_list.txt \\
    --group-file groups.tsv --id-column SRS \\
    -o results/

  # Pipe accession list, no groups
  cut -f1 samples.tsv | python SRS_metadata.py -o results/

Dependencies
------------
  requests; sibling module AccessionResolver (src/download)
"""

import argparse
import csv
import io
import logging
import os
import re
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

_SCRIPT_DIR = Path(__file__).resolve().parent
_SRC_DIR = _SCRIPT_DIR.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from src.common.util.LogUtil import setup_logger
from src.common.util.SepUtil import detect_delimiter
from AccessionResolver import AccessionResolver

logger = logging.getLogger(__name__)


ENA_SAMPLE_XML_URL = "https://www.ebi.ac.uk/ena/browser/api/xml/{acc}"
ENA_FILEREPORT_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run&fields={fields}&format=tsv"
)

RUN_FIELDS = [
    "run_accession",
    "study_accession",
    "study_title",
    "sample_accession",
    "secondary_sample_accession",
    "experiment_accession",
    "experiment_title",
    "sample_title",
    "scientific_name",
    "instrument_platform",
    "instrument_model",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]

SAMPLE_BASE_FIELDS = [
    "srs",
    "group",
    "biosample",
    "sample_alias",
    "sample_title",
    "organism",
    "tax_id",
]

SRS_PATTERN = re.compile(r"^SRS\d+$", re.IGNORECASE)

_thread_local = threading.local()



def _get_session() -> requests.Session:
    """
    Return a thread-local requests.Session (sessions are not thread-safe).
    """
    session = getattr(_thread_local, "session", None)
    if session is None:
        session = requests.Session()
        _thread_local.session = session
    return session


def http_get_text(
    url: str,
    retries: int = 3,
    retry_delay: float = 1.0,
    timeout: float = 60.0,
    verify_ssl: bool = True,
) -> str:
    """
    GET a URL and return the response text, retrying with exponential
    backoff on network errors and 5xx responses.

    Parameters
    ----------
    url : str
        Request URL.
    retries : int
        Maximum number of attempts.
    retry_delay : float
        Base delay in seconds (doubled every attempt).
    timeout : float
        Per-request timeout in seconds.
    verify_ssl : bool
        TLS certificate verification.

    Returns
    -------
    str
        Response body text.

    Raises
    ------
    RuntimeError
        When every attempt failed.
    """
    session = _get_session()
    last_error: Optional[Exception] = None

    for attempt in range(1, retries + 1):
        try:
            response = session.get(url, timeout=timeout, verify=verify_ssl)
            if response.status_code >= 500:
                raise RuntimeError(
                    f"HTTP {response.status_code}: {response.reason}"
                )
            if 400 <= response.status_code < 500 and response.status_code != 429:
                raise requests.HTTPError(
                    f"HTTP {response.status_code}: {response.reason}",
                    response=response,
                )
            response.raise_for_status()
            return response.text

        except Exception as exc:
            last_error = exc
            response = getattr(exc, "response", None)
            status = getattr(response, "status_code", None)
            deterministic = (
                isinstance(exc, requests.HTTPError)
                and status is not None
                and 400 <= status < 500
                and status != 429
            )
            if deterministic:
                logger.error(
                    "Request rejected (HTTP %s) %s: %s",
                    status, url, exc,
                )
                raise
            if attempt < retries:
                delay = retry_delay * (2 ** (attempt - 1))
                logger.warning(
                    "Request failed (%d/%d) %s: %s; retry in %.1fs",
                    attempt, retries, url, exc, delay,
                )
                time.sleep(delay)
            else:
                logger.error(
                    "Request failed (%d/%d) %s: %s",
                    attempt, retries, url, exc,
                )

    raise RuntimeError(f"GET failed after {retries} attempts: {url}: {last_error}")



def atomic_write_text(path: Path, text: str) -> None:
    """
    Write text to path atomically (tmp file + os.replace).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
    os.replace(tmp, path)


def clean_value(value: Optional[str]) -> str:
    """
    Strip a cell value and replace embedded tabs/newlines so the TSV
    stays rectangular.
    """
    if value is None:
        return ""
    return re.sub(r"[\t\r\n]+", " ", value).strip()


def read_cache(path: Path, overwrite: bool) -> Optional[str]:
    """
    Return cached file text, or None when absent/empty/overwrite requested.
    """
    if overwrite:
        return None
    if path.is_file() and path.stat().st_size > 0:
        logger.info("Using cached file: %s", path)
        return path.read_text(encoding="utf-8")
    return None



def fetch_sample_xml(
    srs: str,
    cache_dir: Path,
    retries: int = 3,
    timeout: float = 60.0,
    verify_ssl: bool = True,
    overwrite: bool = False,
) -> str:
    """
    Fetch the ENA sample XML for one SRS accession (cached as
    <cache_dir>/<SRS>.xml).
    """
    cache_path = cache_dir / f"{srs}.xml"
    cached = read_cache(cache_path, overwrite)
    if cached is not None:
        return cached

    url = ENA_SAMPLE_XML_URL.format(acc=srs)
    logger.info("Fetching sample XML: %s", srs)
    text = http_get_text(
        url, retries=retries, timeout=timeout, verify_ssl=verify_ssl
    )
    atomic_write_text(cache_path, text)
    return text


def fetch_run_tsv(
    srs: str,
    cache_dir: Path,
    retries: int = 3,
    timeout: float = 60.0,
    verify_ssl: bool = True,
    overwrite: bool = False,
) -> str:
    """
    Fetch the ENA run table (filereport) for one SRS accession (cached as
    <cache_dir>/<SRS>_runs.tsv).
    """
    cache_path = cache_dir / f"{srs}_runs.tsv"
    cached = read_cache(cache_path, overwrite)
    if cached is not None:
        return cached

    url = ENA_FILEREPORT_URL.format(acc=srs, fields=",".join(RUN_FIELDS))
    logger.info("Fetching run table: %s", srs)
    text = http_get_text(
        url, retries=retries, timeout=timeout, verify_ssl=verify_ssl
    )
    atomic_write_text(cache_path, text)
    return text



def parse_sample_xml(srs: str, xml_text: str) -> Dict[str, str]:
    """
    Parse one ENA sample XML document into a flat record.

    Parameters
    ----------
    srs : str
        The queried SRS accession (used for the `srs` column and errors).
    xml_text : str
        Raw XML returned by the ENA Browser API.

    Returns
    -------
    Dict[str, str]
        Record with keys srs, biosample, sample_alias, sample_title,
        organism, tax_id plus one key per SAMPLE_ATTRIBUTE tag. A repeated
        tag keeps both values joined with ';'. Attribute tags colliding
        with the base keys do not overwrite them.

    Raises
    ------
    ValueError
        When the document is not a SAMPLE_SET / SAMPLE document.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise ValueError(f"invalid XML for {srs}: {exc}") from exc

    sample = root.find("SAMPLE") if root.tag == "SAMPLE_SET" else root
    if sample is None or sample.tag != "SAMPLE":
        raise ValueError(f"no SAMPLE element in ENA response for {srs}")

    record: Dict[str, str] = {
        "srs": srs,
        "biosample": clean_value(sample.get("accession")),
        "sample_alias": clean_value(sample.get("alias")),
        "sample_title": clean_value(
            sample.findtext("TITLE", default="")
        ),
        "organism": clean_value(
            sample.findtext("SAMPLE_NAME/SCIENTIFIC_NAME", default="")
        ),
        "tax_id": clean_value(
            sample.findtext("SAMPLE_NAME/TAXON_ID", default="")
        ),
    }

    base_keys = set(record)
    for attribute in sample.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
        tag = clean_value(attribute.findtext("TAG", default=""))
        value = clean_value(attribute.findtext("VALUE", default=""))
        if not tag or tag in base_keys:
            continue
        if tag in record:
            record[tag] = f"{record[tag]};{value}" if record[tag] else value
        else:
            record[tag] = value

    return record


def parse_run_tsv(srs: str, tsv_text: str) -> List[Dict[str, str]]:
    """
    Parse one ENA filereport TSV into run records.

    Parameters
    ----------
    srs : str
        The queried SRS accession (used for the `srs` column and errors).
    tsv_text : str
        Raw TSV returned by the ENA Portal API (header + rows).

    Returns
    -------
    List[Dict[str, str]]
        One record per run with the returned header columns plus `srs`.
        An empty body (header only or blank) yields an empty list.

    Raises
    ------
    ValueError
        When the returned header columns differ from RUN_FIELDS
        (column order is free: the API hoists run_accession first).
    """
    text = tsv_text.strip()
    if not text:
        return []

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    header = reader.fieldnames or []
    if sorted(header) != sorted(RUN_FIELDS):
        raise ValueError(
            f"unexpected run-table header for {srs}: {header}"
        )

    rows: List[Dict[str, str]] = []
    for raw in reader:
        record: Dict[str, str] = {"srs": srs}
        for key, value in raw.items():
            record[key] = clean_value(value)
        rows.append(record)
    return rows



def parse_group_specs(specs: Optional[List[str]]) -> Dict[str, str]:
    """
    Parse -g/--group entries of the form LABEL=ACC1,ACC2,ACC3.

    Returns
    -------
    Dict[str, str]
        Mapping accession -> group label.

    Raises
    ------
    ValueError
        On a malformed entry or when one accession is assigned to two
        different groups.
    """
    mapping: Dict[str, str] = {}
    for spec in specs or []:
        if "=" not in spec:
            raise ValueError(
                f"malformed --group entry (want LABEL=ACC1,ACC2): {spec}"
            )
        label, accs_text = spec.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError(f"empty group label in --group entry: {spec}")
        for acc in accs_text.split(","):
            acc = acc.strip()
            if not acc:
                continue
            _assign_group(mapping, acc.upper(), label)
    return mapping


def parse_group_file(
    path: str,
    id_column: str,
    group_column: str,
) -> Dict[str, str]:
    """
    Parse a delimited group table (header required) into accession -> label.

    Parameters
    ----------
    path : str
        Table path (TSV/CSV; delimiter auto-detected).
    id_column : str
        Column holding the sample accession.
    group_column : str
        Column holding the group label.

    Raises
    ------
    ValueError
        When a required column is missing.
    """
    with open(path, "r", encoding="utf-8", newline="") as handle:
        sample_text = handle.read(8192)
        handle.seek(0)
        delimiter = detect_delimiter(sample_text)
        if delimiter == "whitespace":
            rows = [
                re.split(r"\s+", line.strip())
                for line in handle
                if line.strip()
            ]
            fieldnames = [name.strip() for name in rows[0]] if rows else []
            raw_records = [
                dict(zip(fieldnames, row)) for row in rows[1:]
            ]
        else:
            reader = csv.DictReader(handle, delimiter=delimiter)
            fieldnames = [name.strip() for name in (reader.fieldnames or [])]
            raw_records = list(reader)
        for wanted in (id_column, group_column):
            if wanted not in fieldnames:
                raise ValueError(
                    f"group file {path} has no column {wanted!r}; "
                    f"found: {fieldnames}"
                )
        mapping: Dict[str, str] = {}
        for raw in raw_records:
            acc = clean_value(raw.get(id_column)).upper()
            label = clean_value(raw.get(group_column))
            if not acc:
                continue
            if not label:
                logger.warning("Empty group label for %s in %s", acc, path)
                continue
            _assign_group(mapping, acc, label)
        return mapping


def _assign_group(mapping: Dict[str, str], acc: str, label: str) -> None:
    """
    Assign acc -> label, rejecting conflicting double assignments.
    """
    existing = mapping.get(acc)
    if existing is not None and existing != label:
        raise ValueError(
            f"{acc} assigned to both {existing!r} and {label!r}"
        )
    mapping[acc] = label



def write_sample_table(
    records: List[Dict[str, str]],
    output_tsv: Path,
) -> None:
    """
    Write sample_metadata.tsv: fixed base columns plus one column per
    sample-attribute tag, in first-seen order.
    """
    attr_keys: List[str] = []
    for record in records:
        for key in record:
            if key not in SAMPLE_BASE_FIELDS and key not in attr_keys:
                attr_keys.append(key)
    fieldnames = SAMPLE_BASE_FIELDS + attr_keys

    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=fieldnames,
        delimiter="\t",
        restval="",
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(records)
    atomic_write_text(output_tsv, buffer.getvalue())
    logger.info(
        "Wrote %d sample rows to %s", len(records), output_tsv
    )


def write_run_table(
    records: List[Dict[str, str]],
    output_tsv: Path,
) -> None:
    """
    Write run_metadata.tsv: `srs`, `group` plus the verbatim ENA
    filereport columns (in the order the API returned them).
    """
    if records:
        run_columns = [
            key for key in records[0] if key not in ("srs", "group")
        ]
    else:
        run_columns = list(RUN_FIELDS)
    fieldnames = ["srs", "group"] + run_columns

    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=fieldnames,
        delimiter="\t",
        restval="",
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(records)
    atomic_write_text(output_tsv, buffer.getvalue())
    logger.info("Wrote %d run rows to %s", len(records), output_tsv)



def collect_srs_metadata(
    srs_list: List[str],
    outdir: str,
    group_map: Optional[Dict[str, str]] = None,
    jobs: int = 4,
    retries: int = 3,
    timeout: float = 60.0,
    verify_ssl: bool = True,
    overwrite: bool = False,
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]], List[str]]:
    """
    Fetch and parse metadata for all SRS accessions.

    Workflow
    --------
    1. Fetch sample XML + run table per SRS in parallel (cached).
    2. Parse into sample records and run records.
    3. Write sample_metadata.tsv and run_metadata.tsv once at the end.

    Returns
    -------
    (sample_records, run_records, failed)
        Parsed rows and the accessions that failed (also written to
        <outdir>/srs_failed.txt by the caller).
    """
    group_map = group_map or {}
    outpath = Path(outdir)
    cache_dir = outpath / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    sample_records: Dict[str, Dict[str, str]] = {}
    run_records: List[Dict[str, str]] = []
    failed: List[str] = []
    errors: Dict[str, str] = {}

    def work(srs: str) -> Tuple[str, Dict[str, str], List[Dict[str, str]]]:
        xml_text = fetch_sample_xml(
            srs, cache_dir,
            retries=retries, timeout=timeout,
            verify_ssl=verify_ssl, overwrite=overwrite,
        )
        tsv_text = fetch_run_tsv(
            srs, cache_dir,
            retries=retries, timeout=timeout,
            verify_ssl=verify_ssl, overwrite=overwrite,
        )
        return srs, parse_sample_xml(srs, xml_text), parse_run_tsv(srs, tsv_text)

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
        futures = {executor.submit(work, srs): srs for srs in srs_list}
        for future in as_completed(futures):
            srs = futures[future]
            try:
                _, sample, runs = future.result()
            except Exception as exc:
                logger.error("Failed to collect metadata for %s: %s", srs, exc)
                failed.append(srs)
                errors[srs] = str(exc)
                continue

            group = group_map.get(srs, "")
            sample["group"] = group
            sample_records[srs] = sample
            for run in runs:
                run["group"] = group
            run_records.extend(runs)
            if not runs:
                logger.warning(
                    "No public runs found for %s (sample metadata kept)",
                    srs,
                )
            logger.info(
                "Collected %s: 1 sample row, %d run rows", srs, len(runs)
            )

    ordered_samples = [
        sample_records[srs] for srs in srs_list if srs in sample_records
    ]
    run_order = {srs: i for i, srs in enumerate(srs_list)}
    run_records.sort(key=lambda row: run_order.get(row.get("srs", ""), 0))

    write_sample_table(ordered_samples, outpath / "sample_metadata.tsv")
    write_run_table(run_records, outpath / "run_metadata.tsv")

    failed = sorted(failed)
    if failed:
        failed_path = outpath / "srs_failed.txt"
        atomic_write_text(
            failed_path, "".join(f"{acc}\t{errors[acc]}\n" for acc in failed)
        )
        logger.error(
            "%d accessions failed, listed in %s", len(failed), failed_path
        )
    else:
        stale = outpath / "srs_failed.txt"
        if stale.exists():
            stale.unlink()
            logger.info("Removed stale failed list: %s", stale)

    return ordered_samples, run_records, failed



def parser_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Fetch sample- and run-level metadata for SRA sample "
            "accessions (SRS...) from ENA into consolidated TSV tables."
        )
    )
    parser.add_argument(
        "-s", "--srs",
        type=str,
        nargs="+",
        action="extend",
        default=None,
        help="Explicit SRS accession list (repeatable).",
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        default=None,
        help=(
            "File with SRS accessions: one per line, or a delimited "
            "table with an ID column (--id-column)."
        ),
    )
    parser.add_argument(
        "-c", "--id-column",
        type=str,
        default="SRS",
        help="ID column name used when --input is a table. Default: SRS.",
    )
    parser.add_argument(
        "-g", "--group",
        type=str,
        nargs="+",
        action="extend",
        default=None,
        help=(
            "Group labels as LABEL=ACC1,ACC2,... (repeatable). "
            "Accessions named here are collected as well. "
            "Example: -g experiment=SRS1,SRS2 -g control=SRS3"
        ),
    )
    parser.add_argument(
        "-G", "--group-file",
        type=str,
        default=None,
        help="Delimited table (header required) with sample + group columns.",
    )
    parser.add_argument(
        "--group-id-col",
        type=str,
        default="SRS",
        help="Accession column of --group-file. Default: SRS.",
    )
    parser.add_argument(
        "--group-name-col",
        type=str,
        default="group",
        help="Group label column of --group-file. Default: group.",
    )
    parser.add_argument(
        "-o", "--outdir",
        type=str,
        required=True,
        help="Output directory for the TSV tables and the cache/ files.",
    )
    parser.add_argument(
        "-l", "--log",
        type=str,
        default=None,
        help="Optional log file (logs always go to stdout as well).",
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=4,
        help="Number of parallel fetch workers. Default: 4.",
    )
    parser.add_argument(
        "-r", "--retries",
        type=int,
        default=3,
        help="Retry attempts per HTTP request. Default: 3.",
    )
    parser.add_argument(
        "-t", "--timeout",
        type=float,
        default=60.0,
        help="Per-request timeout in seconds. Default: 60.",
    )
    parser.add_argument(
        "-w", "--overwrite",
        action="store_true",
        help="Re-fetch even when a cached file exists (default: reuse cache).",
    )
    parser.add_argument(
        "-n", "--dry-run",
        action="store_true",
        help="Log the resolved accessions, groups and output paths; no network.",
    )
    parser.add_argument(
        "--insecure",
        action="store_true",
        help="Disable TLS certificate verification.",
    )
    return parser.parse_args()



def main() -> None:
    """
    Main CLI entry point.
    """
    args = parser_args()

    global logger
    log_file = str(args.log) if args.log else None
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(
        "SRS_metadata", log_file=log_file, stream=sys.stdout
    )

    group_map = parse_group_specs(args.group)
    if args.group_file:
        group_map.update(
            parse_group_file(
                args.group_file,
                id_column=args.group_id_col,
                group_column=args.group_name_col,
            )
        )

    resolver = AccessionResolver(
        acc=(args.srs or []) + list(group_map.keys()),
        acc_file=args.input,
        id_column=args.id_column,
    )
    accs = {acc.upper() for acc in resolver.resolve()}

    bad = sorted(acc for acc in accs if not SRS_PATTERN.match(acc))
    if bad:
        logger.error(
            "Not valid SRS accessions (this tool only queries SRS): %s",
            ", ".join(bad),
        )
        sys.exit(2)
    srs_list = sorted(acc.upper() for acc in accs)

    unknown_groups = sorted(set(group_map) - set(srs_list))
    if unknown_groups:
        logger.warning(
            "Group labels name accessions outside the resolved list: %s",
            ", ".join(unknown_groups),
        )
    unassigned = [srs for srs in srs_list if srs not in group_map]
    if unassigned:
        logger.warning(
            "No group label for %d accessions (empty group column): %s",
            len(unassigned),
            ", ".join(unassigned),
        )
    label_counts: Dict[str, int] = {}
    for label in group_map.values():
        label_counts[label] = label_counts.get(label, 0) + 1

    outdir = str(Path(args.outdir).resolve())

    logger.info("Resolved %d SRS accessions", len(srs_list))
    for label, count in sorted(label_counts.items()):
        logger.info("Group %s: %d accessions", label, count)

    if args.dry_run:
        logger.info("Dry run: planned outputs in %s", outdir)
        logger.info("  sample_metadata.tsv, run_metadata.tsv, cache/<SRS>.xml, cache/<SRS>_runs.tsv")
        for srs in srs_list:
            logger.info(
                "  %s group=%s xml=%s runs=%s",
                srs,
                group_map.get(srs, ""),
                ENA_SAMPLE_XML_URL.format(acc=srs),
                ENA_FILEREPORT_URL.format(acc=srs, fields=",".join(RUN_FIELDS)),
            )
        return

    os.makedirs(outdir, exist_ok=True)

    _, _, failed = collect_srs_metadata(
        srs_list=srs_list,
        outdir=outdir,
        group_map=group_map,
        jobs=args.jobs,
        retries=args.retries,
        timeout=args.timeout,
        verify_ssl=not args.insecure,
        overwrite=args.overwrite,
    )

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
