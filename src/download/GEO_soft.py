#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


USER_AGENT = (
    "Mozilla/5.0 "
    "(compatible; GEO-SOFT-Metadata-Downloader/1.5)"
)

GEO_FAMILY_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
)

EUTILS_BASE = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
)

DEFAULT_TIMEOUT = (20, 120)
SRA_BATCH_SIZE = 100


def normalize_accession(accession: str) -> str:
    """Normalize and validate a GEO Series accession."""
    accession = accession.strip().upper()

    if not re.fullmatch(r"GSE\d+", accession):
        raise ValueError(
            f"Invalid GEO Series accession: {accession}"
        )

    return accession


def build_family_url(accession: str) -> str:
    """Build the official GEO family SOFT URL."""
    accession = normalize_accession(accession)
    digits = accession[3:]

    if len(digits) <= 3:
        group = "GSEnnn"
    else:
        group = f"GSE{digits[:-3]}nnn"

    return (
        f"{GEO_FAMILY_BASE}"
        f"{group}/{accession}/soft/"
        f"{accession}_family.soft.gz"
    )


def create_session(
    proxy: Optional[str],
    retries: int = 3,
) -> requests.Session:
    """Create a configured HTTP session."""
    session = requests.Session()

    session.headers.update({
        "User-Agent": USER_AGENT,
        "Accept": "*/*",
    })

    if proxy:
        session.proxies.update({
            "http": proxy,
            "https": proxy,
        })

    retry = Retry(
        total=retries,
        connect=retries,
        read=retries,
        status=retries,
        backoff_factor=1.0,
        status_forcelist=(
            429,
            500,
            502,
            503,
            504,
        ),
        allowed_methods=frozenset(["GET"]),
        raise_on_status=False,
    )

    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=10,
        pool_maxsize=10,
    )

    session.mount("http://", adapter)
    session.mount("https://", adapter)

    return session


def download_file(
    url: str,
    output_path: Path,
    proxy: Optional[str] = None,
    overwrite: bool = False,
    retries: int = 3,
    timeout: Tuple[int, int] = DEFAULT_TIMEOUT,
) -> Path:
    """Download a file with proxy, retries, and local caching."""
    if output_path.exists() and not overwrite:
        print(f"[CACHE] {output_path}")
        return output_path

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    session = create_session(
        proxy=proxy,
        retries=retries,
    )

    temp_path = output_path.with_suffix(
        output_path.suffix + ".tmp"
    )

    try:
        print(f"[DOWNLOAD] {url}")

        if proxy:
            print(f"[PROXY] {proxy}")

        with session.get(
            url,
            stream=True,
            timeout=timeout,
        ) as response:
            response.raise_for_status()

            total = response.headers.get(
                "Content-Length"
            )

            total_size = (
                int(total)
                if total is not None
                else None
            )

            downloaded = 0

            with temp_path.open("wb") as handle:
                for chunk in response.iter_content(
                    chunk_size=1024 * 1024
                ):
                    if not chunk:
                        continue

                    handle.write(chunk)
                    downloaded += len(chunk)

                    if total_size:
                        percent = (
                            downloaded
                            / total_size
                            * 100
                        )

                        print(
                            f"\r[DOWNLOAD] "
                            f"{downloaded / 1024**2:.1f} / "
                            f"{total_size / 1024**2:.1f} MB "
                            f"({percent:.1f}%)",
                            end="",
                            flush=True,
                        )

            if total_size:
                print()

        temp_path.replace(output_path)

        return output_path

    except Exception:
        if temp_path.exists():
            temp_path.unlink()

        raise

    finally:
        session.close()


def test_proxy(
    proxy: str,
    timeout: int = 10,
) -> None:
    """Test proxy connectivity."""
    print(f"[INFO] Testing proxy: {proxy}")

    session = create_session(proxy)

    try:
        response = session.get(
            "https://www.ncbi.nlm.nih.gov/",
            timeout=timeout,
        )

        print(
            f"[INFO] Proxy test successful: "
            f"HTTP {response.status_code}"
        )

    except Exception as exc:
        raise RuntimeError(
            f"Proxy test failed: {exc}"
        ) from exc

    finally:
        session.close()


def parse_soft_value(
    line: str,
) -> Tuple[str, str]:
    """Parse a GEO SOFT metadata line."""
    if " = " not in line:
        raise ValueError(
            f"Invalid SOFT metadata line: {line}"
        )

    key, value = line.split(
        " = ",
        1,
    )

    return key, value


def append_field(
    record: Dict[str, str],
    key: str,
    value: str,
) -> None:
    """Append a SOFT field while preserving repeated values."""
    if key not in record:
        record[key] = value
    else:
        record[key] += " | " + value


def parse_soft(
    soft_path: Path,
) -> Dict[str, List[Dict[str, str]]]:
    """Parse Series, Sample, and Platform records."""
    records = {
        "SERIES": [],
        "SAMPLE": [],
        "PLATFORM": [],
    }

    current_type: Optional[str] = None
    current_record: Optional[Dict[str, str]] = None

    def flush() -> None:
        nonlocal current_type
        nonlocal current_record

        if (
            current_type in records
            and current_record is not None
        ):
            records[current_type].append(
                current_record
            )

        current_type = None
        current_record = None

    with gzip.open(
        soft_path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for raw_line in handle:
            line = raw_line.rstrip("\r\n")

            if line.startswith("^"):
                flush()

                match = re.match(
                    r"^\^([A-Z]+)\s*=\s*(.+)$",
                    line,
                )

                if match is None:
                    continue

                current_type = match.group(1)

                if current_type not in records:
                    current_type = None
                    continue

                current_record = {
                    "_accession": (
                        match.group(2).strip()
                    )
                }

                continue

            if (
                current_type is None
                or current_record is None
            ):
                continue

            if not line.startswith("!"):
                continue

            try:
                key, value = parse_soft_value(
                    line[1:]
                )
            except ValueError:
                continue

            append_field(
                current_record,
                key,
                value,
            )

    flush()

    return records


def normalize_column_name(
    name: str,
) -> str:
    """Normalize a GEO field name."""
    name = re.sub(
        r"^!+",
        "",
        name,
    )

    name = name.lower()

    name = re.sub(
        r"[^a-z0-9]+",
        "_",
        name,
    )

    return name.strip("_")


def expand_sample_characteristics(
    record: Dict[str, str],
) -> Dict[str, str]:
    """Expand GEO sample characteristics into columns."""
    result = dict(record)

    keys = [
        key
        for key in record
        if key.startswith(
            "Sample_characteristics_ch"
        )
    ]

    for field in keys:
        for item in record[field].split(" | "):
            if ":" not in item:
                continue

            key, value = item.split(
                ":",
                1,
            )

            key = key.strip()
            value = value.strip()

            if not key:
                continue

            column = normalize_column_name(key)

            if column in result:
                result[column] += (
                    " | " + value
                )
            else:
                result[column] = value

    return result


def extract_accessions(
    text: str,
) -> Dict[str, List[str]]:
    """Extract GEO, SRA, BioSample, and BioProject accessions."""
    patterns = {
        "gsm": r"\bGSM\d+\b",
        "biosample": r"\bSAM(?:N|D|E)\d+\b",
        "srx": r"\b(?:SRX|DRX|ERX)\d+\b",
        "srr": r"\b(?:SRR|DRR|ERR)\d+\b",
        "sra_study": r"\b(?:SRP|DRP|ERP)\d+\b",
        "bioproject": r"\bPRJ(?:NA|EB|DB)\d+\b",
    }

    result: Dict[str, List[str]] = {}

    text = text.upper()

    for name, pattern in patterns.items():
        values = re.findall(
            pattern,
            text,
        )

        result[name] = list(
            dict.fromkeys(values)
        )

    return result


def get_sample_relations(
    record: Dict[str, str],
) -> Dict[str, List[str]]:
    """Extract external accessions from GEO sample relations."""
    relation_text = " ".join(
        value
        for key, value in record.items()
        if key.lower().startswith(
            "sample_relation"
        )
    )

    return extract_accessions(
        relation_text
    )


def prepare_records(
    records: Sequence[Dict[str, str]],
    expand_characteristics: bool = False,
    add_relations: bool = False,
    accession_column: Optional[str] = None,
) -> List[Dict[str, str]]:
    """Prepare GEO records for tabular output."""
    prepared = []

    for record in records:

        if expand_characteristics:
            record = expand_sample_characteristics(
                record
            )

        record = dict(record)

        if add_relations:
            relations = get_sample_relations(
                record
            )

            record["GEO_GSM"] = (
                " | ".join(
                    relations["gsm"]
                )
            )

            record["BioSample"] = (
                " | ".join(
                    relations["biosample"]
                )
            )

            record["SRA_Experiment"] = (
                " | ".join(
                    relations["srx"]
                )
            )

            record["SRA_Study"] = (
                " | ".join(
                    relations["sra_study"]
                )
            )

            record["BioProject"] = (
                " | ".join(
                    relations["bioproject"]
                )
            )

        normalized = {
            normalize_column_name(key): value
            for key, value in record.items()
        }

        if accession_column:
            accession = normalized.pop(
                "accession",
                "",
            )

            normalized[
                accession_column
            ] = accession

        prepared.append(normalized)

    return prepared


def collect_columns(
    records: Sequence[Dict[str, str]],
) -> List[str]:
    """Collect columns in first-seen order."""
    columns = []
    seen = set()

    for record in records:
        for key in record:
            if key not in seen:
                columns.append(key)
                seen.add(key)

    return columns


def write_tsv(
    records: Sequence[Dict[str, str]],
    output_path: Path,
) -> None:
    """Write records to a TSV file."""
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not records:
        output_path.write_text(
            "",
            encoding="utf-8",
        )
        return

    columns = collect_columns(records)

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(columns)

        for record in records:
            writer.writerow([
                record.get(column, "")
                for column in columns
            ])


def read_csv(
    path: Path,
) -> List[Dict[str, str]]:
    """Read a CSV file into dictionaries."""
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
        errors="replace",
    ) as handle:

        reader = csv.DictReader(handle)

        return [
            dict(row)
            for row in reader
        ]


def normalize_sra_row(
    row: Dict[str, str],
) -> Dict[str, str]:
    """Normalize SRA RunInfo fields."""
    aliases = {
        "srr": ["Run"],
        "release_date": ["ReleaseDate"],
        "load_date": ["LoadDate"],
        "spots": ["spots"],
        "bases": ["bases"],
        "spots_with_mates": ["spots_with_mates"],
        "avg_length": ["avgLength"],
        "size_mb": ["size_MB"],
        "experiment": ["Experiment"],
        "library_name": ["LibraryName"],
        "library_strategy": ["LibraryStrategy"],
        "library_selection": ["LibrarySelection"],
        "library_source": ["LibrarySource"],
        "library_layout": ["LibraryLayout"],
        "insert_size": ["InsertSize"],
        "insert_dev": ["InsertDev"],
        "platform": ["Platform"],
        "instrument_model": ["Model"],
        "sra_study": ["SRAStudy"],
        "bioproject": ["BioProject"],
        "biosample": [
            "BioSample",
            "BioSampleModel",
        ],
        "sample_name": [
            "SampleName",
            "sample_name",
        ],
        "taxid": ["taxid"],
        "scientific_name": ["ScientificName"],
        "fastq_ftp": ["fastq_ftp"],
        "fastq_aspera": ["fastq_aspera"],
        "fastq_bytes": ["fastq_bytes"],
        "fastq_md5": ["fastq_md5"],
        "download_path": ["download_path"],
    }

    normalized: Dict[str, str] = {}

    for target, candidates in aliases.items():
        value = ""

        for candidate in candidates:
            candidate_value = row.get(candidate)

            if candidate_value:
                value = candidate_value
                break

        normalized[target] = value

    normalized["_raw"] = row

    return normalized


def build_srx_to_gsm(
    gsm_records: Sequence[Dict[str, str]],
) -> Dict[str, List[str]]:
    """Build SRA experiment-to-GSM mapping."""
    mapping: Dict[str, List[str]] = {}

    for record in gsm_records:
        gsm = record.get(
            "gsm",
            "",
        ).strip().upper()

        if not gsm:
            continue

        sra_experiments = (
            record.get(
                "sra_experiment",
                "",
            )
        )

        relations = extract_accessions(
            sra_experiments
        )

        for srx in relations["srx"]:
            mapping.setdefault(
                srx,
                [],
            )

            if gsm not in mapping[srx]:
                mapping[srx].append(gsm)

    return mapping


def build_biosample_to_gsm(
    gsm_records: Sequence[Dict[str, str]],
) -> Dict[str, List[str]]:
    """Build BioSample-to-GSM mapping."""
    mapping: Dict[str, List[str]] = {}

    for record in gsm_records:
        gsm = record.get(
            "gsm",
            "",
        ).strip().upper()

        if not gsm:
            continue

        biosamples = extract_accessions(
            record.get(
                "biosample",
                "",
            )
        )["biosample"]

        for biosample in biosamples:
            mapping.setdefault(
                biosample,
                [],
            )

            if gsm not in mapping[biosample]:
                mapping[biosample].append(
                    gsm
                )

    return mapping


def build_gsm_metadata(
    gsm_records: Sequence[Dict[str, str]],
) -> Dict[str, Dict[str, str]]:
    """Build GSM-to-GEO metadata mapping."""
    mapping: Dict[str, Dict[str, str]] = {}

    for record in gsm_records:
        gsm = record.get(
            "gsm",
            "",
        ).strip().upper()

        if not gsm:
            continue

        mapping[gsm] = dict(record)

    return mapping


def attach_gsm_to_sra(
    sra_rows: Sequence[Dict[str, str]],
    gsm_records: Sequence[Dict[str, str]],
) -> List[Dict[str, str]]:
    """
    Attach GSM accessions and GEO experimental metadata to SRA rows.

    SRA Experiment is used as the primary relationship.
    BioSample is used as a fallback.

    All GSM metadata is added with the ``gsm_`` prefix
    to avoid collisions with SRA fields.
    """
    srx_to_gsm = build_srx_to_gsm(
        gsm_records
    )

    biosample_to_gsm = build_biosample_to_gsm(
        gsm_records
    )

    gsm_metadata = build_gsm_metadata(
        gsm_records
    )

    result = []

    for row in sra_rows:
        row = dict(row)

        raw_row = row.pop(
            "_raw",
            {},
        )

        experiment = row.get(
            "experiment",
            "",
        ).strip().upper()

        gsm_values: List[str] = []

        if experiment:
            gsm_values.extend(
                srx_to_gsm.get(
                    experiment,
                    [],
                )
            )

        if not gsm_values:
            biosample_values = extract_accessions(
                row.get(
                    "biosample",
                    "",
                )
            )["biosample"]

            for biosample in biosample_values:
                gsm_values.extend(
                    biosample_to_gsm.get(
                        biosample,
                        [],
                    )
                )

        if not gsm_values:
            raw_text = "\t".join(
                str(value)
                for value in raw_row.values()
                if value
            )

            detected = extract_accessions(
                raw_text
            )

            gsm_values.extend(
                detected["gsm"]
            )

            if not gsm_values:
                for biosample in detected["biosample"]:
                    gsm_values.extend(
                        biosample_to_gsm.get(
                            biosample,
                            [],
                        )
                    )

        gsm_values = list(
            dict.fromkeys(
                value.strip().upper()
                for value in gsm_values
                if value.strip()
            )
        )

        row["gsm"] = ";".join(
            gsm_values
        )

        if len(gsm_values) == 1:
            gsm_record = gsm_metadata.get(
                gsm_values[0]
            )

            if gsm_record:
                for key, value in gsm_record.items():
                    if key == "gsm":
                        continue

                    if key.startswith("sra_"):
                        continue

                    row[
                        f"gsm_{key}"
                    ] = value

        elif gsm_values:
            merged: Dict[str, List[str]] = {}

            for gsm in gsm_values:
                gsm_record = gsm_metadata.get(
                    gsm
                )

                if not gsm_record:
                    continue

                for key, value in gsm_record.items():
                    if (
                        key == "gsm"
                        or key.startswith("sra_")
                    ):
                        continue

                    if not value:
                        continue

                    merged.setdefault(
                        key,
                        [],
                    )

                    if value not in merged[key]:
                        merged[key].append(
                            value
                        )

            for key, values in merged.items():
                row[
                    f"gsm_{key}"
                ] = " | ".join(values)

        result.append(row)

    return result


def fetch_sra_runinfo_batch(
    session: requests.Session,
    srx_accessions: Sequence[str],
) -> List[Dict[str, str]]:
    """
    Fetch SRA RunInfo for a batch of SRA experiment accessions.

    One SRX may correspond to multiple SRR runs.
    """
    if not srx_accessions:
        return []

    params = {
        "db": "sra",
        "id": ",".join(srx_accessions),
        "rettype": "runinfo",
        "retmode": "text",
    }

    response = session.get(
        f"{EUTILS_BASE}/efetch.fcgi",
        params=params,
        timeout=DEFAULT_TIMEOUT,
    )

    response.raise_for_status()

    text = response.text

    if not text.strip():
        return []

    lines = text.splitlines()

    if not lines:
        return []

    reader = csv.DictReader(
        lines
    )

    return [
        dict(row)
        for row in reader
    ]


def query_sra_runinfo(
    srx_accessions: Sequence[str],
    output_path: Path,
    proxy: Optional[str],
    overwrite: bool = False,
) -> List[Dict[str, str]]:
    """
    Query and cache SRA RunInfo using GEO-linked SRX accessions.

    GEO GSM -> SRA Experiment (SRX) -> SRA RunInfo -> SRR
    """
    if output_path.exists() and not overwrite:
        print(f"[CACHE] {output_path}")

        rows = read_csv(
            output_path
        )

        return [
            normalize_sra_row(row)
            for row in rows
        ]

    unique_srx = list(
        dict.fromkeys(
            accession.strip().upper()
            for accession in srx_accessions
            if accession.strip()
        )
    )

    if not unique_srx:
        print(
            "[SRA] No GEO-linked SRX accessions"
        )

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        output_path.write_text(
            "",
            encoding="utf-8",
        )

        return []

    print(
        f"[SRA] Query {len(unique_srx)} "
        f"SRA experiments"
    )

    session = create_session(
        proxy=proxy
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temp_path = output_path.with_suffix(
        output_path.suffix + ".tmp"
    )

    all_rows: List[Dict[str, str]] = []

    try:
        total_batches = (
            len(unique_srx)
            + SRA_BATCH_SIZE
            - 1
        ) // SRA_BATCH_SIZE

        for batch_index, start in enumerate(
            range(
                0,
                len(unique_srx),
                SRA_BATCH_SIZE,
            ),
            start=1,
        ):
            batch = unique_srx[
                start:start + SRA_BATCH_SIZE
            ]

            print(
                f"[SRA] Fetch batch "
                f"{batch_index}/{total_batches}: "
                f"{len(batch)} SRX"
            )

            rows = fetch_sra_runinfo_batch(
                session,
                batch,
            )

            all_rows.extend(rows)

        if all_rows:
            fieldnames = list(
                all_rows[0].keys()
            )

            with temp_path.open(
                "w",
                encoding="utf-8",
                newline="",
            ) as handle:

                writer = csv.DictWriter(
                    handle,
                    fieldnames=fieldnames,
                    extrasaction="ignore",
                )

                writer.writeheader()
                writer.writerows(all_rows)

            temp_path.replace(
                output_path
            )

            print(
                f"[SRA] Retrieved "
                f"{len(all_rows)} runs"
            )

        else:
            print(
                "[SRA] No RunInfo returned"
            )

            output_path.write_text(
                "",
                encoding="utf-8",
            )

    finally:
        if temp_path.exists():
            temp_path.unlink()

        session.close()

    return [
        normalize_sra_row(row)
        for row in all_rows
    ]


SRA_PARSED_BASE_COLUMNS = [
    "gse",
    "gsm",
    "experiment",
    "srr",
    "sra_study",
    "biosample",
    "bioproject",
    "library_name",
    "library_strategy",
    "library_selection",
    "library_source",
    "library_layout",
    "insert_size",
    "insert_dev",
    "platform",
    "instrument_model",
    "sample_name",
    "taxid",
    "scientific_name",
    "release_date",
    "load_date",
    "spots",
    "bases",
    "spots_with_mates",
    "avg_length",
    "size_mb",
    "fastq_ftp",
    "fastq_aspera",
    "fastq_bytes",
    "fastq_md5",
]


def collect_gsm_columns(
    rows: Sequence[Dict[str, str]],
) -> List[str]:
    """Collect integrated GSM metadata columns."""
    columns = []
    seen = set()

    for row in rows:
        for key in row:
            if not key.startswith("gsm_"):
                continue

            if key not in seen:
                columns.append(key)
                seen.add(key)

    return columns


def write_sra_tsv(
    rows: Sequence[Dict[str, str]],
    output_path: Path,
    gse: str,
) -> None:
    """
    Write integrated GEO-SRA Run records.

    The output granularity is one row per SRA Run (SRR).
    Each row contains SRA RunInfo together with the
    corresponding GSM experimental metadata.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    gsm_columns = collect_gsm_columns(
        rows
    )

    columns = (
        SRA_PARSED_BASE_COLUMNS
        + gsm_columns
    )

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            output = dict(row)
            output["gse"] = gse

            writer.writerow({
                column: output.get(
                    column,
                    "",
                )
                for column in columns
            })


def build_gsm_sra_summary(
    gsm_records: Sequence[Dict[str, str]],
    sra_rows: Sequence[Dict[str, str]],
) -> List[Dict[str, str]]:
    """Add compact SRA information to GEO GSM records."""
    mapping: Dict[
        str,
        Dict[str, List[str]]
    ] = {}

    for row in sra_rows:
        gsm_values = [
            value.strip().upper()
            for value in row.get(
                "gsm",
                "",
            ).split(";")
            if value.strip()
        ]

        for gsm in gsm_values:
            item = mapping.setdefault(
                gsm,
                {
                    "sra_runs": [],
                    "sra_experiments": [],
                },
            )

            srr = row.get(
                "srr",
                "",
            ).strip().upper()

            if (
                srr
                and srr not in item["sra_runs"]
            ):
                item["sra_runs"].append(
                    srr
                )

            experiment = row.get(
                "experiment",
                "",
            ).strip().upper()

            if (
                experiment
                and experiment
                not in item["sra_experiments"]
            ):
                item[
                    "sra_experiments"
                ].append(
                    experiment
                )

    result = []

    for record in gsm_records:
        row = dict(record)

        gsm = row.get(
            "gsm",
            "",
        ).strip().upper()

        item = mapping.get(
            gsm,
            {
                "sra_runs": [],
                "sra_experiments": [],
            },
        )

        row["sra_experiments"] = ";".join(
            item["sra_experiments"]
        )

        row["sra_runs"] = ";".join(
            item["sra_runs"]
        )

        row["sra_run_count"] = str(
            len(item["sra_runs"])
        )

        result.append(row)

    return result


def process_one(
    accession: str,
    output_dir: Path,
    proxy: Optional[str],
    keep_soft: bool,
    overwrite: bool,
    query_sra: bool,
) -> Tuple[str, bool, Optional[str]]:
    """Download and parse one GEO Series."""
    try:
        accession = normalize_accession(
            accession
        )

        raw_dir = output_dir / "raw"
        parsed_dir = output_dir / "parsed"

        raw_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        parsed_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        url = build_family_url(
            accession
        )

        soft_path = (
            raw_dir
            / f"{accession}_family.soft.gz"
        )

        download_file(
            url,
            soft_path,
            proxy=proxy,
            overwrite=overwrite,
        )

        records = parse_soft(
            soft_path
        )

        gse_records = prepare_records(
            records["SERIES"],
            accession_column="gse",
        )

        write_tsv(
            gse_records,
            parsed_dir
            / f"{accession}_GSE.tsv",
        )

        gsm_records = prepare_records(
            records["SAMPLE"],
            expand_characteristics=True,
            add_relations=True,
            accession_column="gsm",
        )

        write_tsv(
            gsm_records,
            parsed_dir
            / f"{accession}_GSM.tsv",
        )

        gpl_records = prepare_records(
            records["PLATFORM"],
            accession_column="gpl",
        )

        write_tsv(
            gpl_records,
            parsed_dir
            / f"{accession}_GPL.tsv",
        )

        if query_sra:
            srx_accessions = []

            for sample in gsm_records:
                sra_experiments = sample.get(
                    "sra_experiment",
                    "",
                )

                relations = extract_accessions(
                    sra_experiments
                )

                srx_accessions.extend(
                    relations["srx"]
                )

            srx_accessions = list(
                dict.fromkeys(
                    srx_accessions
                )
            )

            print(
                f"[SRA] {accession}: "
                f"{len(gsm_records)} GSM, "
                f"{len(srx_accessions)} SRX"
            )

            sra_cache = (
                raw_dir
                / f"{accession}_SRA_RunInfo.csv"
            )

            sra_rows = query_sra_runinfo(
                srx_accessions=srx_accessions,
                output_path=sra_cache,
                proxy=proxy,
                overwrite=overwrite,
            )

            sra_rows = attach_gsm_to_sra(
                sra_rows,
                gsm_records,
            )

            write_sra_tsv(
                sra_rows,
                parsed_dir
                / f"{accession}_SRA.tsv",
                accession,
            )

            gsm_records = (
                build_gsm_sra_summary(
                    gsm_records,
                    sra_rows,
                )
            )

            write_tsv(
                gsm_records,
                parsed_dir
                / f"{accession}_GSM.tsv",
            )

        if not keep_soft:
            soft_path.unlink()

        return (
            accession,
            True,
            None,
        )

    except Exception as exc:
        return (
            accession,
            False,
            str(exc),
        )


def read_accessions(
    values: Sequence[str],
    input_file: Optional[Path],
) -> List[str]:
    """Read and deduplicate GEO accessions."""
    all_values = list(values)

    if input_file:
        with input_file.open(
            "r",
            encoding="utf-8",
        ) as handle:

            all_values.extend(
                line.strip()
                for line in handle
                if line.strip()
                and not line.lstrip().startswith("#")
            )

    result = []
    seen = set()

    for value in all_values:
        accession = normalize_accession(
            value
        )

        if accession not in seen:
            result.append(accession)
            seen.add(accession)

    return result


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Download and parse GEO family SOFT files "
            "and associated SRA RunInfo metadata."
        )
    )

    parser.add_argument(
        "accessions",
        nargs="*",
        help="GEO Series accessions.",
    )

    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        help=(
            "Text file containing one GSE accession "
            "per line."
        ),
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("geo_output"),
    )

    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=2,
        help=(
            "Number of GEO Series processed concurrently. "
            "Default: 2."
        ),
    )

    parser.add_argument(
        "--proxy",
        default=None,
        help=(
            "HTTP/HTTPS proxy, e.g. "
            "http://127.0.0.1:7890"
        ),
    )

    parser.add_argument(
        "--test-proxy",
        action="store_true",
    )

    parser.add_argument(
        "--keep-soft",
        action="store_true",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    parser.add_argument(
        "--no-sra",
        action="store_true",
        help="Disable SRA RunInfo retrieval.",
    )

    return parser.parse_args()


def main() -> int:
    """Run GEO metadata ingestion."""
    args = parse_args()

    if args.threads < 1:
        raise SystemExit(
            "--threads must be >= 1"
        )

    proxy = args.proxy

    if args.test_proxy:
        if not proxy:
            raise SystemExit(
                "--test-proxy requires --proxy"
            )

        test_proxy(proxy)

    accessions = read_accessions(
        args.accessions,
        args.input,
    )

    if not accessions:
        raise SystemExit(
            "No GEO accession provided."
        )

    print(
        f"[INFO] Accessions: "
        f"{len(accessions)}"
    )

    print(
        f"[INFO] Threads: "
        f"{args.threads}"
    )

    print(
        f"[INFO] SRA: "
        f"{'disabled' if args.no_sra else 'enabled'}"
    )

    failures = []

    with ThreadPoolExecutor(
        max_workers=args.threads
    ) as executor:

        futures = {
            executor.submit(
                process_one,
                accession,
                args.output_dir,
                proxy,
                args.keep_soft,
                args.overwrite,
                not args.no_sra,
            ): accession
            for accession in accessions
        }

        for future in as_completed(
            futures
        ):
            accession = futures[future]

            try:
                (
                    accession,
                    success,
                    error,
                ) = future.result()

            except Exception as exc:
                success = False
                error = str(exc)

            if success:
                print(
                    f"[OK] {accession}"
                )

            else:
                print(
                    f"[ERROR] {accession}: "
                    f"{error}",
                    file=sys.stderr,
                )

                failures.append(
                    accession
                )

    print(
        f"[INFO] Success: "
        f"{len(accessions) - len(failures)}"
    )

    print(
        f"[INFO] Failed: "
        f"{len(failures)}"
    )

    if failures:
        print(
            "[INFO] Failed accessions:"
        )

        for accession in failures:
            print(accession)

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

