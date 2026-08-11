#!/usr/bin/env python3

import argparse
import csv
import json
import logging
import os
import random
import re
import subprocess
import time
import xml.etree.ElementTree as ET

from pathlib import Path
from typing import Any, Callable, Dict, List, Literal, Optional

import pandas as pd

try:
    from Bio import Entrez
    from urllib.error import HTTPError, URLError
except ImportError:
    Entrez = None
    HTTPError = URLError = Exception


# ============================================================================
# Logging
# ============================================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logger = logging.getLogger(__name__)


# ============================================================================
# Type definitions
# ============================================================================

EDirectDatabase = Literal["gds", "sra"]
EDirectFormat = Literal["json", "docsum", "runinfo"]

Record = Dict[str, str]


# ============================================================================
# General utilities
# ============================================================================

def get_text_value(
    element: Optional[ET.Element],
    tag: str,
    default: str = "N/A",
) -> str:
    """
    Safely extract text from an XML child element.

    Parameters
    ----------
    element:
        Parent XML element.
    tag:
        Child element tag.
    default:
        Value returned when the element or its text does not exist.

    Returns
    -------
    str
        Stripped text or default value.
    """
    if element is None:
        return default

    node = element.find(tag)

    if node is None or node.text is None:
        return default

    value = node.text.strip()

    return value if value else default


def extract_regex_value(
    pattern: str,
    text: str,
    default: str = "N/A",
    flags: int = re.IGNORECASE,
) -> str:
    """
    Extract the first regex capture group from text.

    Parameters
    ----------
    pattern:
        Regular expression containing a capture group.
    text:
        Input text.
    default:
        Value returned when no match is found.
    flags:
        Regular-expression flags.

    Returns
    -------
    str
        Extracted value or default.
    """
    match = re.search(pattern, text, flags)

    if not match:
        return default

    value = match.group(1).strip()

    return value if value else default


def write_tsv(
    records: List[Record],
    output_tsv: str,
    fieldnames: List[str],
) -> None:
    """
    Write records to a TSV file (always overwrite).

    Parameters
    ----------
    records:
        Records to write.
    output_tsv:
        Output TSV path.
    fieldnames:
        TSV column names.
    """
    if not records:
        logger.warning("没有可写入的记录: %s", output_tsv)
        return

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as tsvfile:

        writer = csv.DictWriter(
            tsvfile,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(records)

    logger.info(
        "已写入 %d 条记录到 %s",
        len(records),
        output_tsv,
    )


def read_nonempty_lines(filepath: str) -> List[str]:
    """
    Read non-empty lines from a text file.
    """
    with open(filepath, "r", encoding="utf-8") as handle:
        return [
            line.strip()
            for line in handle
            if line.strip()
        ]


# ============================================================================
# GEO text parser
# ============================================================================

def extract_geo_info(
    filepath: str,
) -> List[Record]:
    """
    Extract GEO information from a GDS/GEO text file.

    Extracted fields
    ----------------
    Accession
    Title
    FTP_Download
    Organism

    The function only parses data and does not write output files.
    """
    results: List[Record] = []

    if not os.path.exists(filepath):
        logger.error("文件不存在: %s", filepath)
        return results

    with open(filepath, "r", encoding="utf-8") as handle:
        text = handle.read()

    # Split entries such as:
    #
    # 1. ...
    # 2. ...
    entry_blocks = re.split(
        r"^\s*\d+\.",
        text,
        flags=re.MULTILINE,
    )

    for block in entry_blocks:
        if not block.strip():
            continue

        title_match = re.match(
            r"\s*(.+)",
            block,
        )

        title = (
            title_match.group(1).strip()
            if title_match
            else "N/A"
        )

        organism = extract_regex_value(
            r"Organism:\s*([^\r\n]+)",
            block,
        )

        ftp = extract_regex_value(
            r"FTP download:\s*.*?(ftp://[^\s]+)",
            block,
        )

        accession = extract_regex_value(
            r"Series\s+Accession:\s*(\S+)",
            block,
        )

        if accession == "N/A":
            continue

        results.append(
            {
                "Accession": accession.split("\t")[0],
                "Title": title,
                "FTP_Download": ftp,
                "Organism": organism,
            }
        )

    logger.info(
        "从 %s 提取到 %d 条 GEO 记录",
        filepath,
        len(results),
    )

    return results


# ============================================================================
# GEO JSON parser
# ============================================================================

def parse_geo_summary_json(
    filepath: str,
) -> List[Record]:
    """
    Parse GEO summary JSON/text output and extract GSM records.

    Expected fields
    ---------------
    GSM
    Series
    FTP download
    Organism
    Source name

    The parser is intentionally tolerant of the text-like JSON/docsum
    representation produced by EDirect.
    """
    results: List[Record] = []

    if not os.path.exists(filepath):
        logger.error("文件不存在: %s", filepath)
        return results

    with open(filepath, "r", encoding="utf-8") as handle:
        text = handle.read()

    if not text.strip():
        logger.warning("文件为空: %s", filepath)
        return results

    # EDirect 的输出有时表现为：
    #
    # 1. ...
    # 2. ...
    #
    # 因此保留这种 block-based parsing。
    entry_blocks = re.split(
        r"^\s*\d+\.",
        text,
        flags=re.MULTILINE,
    )

    for block in entry_blocks:

        if not re.search(
            r"Accession:\s*GSM\d+",
            block,
            flags=re.IGNORECASE,
        ):
            continue

        data: Record = {}

        data["GSM"] = extract_regex_value(
            r"Accession:\s*(GSM\d+)",
            block,
        )

        data["FTP download"] = extract_regex_value(
            r"FTP download:\s*.*?(ftp://[^\s]+)",
            block,
        )

        data["Series"] = extract_regex_value(
            r"Series:.*?(GSE\d+)",
            block,
        )

        data["Organism"] = extract_regex_value(
            r"Organism:\s*([^\r\n]+)",
            block,
        )

        data["Source name"] = extract_regex_value(
            r"Source name:\s*([^\r\n]+)",
            block,
        )

        if data["GSM"].startswith("GSM"):
            results.append(data)

    logger.info(
        "从 %s 提取到 %d 条 GSM 记录",
        filepath,
        len(results),
    )

    return results


# ============================================================================
# GEO XML parser
# ============================================================================

def parse_geo_summary_xml(
    file_path: str,
) -> List[Record]:
    """
    Parse GEO summary XML and extract GSM records belonging to GSE entries.

    Extracted fields
    ----------------
    GSMId
    GSMDescription
    GSEId
    GSEDescription
    GSEDetailedDescription
    Organism
    GDSType
    PublishDate
    SamplesNum
    BioProject
    FTPLink

    The function only parses data and does not write output files.
    """
    results: List[Record] = []

    if not os.path.exists(file_path):
        logger.error("文件不存在: %s", file_path)
        return results

    try:
        tree = ET.parse(file_path)
        root = tree.getroot()

    except ET.ParseError as exc:
        logger.error(
            "XML 解析错误: %s, 文件: %s",
            exc,
            file_path,
        )
        return results

    except OSError as exc:
        logger.error(
            "无法读取 XML 文件: %s, %s",
            file_path,
            exc,
        )
        return results

    for doc_summary in root.findall(".//DocumentSummary"):

        accession = get_text_value(
            doc_summary,
            "Accession",
        )

        # 这里只处理 GSE
        if not accession.startswith("GSE"):
            continue

        gse_id = accession

        title = get_text_value(
            doc_summary,
            "title",
        )

        gds_type = get_text_value(
            doc_summary,
            "gdsType",
        )

        publish_date = get_text_value(
            doc_summary,
            "PDAT",
        )

        n_samples = get_text_value(
            doc_summary,
            "n_samples",
        )

        ftp_link = get_text_value(
            doc_summary,
            "FTPLink",
        )

        bioproject = get_text_value(
            doc_summary,
            "BioProject",
        )

        taxon = get_text_value(
            doc_summary,
            "taxon",
        )

        summary = get_text_value(
            doc_summary,
            "summary",
        )

        samples = doc_summary.findall(
            "./Samples/Sample"
        )

        for sample_element in samples:

            gsm_id = get_text_value(
                sample_element,
                "Accession",
            )

            gsm_title = get_text_value(
                sample_element,
                "Title",
            )

            # 只保留 GSM
            if not gsm_id.startswith("GSM"):
                continue

            results.append(
                {
                    "GSMId": gsm_id,
                    "GSMDescription": gsm_title,
                    "GSEId": gse_id,
                    "GSEDescription": title,
                    "GSEDetailedDescription": summary,
                    "Organism": taxon,
                    "GDSType": gds_type,
                    "PublishDate": publish_date,
                    "SamplesNum": n_samples,
                    "BioProject": bioproject,
                    "FTPLink": ftp_link,
                }
            )

    logger.info(
        "从 %s 提取到 %d 条 GSM/GSE 记录",
        file_path,
        len(results),
    )

    return results


# ============================================================================
# Batch parsing
# ============================================================================

def batch_extract_info(
    folder_path: str,
    suffix: str,
    callback: Callable[[str], List[Record]],
) -> List[Record]:
    """
    Batch parse files in a directory.

    Parameters
    ----------
    folder_path:
        Directory containing input files.
    suffix:
        File suffix to process.
    callback:
        Parser function.

    Returns
    -------
    list[dict]
        All parsed records.

    Notes
    -----
    This function deliberately does NOT write the output file.
    All records are collected first and written exactly once by the caller.
    """
    folder = Path(folder_path)

    if not folder.exists():
        logger.error("目录不存在: %s", folder_path)
        return []

    if not folder.is_dir():
        logger.error("不是目录: %s", folder_path)
        return []

    all_results: List[Record] = []

    for filepath in sorted(folder.iterdir()):

        if not filepath.is_file():
            continue

        if not filepath.name.endswith(suffix):
            continue

        logger.info(
            "正在处理文件: %s",
            filepath,
        )

        try:
            results = callback(str(filepath))
            all_results.extend(results)

        except Exception:
            logger.exception(
                "处理文件失败: %s",
                filepath,
            )

    logger.info(
        "目录 %s 共提取 %d 条记录",
        folder_path,
        len(all_results),
    )

    return all_results


# ============================================================================
# EDirect command-line implementation
# ============================================================================

def run_edirect_command(
    id: str,
    outfile: str,
    database: EDirectDatabase = "gds",
    format: EDirectFormat = "json",
) -> bool:
    """
    Execute:

        esearch -db DATABASE -query ID |
            efetch -format FORMAT > OUTFILE

    Parameters
    ----------
    id:
        GEO/SRA accession or query.
    outfile:
        Output file.
    database:
        EDirect database.
    format:
        EDirect output format.

    Returns
    -------
    bool
        True when the command succeeds and produces non-empty output.
    """
    output_path = Path(outfile)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(
        "执行 EDirect: db=%s, id=%s, format=%s",
        database,
        id,
        format,
    )

    try:
        with output_path.open(
            "w",
            encoding="utf-8",
        ) as f_out:

            p1 = subprocess.Popen(
                [
                    "esearch",
                    "-db",
                    database,
                    "-query",
                    id,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            p2 = subprocess.Popen(
                [
                    "efetch",
                    "-format",
                    format,
                ],
                stdin=p1.stdout,
                stdout=f_out,
                stderr=subprocess.PIPE,
                text=True,
            )

            if p1.stdout is not None:
                p1.stdout.close()

            _, p2_stderr = p2.communicate()
            p1_stderr = p1.stderr.read() if p1.stderr else ""

            p1_returncode = p1.wait()
            p2_returncode = p2.returncode

        if p1_returncode != 0:
            logger.error(
                "esearch 执行失败 [%s]: %s",
                p1_returncode,
                p1_stderr.strip(),
            )
            return False

        if p2_returncode != 0:
            logger.error(
                "efetch 执行失败 [%s]: %s",
                p2_returncode,
                p2_stderr.strip(),
            )
            return False

    except FileNotFoundError as exc:
        logger.error(
            "未找到 EDirect 命令，请确认 esearch/efetch 已安装: %s",
            exc,
        )
        return False

    except OSError as exc:
        logger.error(
            "执行 EDirect 失败: %s",
            exc,
        )
        return False

    if not output_path.exists():
        logger.error(
            "EDirect 未生成输出文件: %s",
            outfile,
        )
        return False

    if output_path.stat().st_size == 0:
        logger.error(
            "EDirect 输出文件为空: %s",
            outfile,
        )
        return False

    # 保留原来的随机延迟，降低连续请求频率。
    time.sleep(random.uniform(0.5, 1.5))

    logger.info(
        "EDirect 完成: %s",
        outfile,
    )

    return True


# ============================================================================
# Bio.Entrez implementation
# ============================================================================

def run_edirect_command_bio(
    id: str,
    outfile: str,
    database: EDirectDatabase = "gds",
    format: EDirectFormat = "json",
    email: Optional[str] = None,
    max_retries: int = 3,
    retry_delay: float = 1.0,
) -> bool:
    """
    Use Bio.Entrez as a Python implementation of the ESearch/EFetch workflow.

    Equivalent workflow:

        ESearch -> obtain UIDs
        EFetch  -> retrieve records

    Parameters
    ----------
    id:
        Search query.
    outfile:
        Output file.
    database:
        NCBI database.
    format:
        json / docsum / runinfo.
    email:
        Email required/recommended by NCBI.
    max_retries:
        Maximum retry count.
    retry_delay:
        Initial retry delay in seconds.

    Returns
    -------
    bool
        True when retrieval succeeds and output is non-empty.
    """
    if Entrez is None:
        logger.error(
            "BioPython 未安装，无法使用 Bio.Entrez。"
        )
        return False

    if email:
        Entrez.email = email
    elif not getattr(Entrez, "email", None):
        logger.warning(
            "未设置 Entrez.email。建议通过 --email 提供邮箱。"
        )

    output_path = Path(outfile)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    api_params_map = {
        "json": {
            "rettype": "docsum",
            "retmode": "json",
        },
        "docsum": {
            "rettype": "docsum",
            "retmode": "xml",
        },
        "runinfo": {
            "rettype": "runinfo",
            "retmode": "text",
        },
    }

    if format not in api_params_map:
        logger.error(
            "不支持的格式: %s",
            format,
        )
        return False

    params = api_params_map[format]

    logger.info(
        "使用 Bio.Entrez 获取: id=%s, db=%s, format=%s",
        id,
        database,
        format,
    )

    for attempt in range(1, max_retries + 1):

        try:
            # ------------------------------------------------------------
            # ESearch
            # ------------------------------------------------------------
            search_handle = Entrez.esearch(
                db=database,
                term=id,
            )

            try:
                search_results = Entrez.read(
                    search_handle
                )
            finally:
                search_handle.close()

            id_list = search_results.get(
                "IdList",
                [],
            )

            if not id_list:
                logger.warning(
                    "在 %s 数据库中未检索到: %s",
                    database,
                    id,
                )
                return False

            # ------------------------------------------------------------
            # EFetch
            # ------------------------------------------------------------
            fetch_handle = Entrez.efetch(
                db=database,
                id=id_list,
                rettype=params["rettype"],
                retmode=params["retmode"],
            )

            try:
                data = fetch_handle.read()
            finally:
                fetch_handle.close()

            # ------------------------------------------------------------
            # Write output
            # ------------------------------------------------------------
            if isinstance(data, bytes):
                with output_path.open(
                    "wb"
                ) as f_out:
                    f_out.write(data)
            else:
                with output_path.open(
                    "w",
                    encoding="utf-8",
                ) as f_out:
                    f_out.write(data)

            if (
                not output_path.exists()
                or output_path.stat().st_size == 0
            ):
                raise RuntimeError(
                    "NCBI 返回空数据"
                )

            logger.info(
                "Bio.Entrez 获取成功: %s",
                outfile,
            )

            # NCBI 请求之间保留间隔。
            time.sleep(0.35)

            return True

        except (HTTPError, URLError) as exc:

            logger.warning(
                "NCBI 网络请求失败，第 %d/%d 次: %s",
                attempt,
                max_retries,
                exc,
            )

        except Exception as exc:

            logger.warning(
                "NCBI 请求失败，第 %d/%d 次: %s",
                attempt,
                max_retries,
                exc,
            )

        if attempt < max_retries:

            delay = retry_delay * (2 ** (attempt - 1))

            logger.info(
                "%.1f 秒后重试...",
                delay,
            )

            time.sleep(delay)

    logger.error(
        "经过 %d 次尝试后仍然失败: %s",
        max_retries,
        id,
    )

    return False


# ============================================================================
# RunInfo parser
# ============================================================================

RUNINFO_COLUMNS = [
    "Run",
    "SampleName",
    "BioProject",
    "ScientificName",
    "Sex",
    "Disease",
    "Tumor",
    "CenterName",
    "LibraryStrategy",
    "LibraryLayout",
    "avgLength",
]


def parse_runinfo_file(
    filepath: str,
) -> List[Record]:
    """
    Parse one SRA RunInfo CSV file.

    Missing columns are filled with empty strings instead of raising KeyError.
    """
    if not os.path.exists(filepath):
        logger.error(
            "RunInfo 文件不存在: %s",
            filepath,
        )
        return []

    if os.path.getsize(filepath) == 0:
        logger.warning(
            "RunInfo 文件为空，跳过: %s",
            filepath,
        )
        return []

    try:
        df = pd.read_csv(
            filepath,
            sep=",",
            dtype=str,
            keep_default_na=False,
        )

    except Exception as exc:
        logger.error(
            "读取 RunInfo 失败: %s, error=%s",
            filepath,
            exc,
        )
        return []

    # 缺少的字段补为空字符串
    for column in RUNINFO_COLUMNS:
        if column not in df.columns:
            logger.warning(
                "%s 缺少字段 %s，使用空值填充",
                filepath,
                column,
            )
            df[column] = ""

    df = df[RUNINFO_COLUMNS]

    return df.to_dict(
        orient="records"
    )


def runinfo_formatted(
    folder_path: str,
    suffix: str = "_runinfo.csv",
) -> List[Record]:
    """
    Parse all RunInfo CSV files in a directory.

    The function only parses data and does not write output files.
    """
    folder = Path(folder_path)

    if not folder.exists():
        logger.error(
            "目录不存在: %s",
            folder_path,
        )
        return []

    all_results: List[Record] = []

    for filepath in sorted(folder.iterdir()):

        if not filepath.is_file():
            continue

        if not filepath.name.endswith(suffix):
            continue

        logger.info(
            "正在处理 RunInfo 文件: %s",
            filepath,
        )

        results = parse_runinfo_file(
            str(filepath)
        )

        all_results.extend(results)

    logger.info(
        "RunInfo 共提取 %d 条记录",
        len(all_results),
    )

    return all_results


# ============================================================================
# Duplicate removal
# ============================================================================

def remove_duplicate(
    filepath: str,
    output_tsv: str,
) -> None:
    """
    Remove duplicate rows according to the first column.
    """
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        keep_default_na=False,
    )

    if df.empty:
        logger.warning(
            "输入文件为空: %s",
            filepath,
        )
        df.to_csv(
            output_tsv,
            sep="\t",
            index=False,
            header=False,
        )
        return

    df = df.drop_duplicates(
        subset=df.columns[0],
        keep="first",
    )

    df.to_csv(
        output_tsv,
        sep="\t",
        index=False,
        header=False,
    )

    logger.info(
        "去重完成: %s -> %s",
        filepath,
        output_tsv,
    )


# ============================================================================
# GSE batch download + parse
# ============================================================================

def collect_gse_runinfo(
    gse_ids: List[str],
    output_tsv: str,
    download_dir: Optional[str] = None,
    database: EDirectDatabase = "gds",
    format: EDirectFormat = "json",
    email: Optional[str] = None,
    use_bio: bool = True,
) -> List[Record]:
    """
    Download and parse information for a list of GSE IDs.

    Workflow
    --------
    1. Download all GSE records.
    2. Parse all downloaded files.
    3. Aggregate records.
    4. Write the final TSV exactly once.

    This avoids the original implementation's repeated parsing/writing bug.
    """
    if not gse_ids:
        logger.warning("GSE ID 列表为空")
        return []

    if download_dir is None:
        output_path = Path(output_tsv)

        download_dir = str(
            output_path.parent
            / f"{output_path.stem}_download"
        )

    download_path = Path(download_dir)
    download_path.mkdir(
        parents=True,
        exist_ok=True,
    )

    logger.info(
        "GSE 中间文件目录: %s",
        download_path,
    )

    successful_files: List[str] = []

    # ------------------------------------------------------------------------
    # Step 1: download
    # ------------------------------------------------------------------------

    for gse in gse_ids:

        gse = gse.strip()

        if not gse:
            continue

        if not re.fullmatch(
            r"GSE\d+",
            gse,
            flags=re.IGNORECASE,
        ):
            logger.warning(
                "不是标准 GSE ID，仍尝试查询: %s",
                gse,
            )

        if format == "runinfo":
            extension = ".csv"
        elif format == "json":
            extension = ".json"
        elif format == "docsum":
            extension = ".xml"
        else:
            logger.error(
                "不支持的格式: %s",
                format,
            )
            continue

        outfile = download_path / (
            f"{gse}_runinfo{extension}"
        )

        if use_bio:
            success = run_edirect_command_bio(
                gse,
                str(outfile),
                database=database,
                format=format,
                email=email,
            )
        else:
            success = run_edirect_command(
                gse,
                str(outfile),
                database=database,
                format=format,
            )

        if success:
            successful_files.append(
                str(outfile)
            )

    if not successful_files:
        logger.warning(
            "没有成功下载任何 GSE 数据"
        )
        return []

    # ------------------------------------------------------------------------
    # Step 2: parse all downloaded files ONCE
    # ------------------------------------------------------------------------

    if format == "json":

        all_results: List[Record] = []

        for filepath in successful_files:
            all_results.extend(
                parse_geo_summary_json(filepath)
            )

        fieldnames = [
            "GSM",
            "Series",
            "FTP download",
            "Organism",
            "Source name",
        ]

    elif format == "docsum":

        all_results = []

        for filepath in successful_files:
            all_results.extend(
                parse_geo_summary_xml(filepath)
            )

        fieldnames = [
            "GSMId",
            "GSMDescription",
            "GSEId",
            "GSEDescription",
            "GSEDetailedDescription",
            "Organism",
            "GDSType",
            "PublishDate",
            "SamplesNum",
            "BioProject",
            "FTPLink",
        ]

    elif format == "runinfo":

        all_results = []

        for filepath in successful_files:
            all_results.extend(
                parse_runinfo_file(filepath)
            )

        fieldnames = RUNINFO_COLUMNS

    else:
        logger.error(
            "不支持的格式: %s",
            format,
        )
        return []

    # ------------------------------------------------------------------------
    # Step 3: optional de-duplication
    # ------------------------------------------------------------------------

    # 按不同模式使用对应的唯一标识。
    unique_key = None

    if format == "json":
        unique_key = "GSM"

    elif format == "docsum":
        unique_key = "GSMId"

    elif format == "runinfo":
        unique_key = "Run"

    if unique_key:
        seen = set()
        deduplicated_results: List[Record] = []

        for record in all_results:

            key = record.get(
                unique_key,
                "",
            )

            if not key:
                deduplicated_results.append(
                    record
                )
                continue

            if key in seen:
                continue

            seen.add(key)
            deduplicated_results.append(
                record
            )

        logger.info(
            "去重: %d -> %d 条记录",
            len(all_results),
            len(deduplicated_results),
        )

        all_results = deduplicated_results

    # ------------------------------------------------------------------------
    # Step 4: write exactly once
    # ------------------------------------------------------------------------

    write_tsv(
        all_results,
        output_tsv,
        fieldnames,
    )

    return all_results


# ============================================================================
# CLI
# ============================================================================

def parser_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Extract information from GEO GSE/GSM "
            "and SRA RunInfo data."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help=(
            "Input file or directory. "
            "For gse_txt mode, this is a text file containing GSE IDs."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help=(
            "Final output TSV file. "
            "For all modes this argument represents a file, not a directory."
        ),
    )

    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=[
            "json",
            "xml",
            "runinfo",
            "gds",
            "gse_txt",
        ],
        required=True,
        help=(
            "Processing mode: "
            "json, xml, runinfo, gds, or gse_txt."
        ),
    )

    parser.add_argument(
        "--download-dir",
        type=str,
        default=None,
        help=(
            "Directory for intermediate files generated "
            "by gse_txt mode. "
            "Default: <output_stem>_download/"
        ),
    )

    parser.add_argument(
        "--database",
        type=str,
        choices=["gds", "sra"],
        default="gds",
        help="NCBI database used by gse_txt mode.",
    )

    parser.add_argument(
        "--format",
        type=str,
        choices=["json", "docsum", "runinfo"],
        default="json",
        help=(
            "NCBI retrieval format used by gse_txt mode. "
            "Default: json."
        ),
    )

    parser.add_argument(
        "--email",
        type=str,
        default="2530320102@qq.com",
        help=(
            "Email used by Bio.Entrez for NCBI requests."
        ),
    )

    parser.add_argument(
        "--edirect",
        action="store_true",
        help=(
            "Use external esearch/efetch instead of "
            "Bio.Entrez in gse_txt mode."
        ),
    )

    return parser.parse_args()


# ============================================================================
# Main
# ============================================================================

def main() -> None:
    """
    Main CLI entry point.
    """
    args = parser_args()

    input_path = args.input
    output_tsv = args.output
    mode = args.mode

    # ------------------------------------------------------------------------
    # JSON
    # ------------------------------------------------------------------------

    if mode == "json":

        results = batch_extract_info(
            input_path,
            suffix=".json",
            callback=parse_geo_summary_json,
        )

        fieldnames = [
            "GSM",
            "Series",
            "FTP download",
            "Organism",
            "Source name",
        ]

        write_tsv(
            results,
            output_tsv,
            fieldnames,
        )

    # ------------------------------------------------------------------------
    # XML
    # ------------------------------------------------------------------------

    elif mode == "xml":

        results = batch_extract_info(
            input_path,
            suffix=".xml",
            callback=parse_geo_summary_xml,
        )

        fieldnames = [
            "GSMId",
            "GSMDescription",
            "GSEId",
            "GSEDescription",
            "GSEDetailedDescription",
            "Organism",
            "GDSType",
            "PublishDate",
            "SamplesNum",
            "BioProject",
            "FTPLink",
        ]

        write_tsv(
            results,
            output_tsv,
            fieldnames,
        )


    elif mode == "runinfo":

        results = runinfo_formatted(
            input_path,
            suffix="_runinfo.csv",
        )

        write_tsv(
            results,
            output_tsv,
            RUNINFO_COLUMNS,
        )

    elif mode == "gds":

        results = extract_geo_info(
            input_path
        )

        fieldnames = [
            "Accession",
            "Title",
            "FTP_Download",
            "Organism",
        ]

        write_tsv(
            results,
            output_tsv,
            fieldnames,
        )


    elif mode == "gse_txt":

        gse_ids = read_nonempty_lines(
            input_path
        )

        collect_gse_runinfo(
            gse_ids=gse_ids,
            output_tsv=output_tsv,
            download_dir=args.download_dir,
            database=args.database,
            format=args.format,
            email=args.email,
            use_bio=not args.edirect,
        )

    else:
        raise ValueError(
            f"Unsupported mode: {mode}"
        )


if __name__ == "__main__":
    main()

