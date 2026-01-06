#!/usr/bin/env python3
"""
GEO GSM Downloader & Extractor
==============================

This script supports:
- Resolving GSM IDs from flexible sources
- Concurrent downloading with retry & cache
- Extracting "Characteristics" from GEO GSM HTML pages

-------------------------
Supported GSM sources
-------------------------
Priority order:
1. --gsm GSM1 GSM2 ...
2. --gsm-file (txt / csv / tsv, as long as it has the ID column)
3. --from-html-dir (infer from existing HTML files)
4. stdin (pipe)

-------------------------
Run modes
-------------------------
--mode download   : only download HTML
--mode extract    : only extract from existing HTML
--mode both       : download then extract

-------------------------
Examples
-------------------------

1. Download + extract from a sample sheet
-----------------------------------------
python geo_gsm_cli.py \
  --mode both \
  --gsm-file samples.csv \
  --id-column geo_accession \
  --html-dir html \
  --out-csv Characteristics.csv \
  --workers 8 \
  --retries 3

2. Resume extraction from existing HTML
---------------------------------------
python geo_gsm_cli.py \
  --mode extract \
  --from-html-dir \
  --html-dir html \
  --out-csv Characteristics.csv

3. Force re-download
--------------------
python geo_gsm_cli.py \
  --mode download \
  --gsm GSM2786643 \
  --html-dir html \
  --force

4. HPC / pipe usage
-------------------
awk '{logging.info $2}' samples.tsv | \
python geo_gsm_cli.py \
  --mode both \
  --html-dir html \
  --out-csv out.csv
"""

import argparse
import glob
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Literal, Optional
import requests
from bs4 import BeautifulSoup
from GSMresolver import GSMResolver
import logging
import sys
import pandas as pd
from pathlib import Path
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,  # 指定输出到 stdout 而不是 stderr
	datefmt='%Y-%m-%d %H:%M:%S'
)



# ============================================================
# GSM Resolver
# ============================================================


# ============================================================
# Downloader (concurrent + retry + cache)
# ============================================================

def build_ncbi_url(acc: str) -> str:
    """
    Construct the appropriate NCBI webpage URL based on the accession prefix.

    This function determines which NCBI database page should be accessed
    according to the accession type and returns the corresponding URL.

    Supported accession types:
    - GEO accessions: GSM, GSE, GPL
      → https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
    - SRA accessions: SRX, SRR, SRS
      → https://www.ncbi.nlm.nih.gov/sra
    - BioSample accessions: SAMN
      → https://www.ncbi.nlm.nih.gov/biosample

    Parameters
    ----------
    acc : str
        NCBI accession identifier (e.g. "GSM123456", "SRX25503098", "SAMN12345678").

    Returns
    -------
    str
        A fully qualified URL pointing to the corresponding NCBI webpage.

    Raises
    ------
    ValueError
        If the accession prefix is not recognized or supported.
    """
    if acc.startswith(("GSM", "GSE", "GPL")):
        return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"

    if acc.startswith(("SRX", "SRR", "SRS")):
        return f"https://www.ncbi.nlm.nih.gov/sra?term={acc}"

    if acc.startswith("SAMN"):
        return f"https://www.ncbi.nlm.nih.gov/biosample/{acc}"

    raise ValueError(f"Unsupported accession type: {acc}")


def download_one(
    acc: str,
    outdir: str,
    retries: int,
    force: bool
) -> bool:
    """
    Download a single NCBI HTML page corresponding to a given accession.

    The function automatically determines the correct NCBI URL based on
    the accession type (GEO, SRA, or BioSample) and saves the retrieved
    HTML content to the specified output directory.

    If the HTML file already exists, the download will be skipped unless
    `force` is set to True.

    Parameters
    ----------
    acc : str
        NCBI accession identifier (e.g. GSM, GSE, SRX, SAMN).
    outdir : str
        Directory where the downloaded HTML file will be stored.
    retries : int
        Number of retry attempts if the download fails.
    force : bool
        Whether to force re-download even if the file already exists.

    Returns
    -------
    bool
        True if the HTML page is successfully downloaded or already exists;
        False if all retry attempts fail or the accession type is unsupported.
    """
    out_html = os.path.join(outdir, f"{acc}.html")

    if os.path.exists(out_html) and not force:
        logging.info(f"[HTML DOWNLOAD SKIP] {acc} (cached)")
        return True

    try:
        url = build_ncbi_url(acc)
    except ValueError as e:
        logging.error(e)
        return False

    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()

            with open(out_html, "w", encoding="utf-8") as f:
                f.write(r.text)

            logging.info(f"[OK] {acc}")
            return True

        except Exception as e:
            logging.warning(f"[RETRY {attempt}/{retries}] {acc}: {e}")
            time.sleep(2)

    logging.error(f"[FAIL] {acc}")
    return False


def download_batch(
    acc_list: List[str],
    html_dir: str,
    workers: int,
    retries: int,
    force: bool
):
    """
    Batch-download NCBI HTML pages for a list of accessions using multithreading.

    This function concurrently downloads HTML pages for multiple NCBI
    accessions (GEO, SRA, BioSample) and stores them in the specified
    output directory.

    The downloading process is parallelized using a thread pool to
    improve performance when handling large accession lists.

    Parameters
    ----------
    acc_list : List[str]
        A list of NCBI accession identifiers to download.
    html_dir : str
        Directory used to store downloaded HTML files.
    workers : int
        Number of worker threads to use for concurrent downloads.
    retries : int
        Number of retry attempts per accession if download fails.
    force : bool
        Whether to force re-download of existing HTML files.

    Returns
    -------
    None
        This function does not return a value. Download status is reported
        via logging.
    """
    os.makedirs(html_dir, exist_ok=True)

    with ThreadPoolExecutor(max_workers=workers) as exe:
        futures = {
            exe.submit(download_one, acc, html_dir, retries, force): acc
            for acc in acc_list
        }

        for _ in as_completed(futures):
            pass


# ============================================================
# Extract Characteristics
# ============================================================

def extract_feature(html_path: str, item: str = "Characteristics") -> Dict[str, str]:
    """
    Extract key–value metadata from a specified section of a GEO GSM HTML page.

    This function scans all table rows (<tr>) in the HTML document and locates
    the row whose first <td> exactly matches the given `item` name
    (e.g. "Characteristics", "SRA").

    The content of the second <td> is then split by line breaks, and each line
    is parsed as a key–value pair using the first ":" as a separator.

    Example HTML structure expected:
        <tr>
            <td>Characteristics</td>
            <td>
                tissue: liver
                strain: C57BL/6
            </td>
        </tr>

    Parameters
    ----------
    html_path : str
        Path to the GSM HTML file.
    item : str, optional
        The section name to extract from the first table column.
        Default is "Characteristics".

    Returns
    -------
    Dict[str, str]
        A dictionary mapping feature names to their corresponding values.
        If the specified section is not found, an empty dictionary is returned.
    """
    with open(html_path, encoding="utf-8") as f:
        soup = BeautifulSoup(f.read(), "html.parser")

    result: Dict[str, str] = {}

    for tr in soup.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) >= 2 and tds[0].get_text(strip=True) == item:
            for line in tds[1].get_text("\n", strip=True).split("\n"):
                if ":" in line:
                    k, v = line.split(":", 1)
                    result[k.strip()] = v.strip()

    return result


def extract_gse_ids(html_path: str) -> List[str]:
    """
    Extract GSE (Series) accession IDs from a GEO GSM HTML page.

    Parsing logic:
    1. Locate table rows (<tr>) where the first <td> starts with "Series";
    2. In the corresponding second <td>, find all <a> tags;
    3. Collect link texts that start with "GSE".

    Parameters
    ----------
    html_path : str
        Path to the GSM HTML file.

    Returns
    -------
    List[str]
        A list of GSE accession IDs (e.g. ["GSE273331", "GSE273338"]).
        Returns an empty list if no Series information is found.
    """
    with open(html_path, encoding="utf-8") as f:
        soup = BeautifulSoup(f.read(), "html.parser")

    gse_ids: List[str] = []

    for tr in soup.find_all("tr"):
        tds = tr.find_all("td")
        if not tds:
            continue

        header_text = tds[0].get_text(strip=True)

        # Match rows like "Series (2)"
        if header_text.startswith("Series"):
            for a in tds[1].find_all("a"):
                text = a.get_text(strip=True)
                if text.startswith("GSE"):
                    gse_ids.append(text)

    return gse_ids

def extract_relations_txt(
    html_path: str,
    item: Literal["SRA", "BioSample"] = "SRA"
) -> List[str]:
    """
    Extract relation accession IDs (e.g. SRA or BioSample) from a GEO GSM HTML page.

    Parsing logic:
    1. Scan all table rows (<tr>) in the HTML document;
    2. Locate the row whose first <td> exactly matches the given item
       (e.g. "SRA" or "BioSample");
    3. Extract all <a> tag texts from the second <td>;
    4. Filter accession IDs by prefix:
       - "SRX" for SRA
       - "SAMN" for BioSample

    Parameters
    ----------
    html_path : str
        Path to the GSM HTML file.
    item : {"SRA", "BioSample"}, optional
        Relation type to extract. Default is "SRA".

    Returns
    -------
    List[str]
        A list of accession IDs.
        Returns an empty list if the specified relation is not found.
    """
    with open(html_path, encoding="utf-8") as f:
        soup = BeautifulSoup(f.read(), "html.parser")

    prefix_map = {
        "SRA": "SRX",
        "BioSample": "SAMN",
    }

    prefix = prefix_map[item]
    relation_ids: List[str] = []

    for tr in soup.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 2:
            continue

        if tds[0].get_text(strip=True) != item:
            continue

        for a in tds[1].find_all("a"):
            text = a.get_text(strip=True)
            if text.startswith(prefix):
                relation_ids.append(text)

        break  # only one matching row exists

    return relation_ids

def extract_gsm_batch(html_dir: str, outfile: str | None = None) -> pd.DataFrame:
    """
    Batch-extract GSM metadata from GEO HTML files and return as a DataFrame.

    For each GSM HTML file:
    - Extract 'Characteristics' key–value pairs
    - Extract associated GSE IDs
    - Extract SRA and BioSample IDs
    - Merge them into a single row

    If multiple IDs exist, they are joined by commas.

    Parameters
    ----------
    html_dir : str
        Directory containing GSM HTML files.
    outfile : str | None, optional
        Output CSV file path. If None, no file will be written.

    Returns
    -------
    pandas.DataFrame
        A DataFrame where each row corresponds to one GSM record.
    """
    rows: list[dict] = []

    for html in glob.glob(os.path.join(html_dir, "GSM*.html")):
        gsm = os.path.basename(html).replace(".html", "")

        row: Dict[str, str] = {"GSM": gsm}

        # Characteristics
        row.update(extract_feature(html, item="Characteristics"))

        # GSE
        gse_ids = extract_gse_ids(html)
        logging.info(f"{gsm} -> GSE IDs: {gse_ids}")
        if gse_ids:
            row["GSE"] = ",".join(gse_ids)

        # SRA
        sra_ids = extract_relations_txt(html, item="SRA")
        logging.info(f"{gsm} -> SRA IDs: {sra_ids}")
        if sra_ids:
            row["SRA"] = ",".join(sra_ids)

        # BioSample
        biosample_ids = extract_relations_txt(html, item="BioSample")
        logging.info(f"{gsm} -> BioSample IDs: {biosample_ids}")
        if biosample_ids:
            row["BioSample"] = ",".join(biosample_ids)

        rows.append(row)

    df = pd.DataFrame(rows)

    # Optional CSV export
    if outfile:
        sep = "\t" if outfile.endswith(".tsv") else ","
        df.to_csv(outfile, index=False, encoding="utf-8",sep=sep)
        logging.info(f"[DONE] Extracted → {outfile}")

    return df



def extract_sra_library_and_runs(html_path: str) -> Dict[str, object]:
    """
    Extract GSM, Library metadata and Run (SRR) information
    from an NCBI SRA Experiment (SRX) HTML page.

    Priority rules
    --------------
    - GSM is extracted from: Library -> Name
    - Library fields come from the structured 'Library' block
    - Runs are extracted from the SRR table

    Returns
    -------
    {
        "SRX": "SRX25503118",
        "GSM": "GSM8426620",
        "library": {...},
        "runs": [...]
    }
    """
    srx = html_path.split("/")[-1].replace(".html", "")

    with open(html_path, encoding="utf-8") as f:
        soup = BeautifulSoup(f.read(), "html.parser")

    result: Dict[str, object] = {
        "SRX": srx,
        "GSM": None,
        "library": {},
        "runs": []
    }

    # ---------- Library ----------
    for div in soup.find_all("div", class_="expand showed sra-full-data"):
        if div.get_text(strip=True).startswith("Library"):
            for sub in div.find_all("div"):
                text = sub.get_text(strip=True)
                if ":" in text:
                    key, value = text.split(":", 1)
                    key = key.strip()
                    value = value.strip()
                    result["library"][key] = value

                    # GSM is stored as Library -> Name
                    if key == "Name" and value.startswith("GSM"):
                        result["GSM"] = value
            break

    # ---------- Runs ----------
    table = soup.find("table")
    if table:
        rows = table.find_all("tr")[1:]  # skip header
        for tr in rows:
            tds = tr.find_all("td")
            if len(tds) >= 5:
                srr = tds[0].get_text(strip=True)
                if not srr.startswith("SRR"):
                    continue

                result["runs"].append({
                    "SRR": srr,
                    "spots": tds[1].get_text(strip=True),
                    "bases": tds[2].get_text(strip=True),
                    "size": tds[3].get_text(strip=True),
                    "published": tds[4].get_text(strip=True),
                })

    return result


def extract_sra_batch(
    html_dir: str,
    outfile: Optional[str] = None
) -> pd.DataFrame:
    """
    Batch-parse SRX HTML files and extract SRX–GSM–SRR relationships.

    Each row corresponds to one SRR.
    """
    records: List[Dict[str, str]] = []

    for html in glob.glob(os.path.join(html_dir, "SRX*.html")):
        parsed = extract_sra_library_and_runs(html)

        srx = parsed["SRX"]
        gsm = parsed["GSM"]
        library = parsed["library"]
        runs = parsed["runs"]

        if not runs:
            row = {"SRX": srx, "GSM": gsm}
            row.update(library)
            records.append(row)
            continue

        for run in runs:
            row = {
                "SRX": srx,
                "GSM": gsm,
            }
            row.update(library)
            row.update(run)
            records.append(row)

    df = pd.DataFrame(records)

    # Optional CSV export
    if outfile:
        sep = "\t" if outfile.endswith(".tsv") else ","
        df.to_csv(outfile, index=False, encoding="utf-8",sep=sep)
        logging.info(f"[DONE] Extracted → {outfile}")

    return df
# ============================================================
# CLI
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Concurrent GEO GSM downloader and extractor"
    )

    parser.add_argument("--mode", choices=["download", "extract", "both"], required=True)
    parser.add_argument("--gsm", nargs="+")
    parser.add_argument("--gsm-file")
    parser.add_argument("--id-column", default="GSM")
    parser.add_argument("--from-html-dir", action="store_true")

    parser.add_argument("--outdir", required=True)

    parser.add_argument("--workers", type=int, default=4, help="Concurrent downloads")
    parser.add_argument("--retries", type=int, default=3, help="Retry times")
    parser.add_argument("--force", action="store_true", help="Force re-download")

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    html_dir = outdir / "html"
    html_dir.mkdir(parents=True, exist_ok=True)

    resolver = GSMResolver(
        gsm=args.gsm,
        gsm_file=args.gsm_file,
        html_dir=str(html_dir) if args.from_html_dir else None,
        id_column=args.id_column
    )

    gsm_list = resolver.resolve()

    if args.mode in ("download", "both"):
        download_batch(
            gsm_list,
            str(html_dir),
            args.workers,
            args.retries,
            args.force
        )

    if args.mode in ("extract", "both"):
        gsm_outfile = outdir / "gsm_metadata.csv"
        gsm_meta = extract_gsm_batch(str(html_dir), str(gsm_outfile))
        srx_list = gsm_meta["SRA"].dropna().str.split(",").explode().unique().tolist()
        download_batch(
            srx_list,
            str(html_dir),
            args.workers,
            args.retries,
            args.force
        )
        sra_outfile = outdir / "sra_metadata.csv"
        sra_meta = extract_sra_batch(str(html_dir), str(sra_outfile))
        


if __name__ == "__main__":
    main()
