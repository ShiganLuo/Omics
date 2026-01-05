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
awk '{print $2}' samples.tsv | \
python geo_gsm_cli.py \
  --mode both \
  --html-dir html \
  --out-csv out.csv
"""

import argparse
import csv
import glob
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Set

import requests
from bs4 import BeautifulSoup


# ============================================================
# GSM Resolver
# ============================================================
class GSMResolver:
    def __init__(
        self,
        gsm: Optional[List[str]] = None,
        gsm_file: Optional[str] = None,
        html_dir: Optional[str] = None,
        id_column: str = "GSM"
    ):
        self.gsm = gsm
        self.gsm_file = gsm_file
        self.html_dir = html_dir
        self.id_column = id_column.lower()

    def resolve(self) -> List[str]:
        gsms: Set[str] = set()

        if self.gsm:
            gsms.update(self.gsm)

        if self.gsm_file:
            gsms.update(self._from_file(self.gsm_file))

        if self.html_dir:
            gsms.update(self._from_html_dir(self.html_dir))

        if not gsms and not sys.stdin.isatty():
            gsms.update(self._from_stdin())

        if not gsms:
            raise ValueError("No GSM found from any source")

        return sorted(gsms)

    def _from_file(self, path: str) -> Set[str]:
        gsms: Set[str] = set()
        ext = os.path.splitext(path)[1].lower()

        if ext in (".csv", ".tsv"):
            delimiter = "," if ext == ".csv" else "\t"
            with open(path, encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter=delimiter)
                col = self._find_id_column(reader.fieldnames)
                if not col:
                    raise ValueError(f"No column '{self.id_column}' in {path}")
                for row in reader:
                    val = row.get(col)
                    if val:
                        gsms.add(val.strip())
        else:
            with open(path, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        gsms.add(line)

        return gsms

    def _from_html_dir(self, html_dir: str) -> Set[str]:
        return {
            fn.replace(".html", "")
            for fn in os.listdir(html_dir)
            if fn.endswith(".html")
        }

    def _from_stdin(self) -> Set[str]:
        return {line.strip() for line in sys.stdin if line.strip()}

    def _find_id_column(self, columns):
        for c in columns or []:
            if c.lower() == self.id_column:
                return c
        return None


# ============================================================
# Downloader (concurrent + retry + cache)
# ============================================================
def download_one(
    gsm: str,
    outdir: str,
    retries: int,
    force: bool
) -> bool:
    out_html = os.path.join(outdir, f"{gsm}.html")
    if os.path.exists(out_html) and not force:
        print(f"[SKIP] {gsm} (cached)")
        return True

    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"

    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            with open(out_html, "w", encoding="utf-8") as f:
                f.write(r.text)
            print(f"[OK] {gsm}")
            return True
        except Exception as e:
            print(f"[RETRY {attempt}/{retries}] {gsm}: {e}")
            time.sleep(2)

    print(f"[FAIL] {gsm}")
    return False


def download_batch(
    gsm_list: List[str],
    html_dir: str,
    workers: int,
    retries: int,
    force: bool
):
    os.makedirs(html_dir, exist_ok=True)

    with ThreadPoolExecutor(max_workers=workers) as exe:
        futures = {
            exe.submit(download_one, gsm, html_dir, retries, force): gsm
            for gsm in gsm_list
        }
        for _ in as_completed(futures):
            pass


# ============================================================
# Extract Characteristics
# ============================================================
def extract_characteristics(html_path: str) -> Dict[str, str]:
    with open(html_path, encoding="utf-8") as f:
        soup = BeautifulSoup(f.read(), "html.parser")

    result: Dict[str, str] = {}

    for tr in soup.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) >= 2 and tds[0].get_text(strip=True) == "Characteristics":
            for line in tds[1].get_text("\n", strip=True).split("\n"):
                if ":" in line:
                    k, v = line.split(":", 1)
                    result[k.strip()] = v.strip()

    return result


def extract_batch(html_dir: str, outfile: str):
    records: Dict[str, Dict[str, str]] = {}

    for html in glob.glob(os.path.join(html_dir, "*.html")):
        gsm = os.path.basename(html).replace(".html", "")
        records[gsm] = extract_characteristics(html)

    all_keys = sorted({k for v in records.values() for k in v})

    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["GSM"] + all_keys)
        writer.writeheader()
        for gsm, info in records.items():
            row = {"GSM": gsm}
            row.update(info)
            writer.writerow(row)

    print(f"[DONE] Extracted → {outfile}")



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

    parser.add_argument("--html-dir", required=True)
    parser.add_argument("--out-csv")

    parser.add_argument("--workers", type=int, default=4, help="Concurrent downloads")
    parser.add_argument("--retries", type=int, default=3, help="Retry times")
    parser.add_argument("--force", action="store_true", help="Force re-download")

    args = parser.parse_args()

    resolver = GSMResolver(
        gsm=args.gsm,
        gsm_file=args.gsm_file,
        html_dir=args.html_dir if args.from_html_dir else None,
        id_column=args.id_column
    )

    gsm_list = resolver.resolve()

    if args.mode in ("download", "both"):
        download_batch(
            gsm_list,
            args.html_dir,
            args.workers,
            args.retries,
            args.force
        )

    if args.mode in ("extract", "both"):
        if not args.out_csv:
            raise ValueError("--out-csv required for extract/both")
        extract_batch(args.html_dir, args.out_csv)


if __name__ == "__main__":
    main()
