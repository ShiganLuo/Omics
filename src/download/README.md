# Description

all is done

## ftp_download.py

Parameterized FTP batch downloader (refactored from the old hard-coded
`CNP_download` / `CNGB_download` entry points). Collects FTP URLs from one of
three input modes, then downloads them in parallel with retries, atomic
writes (`*.part` then rename) and skip-existing support.

Logging uses `src/common/util/LogUtil.py`:
`<time> | <level> | FTPDownloader | <message>`, written to stdout and, when
`-l/--log` is given, additionally to that log file.

### Input modes

| mode     | input                                            | URL selection                                              |
|----------|--------------------------------------------------|------------------------------------------------------------|
| `urls`   | `-f/--url-file` text file, one record per line    | column `-c/--url-column` (0-based, default 1)              |
| `table`  | `-t/--table` delimited table with header          | FTP columns detected from cell content (`ftp://` URLs and HTTPS FTP-archive links such as `https://ftp.pride.ebi.ac.uk/...`) by default; `-p/--url-pattern` forces columns whose name contains the given string (repeatable) |
| `metadata` | `-m/--meta` experiment table + `-R/--ftp-ref` link list | rows filtered by `-C/--condition-col` / `-V/--condition-value` (repeatable), run accessions matched against `--ftp-ref` |

### Common options

- `-o/--outdir` (required), `-l/--log` (optional; logs go to stdout without it)
- `-j/--jobs` (default 4), `-r/--retries` (default 3), `-T/--timeout` (default 60s)
- `-u/--username`, `-P/--password` (anonymous login when omitted)
- `-w/--overwrite` re-download existing files (default: skip)
- `-n/--dry-run` only log planned transfers

### Outputs

- downloaded files in `--outdir`
- `metadata` mode also writes `matched_urls.txt` and `unmatched_accessions.txt`
- failed URLs are written to `download_failed.txt`; exit code 1 on any failure

### Examples

```bash
# 1. plain URL list
python ftp_download.py urls --url-file links.txt \
    -o fastq/ -l log/ftp_download.log

# 2. table: FTP columns auto-detected from cell content
python ftp_download.py table -t run.tsv \
    -o fastq/ -l log/ftp_download.log --jobs 5

# or force columns by name (e.g. the DownLoad columns of a CNGB run table)
python ftp_download.py table -t run.tsv -p DownLoad \
    -o fastq/ -l log/ftp_download.log --jobs 5

# 3. metadata-driven selection (repeat --condition-value for several groups)
python ftp_download.py metadata \
    --meta metadata_CNP0003135_experiment.tsv \
    --ftp-ref data_download_links_CNP0003135_ftp.txt \
    --condition-col experiment_title \
    --condition-value "B2 CoD2 Human p65KO" \
    --condition-value "B2 CoD3 Human p65KO" \
    -o fastq/ -l log/ftp_download.log --jobs 4
```

The old hard-coded workflows map to: CNP_download -> `metadata` example 3
(experiment metadata + FTP link list, condition filter on `experiment_title`),
CNGB_download -> `table` example 2 (`run.tsv`, `DownLoad` columns).

Note: substring condition matching is literal (no regex); use `--exact-match`
for whole-cell equality.

## 开发日志

ascp.sh 加上 sleep 10s;10s对于下载时间的影响忽略不计，确能避免ascp连续下载突然降速或者下载

2025.04.25:
    - 增加gzipTest.sh 判断ascp下载文件是否压缩正确
    - 增加diff.py 判断那些文件没有下载

## gloubs

ENA collections: 网页搜索EMBL-EBI PUBLIC DATA，在collections找到其详情页，最下面有其UUID
zhang_c2: 76fbc1da-94c8-11f1-9c2b-02ce27bde401