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

## SRS_metadata.py

Fetch sample- and run-level metadata for SRA sample accessions (SRS...)
from the ENA API into consolidated TSV tables. Starting point when the
batch is defined by SRS accessions instead of a GSE (the `GSM_metadata.py`
/ `GEO_soft.py` flows start from GEO).

### Data sources

| source | URL | content |
|--------|-----|---------|
| ENA Browser API (sample XML) | `https://www.ebi.ac.uk/ena/browser/api/xml/<SRS>` | BioSample accession, sample alias/title, organism, all SAMPLE_ATTRIBUTEs (sex, tissue, age, ...) |
| ENA Portal API (filereport) | `.../portal/api/filereport?accession=<SRS>&result=read_run` | one row per run: study/experiment/run accessions, library info, read counts, fastq_ftp + fastq_md5 + fastq_bytes |

### Input sources

1. `-s/--srs ACC1 ACC2 ...` explicit list (repeatable)
2. `-i/--input` file: one accession per line, or a delimited table with
   ID column `-c/--id-column` (default `SRS`)
3. `-g/--group LABEL=ACC1,ACC2,...` — accessions named in group specs are
   collected too (repeatable)
4. stdin (pipe)

### Group labels

`-g/--group LABEL=ACC1,ACC2` (repeatable) and/or `-G/--group-file` table
with `--group-id-col`/`--group-name-col` columns. A `group` column is
always written into both tables; unassigned samples get an empty cell and
a warning. Labels are free text (e.g. `experiment`, `control`).

### Options

- `-o/--outdir` (required), `-l/--log` (optional; stdout always)
- `-j/--jobs` (default 4), `-r/--retries` (default 3), `-t/--timeout` (default 60s)
- `-w/--overwrite` re-fetch cached responses (default: reuse `cache/`)
- `-n/--dry-run` log the plan, no network
- `--insecure` disable TLS verification

### Outputs

- `sample_metadata.tsv` — one row per SRS: `srs, group, biosample,
  sample_alias, sample_title, organism, tax_id` + one column per sample
  attribute
- `run_metadata.tsv` — one row per run: `srs, group, run_accession, ...,
  library_*, read_count, base_count, fastq_*`
- `cache/<SRS>.xml`, `cache/<SRS>_runs.tsv` — raw ENA responses
- `srs_failed.txt` — failed accessions (+error), written only on failure;
  exit code 1 on any failure, 2 on invalid (non-SRS) input

### Examples

```bash
# batch with group labels
python SRS_metadata.py \
    -g experiment=SRS18518057,SRS18518052,SRS18518053 \
    -g control=SRS18527637,SRS18527634 \
    -o /data/16T/Data/RNAseq_Ovary_CD38_PRJNA1002222 \
    --jobs 4 --log /data/16T/Data/RNAseq_Ovary_CD38_PRJNA1002222/srs_metadata.log

# accessions from a file, groups from a table
python SRS_metadata.py -i srs_list.txt --group-file groups.tsv -o results/
```

Note: the ENA filereport endpoint hoists `run_accession` to the first
column regardless of the requested field order; the run table validates
columns as a set and keeps the API's returned column order.

## 开发日志

ascp.sh 加上 sleep 10s;10s对于下载时间的影响忽略不计，确能避免ascp连续下载突然降速或者下载

2025.04.25:
    - 增加gzipTest.sh 判断ascp下载文件是否压缩正确
    - 增加diff.py 判断那些文件没有下载

## gloubs

ENA collections: 网页搜索EMBL-EBI PUBLIC DATA，在collections找到其详情页，最下面有其UUID
zhang_c2: 76fbc1da-94c8-11f1-9c2b-02ce27bde401