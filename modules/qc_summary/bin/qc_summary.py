#!/usr/bin/env python3
"""Cross-sample QC metric aggregation for any Omics subworkflow.

Collects per-sample QC artifacts from the pipeline's QC-producing stages
(raw/trimmed FastQC, alignment flagstat, markdup metrics, optional stage
metrics like MSI scores), aligns them into one sample x metric table, and
renders per-stage QC plots. Tolerates missing stages/files (cells stay NA)
because optional stages (msi/cnvkit/burden/...) are gated by node.py
outfiles registration, not by this script.

Inputs are stage-name -> directory mappings (--qc-dir, repeatable); inside
each directory files are looked up per sample via --sample-path patterns
(defaults match the Omics module layouts). Metrics carried into the table
are declared with --metric stage:column-name (repeatable).

Outputs:
  <outdir>/qc_summary.tsv    sample x metric table (one row per sample)
  <outdir>/qc_summary.html   multi-section visual summary (optional)

Usage:
    python qc_summary.py \
        --sample S1 --sample S2 \
        --qc-dir fastqc_raw=/out/quality/fastqc/raw \
        --qc-dir flagstat=/out/common/3_sorted_bam \
        --qc-dir markdup=/out/common/4_markdup_bam \
        --qc-dir msi=/out/stability/msi \
        --outdir /out/QC/qc_summary
"""

import argparse
import html
import json
import os
import re
import sys
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# setup_logger — verbatim copy of src/common/util/LogUtil.py::setup_logger
# (modules.md 脚本撰写规范 #0: module bin scripts are self-contained)
# ---------------------------------------------------------------------------
from typing import Optional as _Optional, TextIO as _TextIO
import logging as _logging


def setup_logger(
    logger_name: str,
    level: int = _logging.INFO,
    log_file: _Optional[str] = None,
    stream: _Optional[_TextIO] = None,
) -> _logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = _logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = _logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = _logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = _logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


# ============================================================
# Per-stage parsers: stage dir + sample -> metric dict
# ============================================================
# Each parser receives the stage directory and a sample id, and returns a
# flat {metric_name: value} dict (values numeric or "NA"). Missing files
# yield {} so the sample row keeps NA for that stage.

# samtools flagstat summary line -> named counters (order = flagstat output)
_FLAGSTAT_FIELDS = [
    "total", "secondary", "supplementary", "duplicates", "mapped",
    "paired", "read1", "read2", "properly_paired", "with_itself",
    "singletons", "mate_mapped_diff_chr", "mate_mapped_diff_chr_mapq5",
]

# samtools flagstat label -> metric key. The whitelist above uses short
# keys, but real flagstat labels read "in total", "paired in sequencing",
# "with itself and mate mapped", "with mate mapped to a different chr".
_FLAGSTAT_LABEL_ALIASES = {
    "in_total": "total",
    "paired_in_sequencing": "paired",
    "with_itself_and_mate_mapped": "with_itself",
    "with_mate_mapped_to_a_different_chr": "mate_mapped_diff_chr",
}

# GATK MarkDuplicates metrics file -> kept columns (LIBRARY row)
_MARKDUP_FIELDS = [
    "READ_PAIRS_EXAMINED", "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES",
    "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE",
]


def parse_fastqc(stage_dir: str, sample: str) -> Dict[str, str]:
    """Extract per-module stats from a FastQC fastqc_data.txt.

    Reads Total Sequences / Sequence length / %GC and the per-module status
    (basic_statistics, per_base_quality, adapter_content ...). Works for
    both paired reads' zip-extracted data (first _1 file wins) and the
    plain fastqc_data.txt layout.

    Parameters
    ----------
    stage_dir : str
        FastQC output directory (contains <sample>/ subdirs or *_fastqc dirs).
    sample : str
        Sample id.

    Returns
    -------
    Dict[str, str]
        Metric -> value; empty when no fastqc_data.txt is found.
    """
    # candidate data files: <sample>/<sample>_1_fastqc/fastqc_data.txt etc.
    cands: List[str] = []
    sample_dir = os.path.join(stage_dir, sample)
    if os.path.isdir(sample_dir):
        for entry in sorted(os.listdir(sample_dir)):
            if entry.endswith("_fastqc"):
                cands.append(os.path.join(sample_dir, entry, "fastqc_data.txt"))
        cands.append(os.path.join(sample_dir, "fastqc_data.txt"))
    # fallback: flat stage dir with <sample>*_fastqc subdirs
    if not any(os.path.isfile(c) for c in cands) and os.path.isdir(stage_dir):
        for entry in sorted(os.listdir(stage_dir)):
            if entry.startswith(sample) and entry.endswith("_fastqc"):
                cands.append(os.path.join(stage_dir, entry, "fastqc_data.txt"))
    data_file = next((c for c in cands if os.path.isfile(c)), None)
    if not data_file:
        return {}

    metrics: Dict[str, str] = {}
    with open(data_file) as fh:
        in_modules = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("Total Sequences"):
                metrics["total_sequences"] = line.split("\t")[-1]
            elif line.startswith("Sequence length"):
                metrics["sequence_length"] = line.split("\t")[-1]
            elif line.startswith("%GC"):
                metrics["percent_gc"] = line.split("\t")[-1]
            elif line.startswith(">>END_MODULE"):
                in_modules = False
            elif in_modules and line.startswith("#"):
                continue
            elif in_modules and "\t" in line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    key = re.sub(r"[^a-z0-9]+", "_", parts[0].strip().lower()).strip("_")
                    metrics[f"fastqc_{key}"] = parts[1].strip()
            elif line.startswith(">>"):
                in_modules = True
                mod = line[2:].split("\t")[0].strip().lower()
                mod_key = re.sub(r"[^a-z0-9]+", "_", mod).strip("_")
                metrics[f"fastqc_status_{mod_key}"] = line.split("\t")[-1].strip()
    return metrics


def parse_cutadapt(stage_dir: str, sample: str) -> Dict[str, str]:
    """Parse cutadapt full report (<sample>.cutadapt_report.txt paired,
    <sample>.single.cutadapt_report.txt single): read throughput and
    trimming percentages from the Summary block. Unknown trailing lines
    (the wrapper script's echo output) are ignored."""
    path = None
    for cand in (
        os.path.join(stage_dir, sample, f"{sample}.cutadapt_report.txt"),
        os.path.join(stage_dir, sample, f"{sample}.single.cutadapt_report.txt"),
        os.path.join(stage_dir, f"{sample}.cutadapt_report.txt"),
        os.path.join(stage_dir, f"{sample}.single.cutadapt_report.txt"),
    ):
        if os.path.isfile(cand):
            path = cand
            break
    if not path:
        return {}

    metrics: Dict[str, str] = {}
    pct_re = re.compile(r"\(([\d.]+)%\)")
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^Total read pairs processed:\s*([\d,]+)", line)
            if m:
                metrics["cutadapt_read_pairs"] = m.group(1).replace(",", "")
                continue
            m = re.match(r"^Total reads processed:\s*([\d,]+)", line)
            if m:
                metrics["cutadapt_reads"] = m.group(1).replace(",", "")
                continue
            m = re.match(r"^Pairs written \(passing filters\):\s*([\d,]+)", line)
            if m:
                pct = pct_re.search(line)
                metrics["cutadapt_pairs_written_pct"] = pct.group(1) if pct else ""
                continue
            m = re.match(r"^Reads written \(passing filters\):\s*([\d,]+)", line)
            if m:
                pct = pct_re.search(line)
                metrics["cutadapt_reads_written_pct"] = pct.group(1) if pct else ""
                continue
            m = re.match(r"^Quality-trimmed:\s*[\d,]+ bp", line)
            if m:
                pct = pct_re.search(line)
                metrics["cutadapt_quality_trimmed_pct"] = pct.group(1) if pct else ""
                continue
            m = re.match(r"^Total written \(filtered\):\s*[\d,]+ bp", line)
            if m:
                pct = pct_re.search(line)
                metrics["cutadapt_bases_written_pct"] = pct.group(1) if pct else ""
    # drop empty extractions (report present but block missing)
    return {k: v for k, v in metrics.items() if v != ""}


def parse_bwa_mem2(stage_dir: str, sample: str) -> Dict[str, str]:
    """Parse bwa-mem2 alignment log (bwa_mem2-alignment.txt in the rule-log
    dir): PE insert-size distribution inferred from the dominant FR
    orientation — mean/std/median, plus candidate unique-pair count.

    Only the FIRST FR block is taken (deterministic; the first chunk's
    inference already uses ~1M pairs). SE samples have no PE lines and
    keep NA. Non-metric bwa-mem2 chatter is ignored.
    """
    path = None
    for cand in (
        os.path.join(stage_dir, sample, "bwa_mem2-alignment.txt"),
        os.path.join(stage_dir, sample, "bwa-mem2-alignment.txt"),
        os.path.join(stage_dir, sample, "bwaMem2_alignment.log"),
        os.path.join(stage_dir, f"{sample}-bwa_mem2-alignment.txt"),
    ):
        if os.path.isfile(cand):
            path = cand
            break
    if not path:
        return {}

    metrics: Dict[str, str] = {}
    in_fr = False
    with open(path) as fh:
        for line in fh:
            if "analyzing insert size distribution for orientation FR" in line:
                in_fr = True
                continue
            if not in_fr:
                # candidate unique pairs: (FF, FR, RF, RR) — FR is 2nd
                m = re.search(r"\[PE\] # candidate unique pairs.*?: \(\d+, (\d+),", line)
                if m:
                    metrics.setdefault("bwa_unique_pairs_fr", m.group(1))
                continue
            # "(25, 50, 75) percentile: (152, 186, 227)" — median is the
            # SECOND value of the VALUE tuple (186), not the label tuple
            m = re.search(r"percentile: \(\d+, (\d+), \d+\)", line)
            if m:
                metrics["insert_size_median"] = m.group(1)
                continue
            # anchor on "[PE] mean and std" so the "boundaries for computing
            # mean and std.dev: (low, high)" line cannot match
            m = re.search(r"\[PE\] mean and std\.dev: \(([\d.]+), ([\d.]+)\)", line)
            if m:
                metrics["insert_size_mean"] = m.group(1)
                metrics["insert_size_std"] = m.group(2)
                break  # first FR block done
    return metrics


def parse_flagstat(stage_dir: str, sample: str) -> Dict[str, str]:
    """Parse samtools flagstat output (<sample>_flagstat.txt preferred,
    <sample>.flagstat fallback) into named counters plus mapped_rate."""
    path = None
    for cand in (
        os.path.join(stage_dir, sample, f"{sample}_flagstat.txt"),
        os.path.join(stage_dir, sample, f"{sample}.flagstat"),
        os.path.join(stage_dir, f"{sample}_flagstat.txt"),
        os.path.join(stage_dir, f"{sample}.flagstat"),
    ):
        if os.path.isfile(cand):
            path = cand
            break
    if not path:
        return {}

    # summary lines: "<n> + <n> <label> (<pct>%)". Non-percentage lines
    # ("primary", "read1", "read2", "paired in sequencing") have no trailing
    # "(...)" block, so accept the bare "<n> + <n> <label>" form too; keep
    # the FIRST n per label (skips the "primary ..." repeats).
    no_pct = re.compile(r"^(\d+) \+ \d+\s+([A-Za-z0-9_% ()/-]+?)\s*$")
    metrics: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^(\d+) \+ \d+ ([A-Za-z0-9_% ()/-]+?)\s+\(", line)
            if not m:
                m = no_pct.match(line)
            if not m:
                continue
            count, label = m.group(1), m.group(2).strip()
            key = re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")
            key = _FLAGSTAT_LABEL_ALIASES.get(key, key)
            if key in _FLAGSTAT_FIELDS and key not in metrics:
                metrics[key] = count
    if "mapped" in metrics and "total" in metrics and int(metrics["total"]):
        metrics["mapped_rate"] = f"{int(metrics['mapped']) / int(metrics['total']) * 100:.2f}"
    return metrics


def parse_markdup(stage_dir: str, sample: str) -> Dict[str, str]:
    """Parse GATK MarkDuplicates metrics (<sample>.Markdup-metrics.txt),
    keeping the LIBRARY-row duplication metrics."""
    path = os.path.join(stage_dir, sample, f"{sample}.Markdup-metrics.txt")
    if not os.path.isfile(path):
        path = os.path.join(stage_dir, f"{sample}.Markdup-metrics.txt")
    if not os.path.isfile(path):
        return {}

    metrics: Dict[str, str] = {}
    with open(path) as fh:
        in_hist = False
        header: List[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("## HISTOGRAM"):
                break
            if line.startswith("LIBRARY\t"):
                header = line.split("\t")[1:]
                in_hist = False
                continue
            if header and not in_hist and "\t" in line and not line.startswith("#"):
                values = line.split("\t")
                for key, val in zip(header, values[1:]):
                    if key in _MARKDUP_FIELDS and val != "":
                        metrics[key.lower()] = val
                in_hist = True
    return metrics


def parse_msi(stage_dir: str, sample: str) -> Dict[str, str]:
    """Parse msisensor-pro msi output (<sample>_msi_pro.tsv): 3-column
    row (total, evaluable, score) -> msi_score."""
    path = os.path.join(stage_dir, sample, f"{sample}_msi_pro.tsv")
    if not os.path.isfile(path):
        path = os.path.join(stage_dir, f"{sample}_msi_pro.tsv")
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            fields = line.split()
            if len(fields) >= 3:
                return {"msi_score": fields[2]}
    return {}


def parse_sexcheck(stage_dir: str, sample: str) -> Dict[str, str]:
    """Parse the sexcheck module's per-run table (sexcheck.tsv): one row
    per sample (sample_id, x/y coverage metrics, ..., sex). This stage is
    run-level (ONE table for all samples), unlike the per-sample parsers
    above — we read the whole table once per call and return only the
    queried sample's row (call cost is negligible; the table is small).
    """
    path = os.path.join(stage_dir, "sexcheck.tsv")
    if not os.path.isfile(path):
        return {}
    metrics: Dict[str, str] = {}
    with open(path) as fh:
        header: List[str] = []
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if not header:
                header = fields
                continue
            if fields and fields[0] == sample:
                for key, val in zip(header[1:], fields[1:]):
                    metrics[f"sex_{key}"] = val
                break
    return metrics


# stage name -> parser (unknown stages fall back to fastqc layout probing)
PARSERS = {
    "fastqc_raw": parse_fastqc,
    "fastqc_trimmed": parse_fastqc,
    "fastqc": parse_fastqc,
    "cutadapt": parse_cutadapt,
    "flagstat": parse_flagstat,
    "bwa_mem2": parse_bwa_mem2,
    "markdup": parse_markdup,
    "msi": parse_msi,
    "sexcheck": parse_sexcheck,
}


def collect_stage(stage: str, stage_dir: str, samples: List[str]) -> Dict[str, Dict[str, str]]:
    """Run the stage parser for every sample; log and skip broken files.

    Returns {sample: {metric: value}}; samples with no artifact keep {}.
    """
    parser = PARSERS.get(stage, parse_fastqc)
    out: Dict[str, Dict[str, str]] = {}
    for sample in samples:
        try:
            out[sample] = parser(stage_dir, sample)
        except Exception as exc:  # a malformed file must not sink the table
            print(f"[warn] {stage}/{sample}: parse failed ({exc}); leaving NA", file=sys.stderr)
            out[sample] = {}
    return out


# ============================================================
# Table assembly
# ============================================================

def build_table(samples: List[str], stages: Dict[str, Dict[str, Dict[str, str]]]) -> pd.DataFrame:
    """Merge per-stage metric dicts into one sample x metric DataFrame.

    Column order: fixed core columns first (total reads, mapped rate,
    duplication), then remaining columns in stage/first-seen order. Missing
    cells become "NA".
    """
    rows = []
    for sample in samples:
        row = {"sample_id": sample}
        for stage in stages:
            row.update(stages[stage].get(sample, {}))
        rows.append(row)
    table = pd.DataFrame(rows).set_index("sample_id").fillna("NA")

    core_order = [
        "total_sequences", "sequence_length", "percent_gc",
        "cutadapt_read_pairs", "cutadapt_reads",
        "cutadapt_pairs_written_pct", "cutadapt_reads_written_pct",
        "cutadapt_quality_trimmed_pct", "cutadapt_bases_written_pct",
        "total", "mapped", "mapped_rate", "properly_paired",
        "insert_size_mean", "insert_size_median", "insert_size_std",
        "bwa_unique_pairs_fr",
        "percent_duplication", "read_pair_duplicates", "msi_score",
        "sex_sex", "sex_ratio_x_to_y",
    ]
    cols = [c for c in core_order if c in table.columns]
    cols += [c for c in table.columns if c not in cols]
    return table.loc[:, cols]


def to_num(series: pd.Series) -> pd.Series:
    """Coerce a metric column to numeric for plotting (NA -> NaN)."""
    coerced = pd.to_numeric(series, errors="coerce")
    return coerced if isinstance(coerced, pd.Series) else pd.Series([coerced])


# ============================================================
# Plots
# ============================================================

def render_plots(table: pd.DataFrame, outdir: str, plots: List[str]) -> List[str]:
    """Render per-stage QC plots; each plot survives empty data by writing
    a placeholder image so the report keeps its layout.

    Returns the list of written PNG paths (only requested + non-empty ones).
    """
    written: List[str] = []
    for name in plots:
        path = os.path.join(outdir, f"qc_summary_{name}.png")
        if name == "mapped_rate":
            vals = to_num(table["mapped_rate"]) if "mapped_rate" in table.columns else pd.Series(dtype=float)
            if vals.notna().any():
                _bar_plot(vals.sort_index(), "Mapped reads (%)", "mapped rate", path)
                written.append(path)
        elif name == "duplication":
            vals = to_num(table["percent_duplication"]) if "percent_duplication" in table.columns else pd.Series(dtype=float)
            if vals.notna().any():
                _bar_plot((vals * (100 if vals.max() <= 1.0 else 1.0)).sort_index(),
                          "Duplication rate (%)", "duplication", path)
                written.append(path)
        elif name == "total_sequences":
            vals = to_num(table["total_sequences"]) if "total_sequences" in table.columns else pd.Series(dtype=float)
            if vals.notna().any():
                _bar_plot((vals / 1e6).sort_index(), "Raw reads (M)", "total sequences", path)
                written.append(path)
    return written


def _bar_plot(vals: pd.Series, ylabel: str, title: str, path: str) -> None:
    """Render one horizontal bar chart of per-sample values."""
    fig, ax = plt.subplots(figsize=(7, max(2.0, 0.45 * len(vals) + 1.0)))
    vals.plot.barh(ax=ax, color="#0094D8")
    ax.set_xlabel(ylabel)
    ax.set_title(f"QC summary — {title}")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ============================================================
# HTML rendering
# ============================================================

_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body {{ font-family: "Noto Sans CJK SC", "PingFang SC", Arial, sans-serif; margin: 24px; color: #333; }}
h1 {{ font-size: 20px; border-bottom: 2px solid #182543; padding-bottom: 6px; }}
h2 {{ font-size: 15px; color: #182543; margin-top: 28px; }}
table {{ border-collapse: collapse; font-size: 12px; margin-top: 8px; }}
th, td {{ border: 1px solid #d5dbe3; padding: 4px 8px; text-align: left; }}
th {{ background: #182543; color: #fff; }}
tr:nth-child(even) td {{ background: #f1f5fa; }}
img {{ max-width: 900px; margin: 8px 0; display: block; border: 1px solid #e3e8ef; }}
.pass {{ color: #00A878; font-weight: 600; }}
.warn {{ color: #F09A36; font-weight: 600; }}
.fail {{ color: #D64545; font-weight: 600; }}
</style>
</head>
<body>
<h1>{title}</h1>
{body}
</body>
</html>
"""


def render_html(table: pd.DataFrame, pngs: List[str], out_path: str,
                title: str, meta: Optional[dict] = None) -> None:
    """Write the visual summary: metric table + embedded PNG sections."""
    parts: List[str] = []
    if meta:
        rows = "".join(
            f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>"
            for k, v in meta.items())
        parts.append(f"<h2>Run info</h2><table>{rows}</table>")
    html_table = table.to_html(escape=True)
    parts.append(f"<h2>Metric table</h2>{html_table}")
    if pngs:
        imgs = "".join(
            f'<h2>{html.escape(os.path.basename(p))}</h2><img src="{html.escape(os.path.basename(p))}">'
            for p in pngs)
        parts.append(imgs)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(_HTML_TEMPLATE.format(title=html.escape(title), body="\n".join(parts)))


# ============================================================
# CLI
# ============================================================

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-s", "--sample", action="append", default=[],
                   help="Sample id (repeatable)")
    p.add_argument("-q", "--qc-dir", action="append", default=[],
                   metavar="STAGE=DIR",
                   help="Stage name -> QC artifact directory (repeatable)")
    p.add_argument("-o", "--outdir", required=True,
                   help="Output directory for qc_summary.tsv/html/png")
    p.add_argument("-p", "--plot", action="append", default=["mapped_rate", "duplication", "total_sequences"],
                   choices=["mapped_rate", "duplication", "total_sequences"],
                   help="QC plots to render (repeatable)")
    p.add_argument("-t", "--title", default="QC Summary",
                   help="Report title shown in the HTML header")
    p.add_argument("-m", "--meta", default=None,
                   help="Optional JSON file with run metadata shown as a table")
    p.add_argument("-l", "--log", default=None,
                   help="Optional log file path (also logs to stdout)")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    logger = setup_logger("qc_summary", stream=sys.stdout,
                          log_file=args.log)
    if args.log:
        os.makedirs(os.path.dirname(os.path.abspath(args.log)), exist_ok=True)

    if not args.sample:
        logger.error("No samples given; pass --sample for every sample id")
        return 1
    stages: Dict[str, Dict[str, Dict[str, str]]] = {}
    for pair in args.qc_dir:
        if "=" not in pair:
            logger.error(f"Bad --qc-dir '{pair}'; expected STAGE=DIR")
            return 1
        stage, stage_dir = pair.split("=", 1)
        if not os.path.isdir(stage_dir):
            # optional stage not run in this DAG — keep the column set empty,
            # do not fail (node.py gates optional stages via outfiles)
            logger.warning(f"stage '{stage}' dir missing ({stage_dir}); skipping")
            continue
        stages[stage] = collect_stage(stage, stage_dir, args.sample)
        logger.info(f"stage '{stage}': "
                    + ", ".join(f"{s}={len(m)} metrics" for s, m in stages[stage].items()))

    table = build_table(args.sample, stages)
    os.makedirs(args.outdir, exist_ok=True)
    tsv_path = os.path.join(args.outdir, "qc_summary.tsv")
    table.to_csv(tsv_path, sep="\t")
    logger.info(f"wrote {tsv_path} ({table.shape[0]} samples x {table.shape[1]} metrics)")

    pngs = render_plots(table, args.outdir, args.plot)
    html_path = os.path.join(args.outdir, "qc_summary.html")
    meta = json.load(open(args.meta)) if args.meta and os.path.isfile(args.meta) else None
    render_html(table, pngs, html_path, args.title, meta)
    logger.info(f"wrote {html_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
