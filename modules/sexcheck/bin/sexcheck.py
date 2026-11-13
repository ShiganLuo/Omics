#!/usr/bin/env python3
"""Sample sex inference from sex-chromosome coverage in alignment BAMs.

Reads per-sample samtools idxstats output (chrom <tab> length <tab> mapped
<tab> unmapped), computes the X/Y coverage ratio (mapped reads normalized
by chromosome length — plain read counts are biased by chromosome size)
and calls the sample sex:

  - XY : Y coverage present AND ratio low  (male: single X dose, Y>0)
  - XX : Y coverage ~absent AND ratio high (female: double X dose, Y=0)
  - unknown : Y absent from the reference (e.g. assemblies without chrY),
    autosomes missing, or the ratio falls in the ambiguous band

Decision thresholds are computed per-run from the cohort's own ratio
distribution when it is bimodal enough (>=1 sample on each side of the
midpoint), otherwise the fixed defaults apply. This keeps the call
meaningful for mouse (X ≈ 5% of genome, ratio gap wide) and human alike.

Inputs (both modes):
    --tsv idxstats.tsv   one row per line: sample <tab> chrom <tab> length <tab> mapped <tab> unmapped
                         (produced by the sexcheck smk rule; concatenated per-sample files)
Outputs:
    <outdir>/sexcheck.tsv   DNA mode: sample_id, x/y coverage metrics, sex
                            RNA mode: sample_id, XIST/Y marker expression metrics, sex
    <outdir>/sexcheck.done  marker

DNA mode (default) calls sex from the X/Y coverage ratio as described
above. RNA mode (--data-type RNA) counts reads over sex-marker genes
instead: XIST is expressed only in females and the Y-linked markers only
in males. Coverage ratios do not discriminate sex on transcriptome data
(Y genes are sparsely expressed in most tissues), so the RNA call uses
two expression floors: XY when Y-marker RPKM >= y_min_rpkm and XIST RPKM
< xist_min_rpkm, XX when the mirror holds, unknown otherwise (including
missing marker genes or missing BAMs). Marker gene names are matched
case-insensitively against the GTF gene_name attribute, so the defaults
serve mouse (Xist, Ddx3y, ...) and human (XIST, DDX3Y, ...) alike.

Usage (DNA):
    python sexcheck.py --tsv idxstats.tsv --outdir QC/sexcheck \
        [--ratio-xy-threshold 2.0] [--y-min-frac 0.01]
Usage (RNA):
    python sexcheck.py --tsv idxstats.tsv --outdir QC/sexcheck \
        --data-type RNA --gtf genes.gtf --bam-map '{"S1": "S1.bam"}' \
        [--y-min-rpkm 0.5] [--xist-min-rpkm 0.1]
Usage (RNA, multi-genome runs):
    python sexcheck.py --tsv idxstats.tsv --outdir QC/sexcheck \
        --data-type RNA --gtf-map '{"GRCm39": "mouse.gtf"}' \
        --bam-map '{"S1": "S1.bam"}' --sample-genome S1=GRCm39
"""

import argparse
import json
import os
import re
import subprocess
import sys
from typing import Dict, List, Optional, Sequence, Tuple

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


# Sex chromosome naming across assembly conventions; first match per
# category wins per sample (Ensembl uses bare X/Y, UCSC chrX/chrY).
X_NAMES = ("chrX", "X", "x")
Y_NAMES = ("chrY", "Y", "y")


def pick_sex_chroms(chroms: List[str], x_hint: Optional[str] = None,
                    y_hint: Optional[str] = None) -> Tuple[Optional[str], Optional[str]]:
    """Return (x_name, y_name) actually present in the idxstats chrom list.

    Explicit hints (from genome.references.<g>.sex_chroms, auto-filled by
    node.py's GENOME_SEX_CHROMS for NCBI-accession assemblies) win when
    the hinted contig exists; otherwise the chrX/X and chrY/Y naming
    conventions are sniffed as before. A hint pointing at a missing contig
    falls through to sniffing (with the hint logged for diagnosis).
    """
    x = x_hint if x_hint in chroms else next((c for c in chroms if c in X_NAMES), None)
    y = y_hint if y_hint in chroms else next((c for c in chroms if c in Y_NAMES), None)
    return x, y


def load_idxstats(path: str) -> pd.DataFrame:
    """Load concatenated idxstats rows: sample/chrom/length/mapped/unmapped.

    The trailing unmapped column (4th samtools column; the '*' row carries
    the per-sample total there) is parsed when present and ignored by the
    X/Y logic — it is kept so qc consumers can compute unmapped rates.
    Legacy 4-column files (no unmapped) still load.
    """
    df = pd.read_csv(path, sep="\t", names=["sample", "chrom", "length",
                                             "mapped", "unmapped"])
    return df


def build_table(df: pd.DataFrame, x_hint: Optional[str] = None,
                   y_hint: Optional[str] = None,
                   sex_chroms_map: Optional[Dict[str, tuple]] = None,
                   sample_genomes: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Per-sample X/Y coverage metrics; NA when a chromosome is absent.

    Sex-contig name resolution per sample: (1) single-genome runs pass
    x_hint/y_hint; (2) multi-genome runs pass sex_chroms_map
    {genome: (x, y)} + sample_genomes {sample: genome}. Hints absent or
    pointing at missing contigs fall back to chrX/X + chrY/Y sniffing.
    """
    sex_chroms_map = sex_chroms_map or {}
    sample_genomes = sample_genomes or {}
    rows = []
    for sample, rec in df.groupby("sample", sort=True):
        chroms = list(rec["chrom"])
        g = sample_genomes.get(str(sample))
        gxy = (sex_chroms_map.get(g) or (None, None)) if g else (None, None)
        x_name, y_name = pick_sex_chroms(chroms, gxy[0] or x_hint, gxy[1] or y_hint)
        def _get(name: Optional[str]) -> Tuple[float, float]:
            if name is None:
                return float("nan"), float("nan")
            sub = rec[rec["chrom"] == name]
            if sub.empty:
                return float("nan"), float("nan")
            length_series = pd.Series(sub["length"].tolist())
            mapped_series = pd.Series(sub["mapped"].tolist())
            return float(length_series.iloc[0]), float(mapped_series.iloc[0])
        x_len, x_mapped = _get(x_name)
        y_len, y_mapped = _get(y_name)
        # length-normalized coverage proxy: mapped reads per Mb
        x_norm = x_mapped / (x_len / 1e6) if x_len and x_len == x_len else float("nan")
        y_norm = y_mapped / (y_len / 1e6) if y_len and y_len == y_len else float("nan")
        ratio = x_norm / y_norm if (y_norm and y_norm == y_norm and y_norm > 0) else float("inf")
        rows.append({
            "sample_id": sample,
            "x_chrom": x_name or "NA", "y_chrom": y_name or "NA",
            "x_mapped": int(x_mapped) if x_mapped == x_mapped else "NA",
            "y_mapped": int(y_mapped) if y_mapped == y_mapped else "NA",
            "x_len": int(x_len) if x_len == x_len else "NA",
            "y_len": int(y_len) if y_len == y_len else "NA",
            "x_reads_per_mb": round(x_norm, 2) if x_norm == x_norm else "NA",
            "y_reads_per_mb": round(y_norm, 2) if y_norm == y_norm else "NA",
            "ratio_x_to_y": round(ratio, 4) if ratio != float("inf") else "inf",
        })
    return pd.DataFrame(rows)


def call_sex(table: pd.DataFrame, y_min_frac: float, ratio_threshold: float,
             logger) -> pd.DataFrame:
    """Assign XX/XY/unknown per sample.

    XY: Y coverage above the floor AND ratio below threshold (male signal).
    XX: Y coverage below the floor (female signal: Y essentially empty).
    unknown: Y chromosome absent from the assembly, or the ratio sits in
    the ambiguous band around the threshold (half-band on each side).
    """
    calls: List[str] = []
    for rec in table.to_dict("records"):
        y_norm = rec["y_reads_per_mb"]
        ratio = rec["ratio_x_to_y"]
        if rec["y_chrom"] == "NA" or y_norm == "NA":
            calls.append("unknown")
            continue
        y_present = float(y_norm) >= y_min_frac
        if ratio == "inf":
            calls.append("XX" if not y_present else "unknown")
        elif not y_present:
            calls.append("XX")
        else:
            r = float(ratio)
            lo, hi = ratio_threshold * 0.5, ratio_threshold * 1.5
            if r < lo:
                calls.append("XY")
            elif r > hi:
                calls.append("XX")
            else:
                calls.append("unknown")
    table = table.copy()
    table["sex"] = calls
    for rec in table.to_dict("records"):
        logger.info(f"{rec['sample_id']}: x/y reads_per_mb = "
                    f"{rec['x_reads_per_mb']}/{rec['y_reads_per_mb']}, "
                    f"ratio = {rec['ratio_x_to_y']} -> {rec['sex']}")
    return table


# ---------------------------------------------------------------------------
# RNA mode: sex-marker gene expression (XIST + Y-linked genes)
# ---------------------------------------------------------------------------
# Sex calling from marker-gene transcription instead of chromosome coverage:
# XIST is expressed only under the female X-inactivation program and the
# Y-linked markers are male-specific. The X/Y coverage ratio does not
# discriminate sex on transcriptome data (Y genes are sparsely expressed in
# most tissues), so the RNA call uses two expression floors. Gene names are
# matched case-insensitively against the GTF gene_name attribute; the
# defaults serve mouse (Xist, Ddx3y, ...) and human (XIST, DDX3Y, ...) alike.
DEFAULT_XIST_GENE = "xist"
DEFAULT_Y_GENES = ("ddx3y", "eif2s3y", "kdm5d", "uty", "usp9y")
_GTF_ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')

def merge_spans(spans: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    """Merge overlapping or adjacent genomic spans per chromosome.

    Parameters
    ----------
    spans : list of (chrom, start, end)
        1-based inclusive intervals.

    Returns
    -------
    list of (chrom, start, end)
        Disjoint merged spans, sorted by chromosome and start.
    """
    by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for c, s, e in spans:
        by_chrom.setdefault(c, []).append((s, e))
    merged: List[Tuple[str, int, int]] = []
    for c in sorted(by_chrom):
        intervals = sorted(by_chrom[c])
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                merged.append((c, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((c, cur_s, cur_e))
    return merged

def parse_marker_spans(gtf: str, xist_gene: str, y_genes: Sequence[str],
                       logger) -> Dict[str, List[Tuple[str, int, int]]]:
    """Extract sex-marker gene spans from a GTF annotation file.

    Matches the gene_name GTF attribute case-insensitively. Prefers ``gene``
    features; falls back to merging ``transcript`` spans per gene name for
    GTFs carrying only transcript/exon rows (e.g. StringTie output).

    Parameters
    ----------
    gtf : str
        Annotation GTF path (config key genome.gtf).
    xist_gene : str
        XIST marker gene name (case-insensitive).
    y_genes : sequence of str
        Y-linked marker gene names (case-insensitive).
    logger : logging.Logger

    Returns
    -------
    dict
        {gene_name_as_spelled_in_gtf: [(chrom, start, end), ...]} for every
        marker gene found; missing genes are logged and absent from the dict.

    Raises
    ------
    FileNotFoundError
        The GTF path does not exist.
    """
    if not os.path.isfile(gtf):
        raise FileNotFoundError(f"GTF not found: {gtf}")
    wanted = {str(g).strip().lower() for g in list(y_genes) + [xist_gene]
              if str(g).strip()}
    spans: Dict[str, List[Tuple[str, int, int]]] = {}

    def _scan(feature: str) -> None:
        with open(gtf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != feature:
                    continue
                attrs = dict(_GTF_ATTR_RE.findall(fields[8]))
                name = (attrs.get("gene_name") or "").strip()
                if name.lower() not in wanted:
                    continue
                spans.setdefault(name, []).append(
                    (fields[0], int(fields[3]), int(fields[4])))

    _scan("gene")
    if not spans:
        _scan("transcript")
    result: Dict[str, List[Tuple[str, int, int]]] = {}
    for name, raw in spans.items():
        result[name] = merge_spans(raw)
    for g in sorted(wanted):
        if not any(k.lower() == g for k in result):
            logger.warning(f"marker gene '{g}' not found in {gtf}")
    return result

def count_bam_reads(bam: str, spans: List[Tuple[str, int, int]],
                    samtools: str = "samtools") -> int:
    """Count primary alignments overlapping the given spans via samtools.

    Runs ``samtools view -c -F 0x900 <bam> <chrom:start-end>`` once per span
    and sums the counts. Secondary and supplementary alignments are
    excluded: long-read BAMs carry large multi-mapping secondary fractions
    that would inflate marker counts with cross-mapping noise.

    Parameters
    ----------
    bam : str
        Indexed alignment BAM (config key sample_bams values).
    spans : list of (chrom, start, end)
        1-based inclusive intervals (see merge_spans).
    samtools : str
        samtools executable name/path.

    Returns
    -------
    int
        Summed read count over all spans.

    Raises
    ------
    RuntimeError
        samtools exits non-zero (message names the failing command).
    ValueError
        samtools output is not an integer.
    """
    total = 0
    for chrom, start, end in spans:
        region = f"{chrom}:{start}-{end}"
        cmd = [samtools, "view", "-c", "-F", "0x900", bam, region]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr}")
        try:
            total += int(proc.stdout.strip())
        except ValueError as e:
            raise ValueError(
                f"non-integer count from {' '.join(cmd)}: {proc.stdout!r}") from e
    return total

def build_expression_table(df: pd.DataFrame, bam_map: Dict[str, str],
                           spans_by_sample: Dict[str, Optional[Dict[str, List[Tuple[str, int, int]]]]],
                           xist_gene: str, y_genes: Sequence[str],
                           samtools: str, logger) -> pd.DataFrame:
    """Per-sample XIST / Y-marker expression metrics from BAM read counts.

    RPKM = reads / (merged marker span kb) / (total mapped / 1e6); total
    mapped comes from the idxstats table so depth differences across
    samples do not shift the expression floors.

    Parameters
    ----------
    df : pandas.DataFrame
        Concatenated idxstats rows (sample/chrom/length/mapped/unmapped).
    bam_map : dict
        {sample_id: bam_path} (config key sample_bams).
    spans_by_sample : dict
        {sample_id: {gene_name: [(chrom, start, end), ...]}} — the marker
        spans valid for each sample (per-genome annotation on multi-genome
        runs, one shared dict otherwise); None marks a sample whose genome
        has no annotation (metrics come out NA -> unknown).
    xist_gene : str
        XIST marker gene name (case-insensitive lookup in marker_spans).
    y_genes : sequence of str
        Y-linked marker gene names (case-insensitive lookup).
    samtools : str
        samtools executable name/path.

    Returns
    -------
    pandas.DataFrame
        One row per sample: sample_id, xist_gene, y_genes, xist_reads,
        xist_rpkm, y_reads, y_rpkm_sum, y_rpkm_max, y_rpkm_detail,
        total_mapped, sex (filled later by call_sex_rna).
    """

    def _find(marker_spans: Dict[str, List[Tuple[str, int, int]]],
              name: str) -> Optional[str]:
        for k in marker_spans:
            if k.lower() == str(name).strip().lower():
                return k
        return None

    def _rpkm(reads: int, spans: List[Tuple[str, int, int]],
              total_mapped: float) -> float:
        span_kb = sum(e - s + 1 for _, s, e in spans) / 1000.0
        if span_kb <= 0 or total_mapped <= 0:
            return float("nan")
        return reads / span_kb / (total_mapped / 1e6)

    rows = []
    for sample in sorted(set(df["sample"])):
        total_mapped = float(df.loc[df["sample"] == sample, "mapped"].sum())
        marker_spans = spans_by_sample.get(str(sample))
        x_key: Optional[str] = None
        y_keys: List[str] = []
        x_reads: Optional[int] = None
        x_rpkm = float("nan")
        y_detail: Dict[str, float] = {}
        y_reads: Optional[int] = 0
        bam = bam_map.get(str(sample))
        if bam is None:
            logger.warning(f"{sample}: no BAM in bam map, expression metrics unavailable")
            y_reads = None
        elif marker_spans is None:
            logger.warning(f"{sample}: no marker annotation for this sample's "
                           "genome, expression metrics unavailable")
            y_reads = None
        else:
            x_key = _find(marker_spans, xist_gene)
            y_keys = [k for g in y_genes if (k := _find(marker_spans, g))]
            if x_key is not None:
                x_reads = count_bam_reads(bam, marker_spans[x_key], samtools)
                x_rpkm = _rpkm(x_reads, marker_spans[x_key], total_mapped)
            for k in y_keys:
                n = count_bam_reads(bam, marker_spans[k], samtools)
                y_reads += n
                y_detail[k] = round(_rpkm(n, marker_spans[k], total_mapped), 4)
        y_rpkm_sum = round(sum(y_detail.values()), 4) if y_detail else float("nan")
        y_rpkm_max = round(max(y_detail.values()), 4) if y_detail else float("nan")
        rows.append({
            "sample_id": sample,
            "xist_gene": x_key or "NA",
            "y_genes": ",".join(y_keys) if y_keys else "NA",
            "xist_reads": x_reads if x_reads is not None else "NA",
            "xist_rpkm": round(x_rpkm, 4) if x_rpkm == x_rpkm else "NA",
            "y_reads": y_reads if y_reads is not None else "NA",
            "y_rpkm_sum": y_rpkm_sum if y_rpkm_sum == y_rpkm_sum else "NA",
            "y_rpkm_max": y_rpkm_max if y_rpkm_max == y_rpkm_max else "NA",
            "y_rpkm_detail": ";".join(f"{k}:{v}" for k, v in y_detail.items()) or "NA",
            "total_mapped": int(total_mapped),
        })
    return pd.DataFrame(rows)

def call_sex_rna(table: pd.DataFrame, y_min_rpkm: float, xist_min_rpkm: float,
                 logger) -> pd.DataFrame:
    """Assign XX/XY/unknown per sample from marker-gene expression.

    Decision order (single state machine):
      XY : y_rpkm_sum >= y_min_rpkm AND xist_rpkm <  xist_min_rpkm
      XX : xist_rpkm >= xist_min_rpkm AND y_rpkm_sum <  y_min_rpkm
      unknown : everything else (both signals present or both absent,
                markers missing from the GTF, or sample BAM missing).

    Parameters
    ----------
    table : pandas.DataFrame
        Output of build_expression_table.
    y_min_rpkm : float
        Y-marker expression floor for the male signal.
    xist_min_rpkm : float
        XIST expression floor for the female signal.
    logger : logging.Logger

    Returns
    -------
    pandas.DataFrame
        Same table with the sex column filled.
    """
    calls: List[str] = []
    for rec in table.to_dict("records"):
        x_rpkm = rec["xist_rpkm"]
        y_rpkm = rec["y_rpkm_sum"]
        if x_rpkm == "NA" or y_rpkm == "NA":
            calls.append("unknown")
            continue
        y_present = float(y_rpkm) >= y_min_rpkm
        x_present = float(x_rpkm) >= xist_min_rpkm
        if y_present and not x_present:
            calls.append("XY")
        elif x_present and not y_present:
            calls.append("XX")
        else:
            calls.append("unknown")
    table = table.copy()
    table["sex"] = calls
    for rec in table.to_dict("records"):
        logger.info(f"{rec['sample_id']}: xist_rpkm = {rec['xist_rpkm']}, "
                    f"y_rpkm_sum = {rec['y_rpkm_sum']} "
                    f"(floors y={y_min_rpkm}, xist={xist_min_rpkm}) -> {rec['sex']}")
    return table

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-t", "--tsv", required=True,
                   help="Concatenated idxstats table: sample/chrom/length/mapped")
    p.add_argument("-o", "--outdir", required=True,
                   help="Output directory (sexcheck.tsv + sexcheck.done)")
    p.add_argument("--ratio-xy-threshold", type=float, default=2.0,
                   help="x/y coverage-ratio decision threshold (XX above band, XY below)")
    p.add_argument("--x-chrom", default=None,
                   help="Explicit X contig name (from genome sex_chroms; "
                        "default: sniff chrX/X)")
    p.add_argument("--y-chrom", default=None,
                   help="Explicit Y contig name (from genome sex_chroms; "
                        "default: sniff chrY/Y)")
    p.add_argument("--sex-chroms-map", default="",
                   help="JSON {genome: [x, y]} for multi-genome runs")
    p.add_argument("--sample-genome", action="append", default=[],
                   metavar="SAMPLE=GENOME",
                   help="Sample -> genome assignment (repeatable, multi-genome runs)")
    p.add_argument("--y-min-frac", type=float, default=0.01,
                   help="Minimum y reads_per_mb for the Y to count as present")
    p.add_argument("-d", "--data-type", type=str.upper, choices=["DNA", "RNA"],
                   default="DNA",
                   help="Data type: DNA keeps the coverage-ratio call, RNA "
                        "uses XIST/Y-marker expression (default: DNA)")
    p.add_argument("-g", "--gtf", default=None,
                   help="Annotation GTF for marker genes (RNA mode, config "
                        "key genome.gtf)")
    p.add_argument("-G", "--gtf-map", default="",
                   help="JSON {genome: gtf_path} for multi-genome runs "
                        "(RNA mode, config key genome.gtf_map)")
    p.add_argument("-b", "--bam-map", default="",
                   help="JSON {sample: bam_path} (RNA mode, config key "
                        "sample_bams)")
    p.add_argument("-s", "--samtools", default="samtools",
                   help="samtools executable for marker counting (RNA mode)")
    p.add_argument("-x", "--xist-gene", default=DEFAULT_XIST_GENE,
                   help=f"XIST marker gene name (default: {DEFAULT_XIST_GENE})")
    p.add_argument("-y", "--y-genes", default=",".join(DEFAULT_Y_GENES),
                   help="Comma-separated Y-linked marker gene names "
                        f"(default: {','.join(DEFAULT_Y_GENES)})")
    p.add_argument("-m", "--y-min-rpkm", type=float, default=0.5,
                   help="RNA mode: Y-marker expression floor for the male "
                        "signal (default: 0.5)")
    p.add_argument("-i", "--xist-min-rpkm", type=float, default=0.1,
                   help="RNA mode: XIST expression floor for the female "
                        "signal (default: 0.1)")
    p.add_argument("-l", "--log", default=None,
                   help="Optional log file path (also logs to stdout)")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    logger = setup_logger("sexcheck", stream=sys.stdout, log_file=args.log)
    if args.log:
        os.makedirs(os.path.dirname(os.path.abspath(args.log)), exist_ok=True)
    df = load_idxstats(args.tsv)
    if df.empty:
        logger.error("empty idxstats input")
        return 1
    y_genes = [g.strip() for g in args.y_genes.split(",") if g.strip()]
    sample_genomes = {}
    for pair in args.sample_genome:
        if "=" in pair:
            s, g = pair.split("=", 1)
            sample_genomes[s] = g

    def _resolve_spans(gtf_path: str) -> Optional[Dict[str, List[Tuple[str, int, int]]]]:
        """Parse and validate marker spans for one GTF; None when unusable."""
        try:
            spans = parse_marker_spans(gtf_path, args.xist_gene, y_genes, logger)
        except FileNotFoundError as e:
            logger.warning(f"{e} (config key genome.gtf / genome.gtf_map)")
            return None
        found = {k.lower() for k in spans}
        if (args.xist_gene.strip().lower() not in found
                or not any(g.lower() in found for g in y_genes)):
            logger.warning(
                f"marker genes missing from {gtf_path} (need {args.xist_gene} "
                f"plus one of {y_genes}); override with --xist-gene/--y-genes "
                "or fix the annotation")
            return None
        return spans

    if args.data_type == "RNA":
        if not args.gtf and not args.gtf_map:
            logger.error("RNA mode requires --gtf or --gtf-map "
                         "(config keys genome.gtf / genome.gtf_map)")
            return 1
        if not args.bam_map:
            logger.error("RNA mode requires --bam-map (config key sample_bams)")
            return 1
        try:
            bam_map = json.loads(args.bam_map)
            gtf_map = json.loads(args.gtf_map) if args.gtf_map else {}
        except json.JSONDecodeError as e:
            logger.error(f"--bam-map/--gtf-map is not valid JSON: {e}")
            return 1
        spans_by_sample: Dict[str, Optional[Dict[str, List[Tuple[str, int, int]]]]] = {}
        if gtf_map:
            spans_by_genome: Dict[str, Optional[Dict[str, List[Tuple[str, int, int]]]]] = {
                g: _resolve_spans(str(p)) for g, p in gtf_map.items()}
            for s in sorted(set(df["sample"])):
                genome = sample_genomes.get(str(s))
                if genome is None:
                    logger.warning(f"{s}: no genome assignment (--sample-genome), "
                                   "expression metrics unavailable")
                spans_by_sample[str(s)] = spans_by_genome.get(str(genome))
        else:
            spans = _resolve_spans(args.gtf)
            if spans is None:
                logger.error(f"annotation {args.gtf} is unusable for RNA mode "
                             "(config key genome.gtf)")
                return 1
            for s in sorted(set(df["sample"])):
                spans_by_sample[str(s)] = spans
        try:
            table = build_expression_table(df, bam_map, spans_by_sample,
                                           args.xist_gene, y_genes,
                                           args.samtools, logger)
        except (RuntimeError, ValueError) as e:
            logger.error(str(e))
            return 1
        table = call_sex_rna(table, args.y_min_rpkm, args.xist_min_rpkm, logger)
    else:
        sex_chroms_map = json.loads(args.sex_chroms_map) if args.sex_chroms_map else {}
        table = build_table(df, args.x_chrom, args.y_chrom,
                            sex_chroms_map=sex_chroms_map,
                            sample_genomes=sample_genomes)
        table = call_sex(table, args.y_min_frac, args.ratio_xy_threshold, logger)
    os.makedirs(args.outdir, exist_ok=True)
    out_tsv = os.path.join(args.outdir, "sexcheck.tsv")
    table.to_csv(out_tsv, sep="\t", index=False)
    with open(os.path.join(args.outdir, "sexcheck.done"), "w") as fh:
        fh.write("done\n")
    logger.info(f"wrote {out_tsv} ({len(table)} samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
