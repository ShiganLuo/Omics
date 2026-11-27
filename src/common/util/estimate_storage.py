#!/usr/bin/env python
"""Estimate the on-disk footprint of an Omics workflow run (final kept /
transient / peak), or measure and calibrate against a finished run.

Two kinds of run JSON exist (docs/sample-params-split.md):

  * INSTANTIATED json (raw.json, written by node.run<Workflow>()): one full
    run's parameters -- module switches (Params), genome, samples and every
    input file (raw_files) plus the registered deliverable list (outfiles).
    Passing it is enough; no meta is needed:
        python src/common/util/estimate_storage.py estimate --run-json OUTDIR/raw.json

  * UNINSTANTIATED json (the config/<Workflow>.json template): it holds no
    samples and no input files, so meta + workflow are required. The tool
    then reads config/<Workflow>.json (or --config) and INSTANTIATES it with
    the very same code path run.py uses (MetadataUtils + node.run<Workflow>),
    producing an instantiated json which drives the estimate:
        python src/common/util/estimate_storage.py estimate \\
            --meta data/meta_xxx.tsv --workflow PacVar
        # keep the generated instance for inspection/reuse:
        ... estimate --meta data/meta_xxx.tsv --workflow PacVar \\
            --outdir OUTDIR [--save-json PATH]

Other subcommands:
    # measure a run dir after (or mid-) execution; no manifest needed
    python src/common/util/estimate_storage.py measure --root OUTDIR

    # calibrate the ratio table from a real run (append-only history);
    # inputs read from the run's raw.json or from --meta
    python src/common/util/estimate_storage.py calibrate --root OUTDIR \\
        --run-json OUTDIR/raw.json
    python src/common/util/estimate_storage.py calibrate --root OUTDIR \\
        --meta data/meta_xxx.tsv --workflow PacVar

Behavior (estimate):
  1. Every input file is stat()ed for real size (missing file -> error, we
     never guess sizes); files are attributed to samples by sample-id match
     in the path, and split into input_bam / input_fastq / input_total.
  2. Product classes come from two sources:
       a. a curated per-workflow MANIFEST (intermediates included) when one
          exists (currently PacVar), with kept/temp semantics matching the
          modules: maybe_temp(...) products are temp while
          Params.<mod>.keep_intermediates == False and are auto-deleted by
          Snakemake after their last consumer;
       b. the instantiated json's registered outfiles: any outfile not
          covered by the manifest is sized with the generic per-file model
          (FALLBACK_* tables below). This is what makes EVERY workflow
          estimable -- node.run<Workflow>() registers the full deliverable
          set for all of them.
  3. Class size = basis (input bam / fastq / total / genome fasta / fixed)
     x ratio (output/basis). Manifest ratios are per-class totals; generic
     ratios are per FILE times the file's attributed basis.
  4. Peak = max simultaneous bytes over the stage timeline: a class is on
     disk for born <= s <= die (temp classes release after their last
     consumer stage; kept classes have die=None). Generic outfile classes
     are kept, born at the final stage. Workflows without a curated
     manifest therefore report final + overhead but no modeled temp peak.
  5. Ratios prefer real calibration: storage_ratio_history.tsv entries with
     the same workflow+class and >= MIN_POINTS points give p50 with p10/p90
     as the range; otherwise the empirical low/mid/high is used.

Output: per-class table + final kept / transient / peak totals with ranges;
--report-json saves the same numbers as JSON.

Principles (shared with calibrate_resources.py):
  - no guessing: stat inputs for real, error on missing files;
  - ranges instead of point estimates (low/mid/high, or p10/p50/p90 once
    calibrated);
  - append-only calibration history so every correction is auditable.

Explicitly NOT modeled: GATK/Snakemake java tmp-dir churn (under out/tmp,
cleaned after MarkDuplicates) and filesystem block overhead (du exceeds the
byte sum when many small files exist).

Adding a curated workflow manifest: add a class list to MANIFESTS (fields in
the ClassSpec comment) or merge one at runtime with --manifest extra.json
({"WorkflowName": [class dict, ...]}).
"""
import argparse
import copy
import csv
import datetime
import json
import os
import re
import shutil
import sys
import tempfile

MIN_POINTS = 3
HISTORY_NAME = "storage_ratio_history.tsv"
GENERIC_BORN = 6

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))


PACVAR_CLASSES = [
    dict(id="pbmm2_sorted_bam", desc="pbmm2 aligned BAM (bare output, never temp)",
         scope="sample", basis="input_bam", ratio=(0.95, 1.00, 1.10),
         gate=[], temp_when=None, born=2, die=None,
         patterns=["common/2_sorted_bam/**"]),
    dict(id="addReadsGroup_bam", desc="GATK add-read-group BAM (maybe_temp, deleted when keep_intermediates=false)",
         scope="sample", basis="input_bam", ratio=(0.95, 1.00, 1.10),
         gate=[], temp_when=("Params.gatk_prepare.keep_intermediates", False, True),
         born=3, die=4,
         patterns=["common/3_markdup_bam/**/*.addReadsGroup*"]),
    dict(id="markdup_bam", desc="MarkDuplicates final BAM + metrics (main deliverable)",
         scope="sample", basis="input_bam", ratio=(0.95, 1.00, 1.10),
         gate=[], temp_when=None, born=4, die=None,
         patterns=["common/3_markdup_bam/**/*.sorted_markdup*"]),
    dict(id="qc_sexcheck", desc="sex-check idxstats + table",
         scope="run", basis="fixed", ratio=(1e6, 5e7, 5e8),
         gate=[], temp_when=None, born=4, die=None,
         patterns=["QC/sexcheck/**"]),
    dict(id="qc_summary", desc="QC summary table + HTML",
         scope="run", basis="fixed", ratio=(1e6, 2e7, 2e8),
         gate=[], temp_when=None, born=4, die=None,
         patterns=["QC/qc_summary/**"]),
    dict(id="deepvariant_vcf", desc="DeepVariant SNV/indel VCF",
         scope="sample", basis="input_bam", ratio=(0.0002, 0.001, 0.003),
         gate=[("Params.snv.enabled", True, True),
               ("Params.workflow.snv_caller", "deepvariant", "deepvariant")],
         temp_when=None, born=5, die=None,
         patterns=["variation/germline_snv_indel/**/*.vcf.gz*"]),
    dict(id="deepvariant_examples", desc="DeepVariant make_examples intermediates (deleted when keep_intermediates=false)",
         scope="sample", basis="input_bam", ratio=(0.005, 0.02, 0.06),
         gate=[("Params.snv.enabled", True, True),
               ("Params.workflow.snv_caller", "deepvariant", "deepvariant")],
         temp_when=("Params.deepvariant.keep_intermediates", False, True),
         born=5, die=5,
         patterns=["variation/germline_snv_indel/**/make_examples*",
                   "variation/germline_snv_indel/**/*.tfrecord*"]),
    dict(id="pbsv_svsig", desc="pbsv discover signal files (registered deliverable, kept)",
         scope="sample", basis="fixed", ratio=(5e6, 3e7, 2e8),
         gate=[("Params.sv.enabled", True, True)],
         temp_when=None, born=5, die=None,
         patterns=["variation/germline_sv/**/*.svsig.gz*"]),
    dict(id="pbsv_vcf", desc="pbsv call SV VCF",
         scope="sample", basis="fixed", ratio=(2e6, 2e7, 1e8),
         gate=[("Params.sv.enabled", True, True)],
         temp_when=None, born=5, die=None,
         patterns=["variation/germline_sv/**/*.vcf.gz*"]),
    dict(id="hiphase_phased_vcf", desc="HiPhase phased VCFs (SNV + SV)",
         scope="sample", basis="fixed", ratio=(1e7, 1e8, 5e8),
         gate=[("Params.phase.enabled", True, True),
               ("Params.snv.enabled", True, True),
               ("Params.sv.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["variation/**/*phased.vcf.gz*"]),
    dict(id="sawfish_joint", desc="sawfish joint CNV call (cohort level)",
         scope="run", basis="fixed", ratio=(1e8, 5e8, 2e9),
         gate=[("Params.sawfish.enabled", True, False)],
         temp_when=None, born=6, die=None,
         patterns=["variation/germline_cnv/cnv/**"]),
    dict(id="sawfish_workdir", desc="sawfish per-sample work dirs (deleted when keep_intermediates=false)",
         scope="sample", basis="fixed", ratio=(1e8, 1e9, 5e9),
         gate=[("Params.sawfish.enabled", True, False)],
         temp_when=("Params.sawfish.keep_intermediates", False, True),
         born=5, die=6,
         patterns=["variation/germline_cnv/**"]),
    dict(id="bamCoverage_bw", desc="coverage bigWig (shared by decoil/regionplot)",
         scope="sample", basis="input_bam", ratio=(0.0005, 0.002, 0.006),
         gate=[],
         temp_when=None, born=5, die=None,
         patterns=["eccdna/coverage/**"]),
    dict(id="sv_go_kegg", desc="SV gene GO/KEGG analysis",
         scope="run", basis="fixed", ratio=(1e7, 5e7, 2e8),
         gate=[("Params.function.sv_go_kegg.enabled", True, True),
               ("Params.sv.enabled", True, True)],
         need_comparisons=True, temp_when=None, born=6, die=None,
         patterns=["**/sv_go_kegg/**"]),
    dict(id="decoil_out", desc="Decoil eccDNA reconstruction + visualization (incl. IGV HTML)",
         scope="sample", basis="fixed", ratio=(2e7, 1e8, 1e9),
         gate=[("Params.decoil.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["eccdna/**"]),
    dict(id="sv_analysis", desc="SV comparison analysis (merged/annotated/MEI/RepeatMasker/circos)",
         scope="run", basis="fixed", ratio=(2e8, 1e9, 5e9),
         gate=[("Params.sv.enabled", True, True)],
         need_comparisons=True, temp_when=None, born=6, die=None,
         patterns=["variation/germline_sv_analysis/**"]),
    dict(id="trgt_out", desc="TRGT repeat genotyping (needs repeat_bed)",
         scope="sample", basis="fixed", ratio=(1e7, 1e8, 1e9),
         gate=[("Params.repeat.enabled", True, True)],
         require=["genome_repeat_bed"], temp_when=None, born=5, die=None,
         patterns=["repeat/trgt/**"]),
    dict(id="hifiasm_assembly", desc="hifiasm assembly dirs (contigs + undeclared .bin/.gfa intermediates, counted together)",
         scope="sample", basis="input_bam", ratio=(0.05, 0.15, 0.30),
         gate=[],
         temp_when=None, born=2, die=None,
         patterns=["repeat/centromere/**/assembly/**"]),
    dict(id="centromere_out", desc="centromere analysis (RepeatMasker + stats)",
         scope="sample", basis="fixed", ratio=(1e7, 1e8, 1e9),
         gate=[("Params.telomere.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["repeat/centromere/**"]),
    dict(id="telomere_out", desc="telomere length/composition/region outputs",
         scope="run", basis="fixed", ratio=(5e7, 5e8, 2e9),
         gate=[("Params.telomere.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["telomere/**"]),
    dict(id="regionplot_out", desc="region plots",
         scope="run", basis="fixed", ratio=(5e6, 5e7, 5e8),
         gate=[("Params.regionplot.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["regionplot/**"]),
    dict(id="run_overhead", desc="raw.json/log/tmp/.snakemake/input-symlink dirs fixed overhead",
         scope="run", basis="fixed", ratio=(5e8, 2e9, 1e10),
         gate=[], temp_when=None, born=1, die=None,
         patterns=["log/**", "tmp/**", ".snakemake/**", "raw.json",
                   "common/raw.json", "common/1_raw_data/**"]),
]

PACVAR_GATE_ANY = {
    "bamCoverage_bw": [("Params.decoil.enabled", True, True),
                       ("Params.regionplot.enabled", True, True)],
    "hifiasm_assembly": [("Params.telomere.enabled", True, True),
                         ("Params.hifiasm.qc.enabled", True, True)],
}

MANIFESTS = {"PacVar": PACVAR_CLASSES}
GATE_ANY = {"PacVar": PACVAR_GATE_ANY}

FALLBACK_PATTERNS = [
    ("fb_bam", ["**/*.bam", "**/*.bam.bai", "**/*.bai", "**/*.cram", "**/*.crai", "**/*.pbi"]),
    ("fb_vcf", ["**/*.vcf", "**/*.vcf.gz", "**/*.vcf.gz.tbi", "**/*.vcf.gz.csi",
                "**/*.bcf", "**/*.bcf.csi"]),
    ("fb_bigwig", ["**/*.bw", "**/*.bigWig", "**/*.bigwig", "**/*.wig", "**/*.wig.gz"]),
    ("fb_fasta", ["**/*.fa", "**/*.fasta", "**/*.fa.gz", "**/*.fasta.gz",
                  "**/*.fa.fai", "**/*.fai", "**/*.gfa", "**/*.gfa.gz"]),
    ("fb_fastq", ["**/*.fastq", "**/*.fastq.gz", "**/*.fq", "**/*.fq.gz"]),
    ("fb_bed", ["**/*.bed", "**/*.bed.gz", "**/*.gtf", "**/*.gtf.gz",
                "**/*.gff", "**/*.gff3", "**/*.gff.gz"]),
    ("fb_h5", ["**/*.h5ad", "**/*.h5", "**/*.hdf5", "**/*.loom", "**/*.rds", "**/*.rdata"]),
    ("fb_office", ["**/*.pptx", "**/*.docx", "**/*.xlsx"]),
    ("fb_table", ["**/*.tsv", "**/*.csv", "**/*.txt", "**/*.json", "**/*.yaml", "**/*.yml"]),
    ("fb_figure", ["**/*.png", "**/*.pdf", "**/*.svg", "**/*.html", "**/*.htm",
                   "**/*.jpg", "**/*.jpeg", "**/*.gif"]),
    ("fb_gz", ["**/*.gz", "**/*.gzip", "**/*.zip", "**/*.tar", "**/*.tgz"]),
]

FALLBACK_MODEL = {
    "fb_fastq":  ("input_fastq", (0.08, 0.25, 0.50)),
    "fb_bam":    ("input_total", (0.60, 1.10, 1.60)),
    "fb_bigwig": ("input_total", (2e-4, 1.5e-3, 6e-3)),
    "fb_vcf":    ("fixed", (2e6, 3e7, 3e8)),
    "fb_fasta":  ("fixed", (1e7, 2e8, 5e9)),
    "fb_bed":    ("fixed", (1e5, 5e6, 2e8)),
    "fb_h5":     ("fixed", (1e7, 3e8, 5e9)),
    "fb_office": ("fixed", (1e6, 1e7, 1e8)),
    "fb_table":  ("fixed", (1e5, 3e6, 3e8)),
    "fb_figure": ("fixed", (5e4, 2e6, 5e7)),
    "fb_gz":     ("fixed", (1e6, 5e7, 2e9)),
    "fb_other":  ("fixed", (1e3, 1e5, 1e7)),
}
FALLBACK_BASIS = {cid: model[0] for cid, model in FALLBACK_MODEL.items()}
INDEX_SUFFIXES = (".bai", ".crai", ".pbi", ".tbi", ".csi", ".fai")
INDEX_FIXED = (1e5, 5e7, 5e8)



def glob_to_re(glob_pat):
    """Convert a glob (with cross-directory **) to an anchored regex.

    fnmatch does not support ** across directories, so translate manually.

    Args:
        glob_pat: glob pattern string.

    Returns:
        Compiled regex matching whole relative paths.
    """
    out = []
    i = 0
    while i < len(glob_pat):
        c = glob_pat[i]
        if glob_pat.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif glob_pat.startswith("**", i):
            out.append(".*")
            i += 2
        elif c == "*":
            out.append("[^/]*")
            i += 1
        elif c == "?":
            out.append("[^/]")
            i += 1
        else:
            out.append(re.escape(c))
            i += 1
    return re.compile("^" + "".join(out) + "$")


def human(nbytes):
    """Format a byte count with decimal units (KB/MB/GB...)."""
    x = float(nbytes)
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(x) < 1000.0 or unit == "PB":
            return f"{x:.1f} {unit}" if unit != "B" else f"{x:.0f} B"
        x /= 1000.0


def cfg_get(cfg, dotted, default=None):
    """Fetch cfg.a.b.c by dotted path; return default when any step is absent."""
    obj = cfg
    for part in dotted.split("."):
        if isinstance(obj, dict) and part in obj:
            obj = obj[part]
        else:
            return default
    return obj


def input_kind(path):
    """Classify an input file path as 'input_bam' / 'input_fastq' / None."""
    low = path.lower()
    if low.endswith((".bam", ".cram")):
        return "input_bam"
    if low.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
        return "input_fastq"
    return None



def instance_samples(datajson):
    """Collect sample ids from an instantiated json.

    node.run<Workflow>() stores sample lists under sample_params (new
    layout) and/or Params.workflow (legacy runs); known list keys on both
    levels are unioned.

    Returns:
        Sorted list of sample ids (possibly empty).
    """
    ids = set()

    def add(value):
        if isinstance(value, str):
            ids.add(value)
        elif isinstance(value, list):
            ids.update(v for v in value if isinstance(v, str))
        elif isinstance(value, dict):
            for v in value.values():
                add(v)

    sample_keys = ("samples", "paired_samples", "single_samples",
                   "genome_paired_samples", "genome_single_samples",
                   "dna_paired_samples", "dna_single_samples",
                   "ip_samples", "input_samples", "treated_ip_samples",
                   "treated_input_samples", "control_samples")
    for src in (datajson.get("sample_params") or {},
                cfg_get(datajson, "Params.workflow", {}) or {}):
        for key in sample_keys:
            if key in src:
                add(src[key])
    return sorted(ids)


def instance_raw_files(datajson):
    """Return the run's input file list (sample_params.raw_files or legacy
    Params.workflow.raw_files)."""
    sp = datajson.get("sample_params") or {}
    wf = cfg_get(datajson, "Params.workflow", {}) or {}
    return list(sp.get("raw_files") or wf.get("raw_files") or [])


def instance_input_dirs(datajson):
    """Input DIRECTORIES for runs whose inputs are dirs, not flat file lists.

    scRNAseq in cellranger_dir mode carries Params.workflow.cellranger_dirs
    (or sample_params.cellranger_input_dict fastq_dirs) with raw_files empty.

    Returns:
        List of (dir_path, sample_id_or_None); never mixed with raw_files.
    """
    out, seen = [], set()
    wf = cfg_get(datajson, "Params.workflow", {}) or {}
    for src in (datajson.get("sample_params") or {}, wf):
        cid = src.get("cellranger_input_dict")
        if isinstance(cid, dict):
            for sid, info in cid.items():
                if isinstance(info, dict) and info.get("fastq_dir"):
                    p = info["fastq_dir"]
                    if p not in seen:
                        seen.add(p)
                        out.append((p, sid if isinstance(sid, str) else None))
        for d in (src.get("cellranger_dirs") or []):
            if isinstance(d, str) and d not in seen:
                seen.add(d)
                out.append((d, None))
    return out


def attribute_to_sample(path, samples):
    """Attribute a file path to one of the sample ids.

    Longest id first so 'DMSO_1' wins over 'DMSO'; a match is either a whole
    path component or a boundary-delimited substring of the basename
    (input files are usually named <lane>-<sample>.R1.raw.fastq.gz).

    Returns:
        Sample id or None when no sample id occurs in the path.
    """
    base = os.path.basename(path)
    parts = set(re.split(r"[\\/]", path))
    for sid in sorted(samples, key=len, reverse=True):
        if sid in parts:
            return sid
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(sid)}(?![A-Za-z0-9])", base):
            return sid
    return None


def instance_inputs(datajson, samples, allow_missing=False):
    """Stat the run's input data -> per-sample basis bytes.

    Normally the flat raw_files list; when it is empty (scRNAseq
    cellranger_dir mode) the input DIRECTORIES are walked and every real
    file sized instead. Symlinks are followed (size of the real data).
    Files not attributable to a sample are split evenly across samples and
    reported as a note.

    Args:
        datajson: instantiated run json.
        samples: sample ids from instance_samples().
        allow_missing: when False, missing inputs abort.

    Returns:
        (per_sample, basis_all) where per_sample maps sample_id to
        {"input_bam", "input_fastq", "input_total", "n_files"} and basis_all
        holds the run-wide sums.
    """
    raw_files = instance_raw_files(datajson)
    input_dirs = instance_input_dirs(datajson) if not raw_files else []
    if not raw_files and not input_dirs:
        raise SystemExit("instantiated json lists no input files "
                         "(sample_params.raw_files / Params.workflow.raw_files)")
    per_sample = {}
    totals = {"input_bam": 0, "input_fastq": 0, "input_total": 0}
    unattributed = []
    missing = []

    def add_file(path, size, sid):
        kind = input_kind(path)
        totals["input_total"] += size
        if kind:
            totals[kind] += size
        if sid is None:
            unattributed.append((path, size))
            return
        d = per_sample.setdefault(sid, {"input_bam": 0, "input_fastq": 0,
                                        "input_total": 0, "n_files": 0})
        d["input_total"] += size
        if kind:
            d[kind] += size
        d["n_files"] += 1

    for path in raw_files:
        if not os.path.exists(path):
            missing.append(path)
            continue
        add_file(path, os.path.getsize(path),
                 attribute_to_sample(path, samples) if samples else None)
    for dirpath, owner in input_dirs:
        if not os.path.isdir(dirpath):
            missing.append(dirpath)
            continue
        for sub, _dirs, files in os.walk(dirpath):
            for fn in files:
                p = os.path.join(sub, fn)
                try:
                    size = os.path.getsize(p)
                except OSError:
                    continue
                sid = owner or (attribute_to_sample(p, samples) if samples else None)
                add_file(p, size, sid)
    if missing and not allow_missing:
        raise SystemExit("input files missing (sizes are never guessed), "
                         f"{len(missing)} total, e.g. {missing[0]}")
    if unattributed and samples:
        share = {k: 0 for k in ("input_bam", "input_fastq", "input_total")}
        for path, size in unattributed:
            share["input_total"] += size
            kind = input_kind(path)
            if kind:
                share[kind] += size
        n = len(samples)
        for sid in samples:
            d = per_sample.setdefault(sid, {"input_bam": 0, "input_fastq": 0,
                                            "input_total": 0, "n_files": 0})
            for k in share:
                d[k] += share[k] / n
        print(f"note: {len(unattributed)} input file(s) not attributable to a "
              f"sample id; {human(sum(s for _, s in unattributed))} split evenly "
              f"across {n} samples")
    return per_sample, totals


def instance_comparisons(datajson):
    """Count derived comparison pairs recorded in the instantiated json.

    sample_params.design ({scope: {groups, comparisons}}) is the single
    design model; each comparison expands per coding (consecutive / vs_ref /
    all_pairs / helmert) into two-group contrasts. Legacy flat /
    organism-nested comparison fields are understood too.

    Returns:
        Number of comparison pairs (0 when none).
    """
    design = datajson.get("sample_params", {}).get("design") or {}
    n_design = 0
    for block in design.values():
        for entry in ((block or {}).get("comparisons") or {}).values():
            names = (entry or {}).get("groups") or []
            k = len(names)
            if k < 2:
                continue
            coding = (entry or {}).get("coding") or "consecutive"
            n_design += (k - 1) if coding in ("consecutive", "helmert", "vs_ref") \
                else k * (k - 1) // 2
    if n_design:
        return n_design
    candidates = [datajson.get("sample_params", {}).get("comparisons"),
                  cfg_get(datajson, "Params.sv.comparisons"),
                  cfg_get(datajson, "Params.DESeq2.group_pairs"),
                  cfg_get(datajson, "Params.function.group_pairs")]
    pair_keys = {"control_samples", "experimental_samples",
                 "control_group_name", "experimental_group_name"}
    for cmps in candidates:
        if not isinstance(cmps, dict) or not cmps:
            continue
        n = 0
        for v in cmps.values():
            if isinstance(v, dict) and (set(v) & pair_keys):
                n += 1
            elif isinstance(v, dict):
                n += len(v)
        if n:
            return n
    return 0


def genome_size(datajson):
    """Size of genome.<default>'s fasta (0 when unset or missing)."""
    g = datajson.get("genome") or {}
    name = g.get("default")
    ref = (g.get("references") or {}).get(name) or {}
    fasta = ref.get("fasta")
    if fasta and os.path.exists(fasta):
        return os.path.getsize(fasta)
    return 0


def detect_workflow(datajson, config_dir=None):
    """Best-effort workflow name for an instantiated json.

    raw.json does not record its workflow name, so score the Params +
    resources key signature against every config/<wf>.json template (Dice
    coefficient) and accept a clear winner.

    Args:
        datajson: instantiated run json.
        config_dir: template directory (default <project root>/config).

    Returns:
        Workflow name or None when no template matches confidently.
    """
    config_dir = config_dir or os.path.join(_PROJECT_ROOT, "config")

    def dotted_keys(d, prefix):
        keys = set()
        for k, v in (d or {}).items():
            p = f"{prefix}.{k}" if prefix else k
            keys.add(p)
            if isinstance(v, dict):
                keys |= dotted_keys(v, p)
        return keys

    inst = (dotted_keys(datajson.get("Params"), "Params")
            | dotted_keys(datajson.get("resources"), "res"))
    scores = {}
    for fn in sorted(os.listdir(config_dir)):
        if not fn.endswith(".json") or fn.endswith(".schema.json"):
            continue
        try:
            with open(os.path.join(config_dir, fn)) as fh:
                tpl = json.load(fh)
        except (OSError, ValueError):
            continue
        tk = dotted_keys(tpl.get("Params"), "Params") | dotted_keys(tpl.get("resources"), "res")
        if not tk:
            continue
        scores[fn[:-5]] = 2.0 * len(tk & inst) / (len(tk) + len(inst))
    if not scores:
        return None
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    top_name, top = ranked[0]
    second = ranked[1][1] if len(ranked) > 1 else 0.0
    if top >= 0.2 and top >= 1.25 * max(second, 1e-9):
        return top_name
    return None


def is_instantiated(datajson):
    """A json is instantiated when node.run<Workflow>() registered outfiles
    (the template has none). Inputs may still live in input dirs instead of
    raw_files (scRNAseq cellranger_dir mode)."""
    return bool(datajson.get("outfiles"))



def instantiate_from_meta(meta, workflow, config_path, outdir, comparison=None):
    """Generate the instantiated run json for meta + workflow.

    Args:
        meta: meta tsv/csv path.
        workflow: workflow name (must be in run.WORKFLOW_DISPATCH).
        config_path: uninstantiated config template (config/<wf>.json).
        outdir: run output dir; node.run<Workflow>() writes raw.json here
            and pre-creates log/ and tmp/ exactly like a real run.
        comparison: optional comparison file (--comparison of run.py).

    Returns:
        (datajson, raw_json_path) of the generated instance.
    """
    if _PROJECT_ROOT not in sys.path:
        sys.path.insert(0, _PROJECT_ROOT)
    from src.common.util.MetaUtil import MetadataUtils
    from run import WORKFLOW_DISPATCH

    if workflow not in WORKFLOW_DISPATCH:
        raise SystemExit(f"unknown workflow {workflow!r} "
                         f"(known: {sorted(WORKFLOW_DISPATCH)})")
    with open(config_path) as fh:
        template = json.load(fh)
    os.makedirs(outdir, exist_ok=True)
    meta_util = MetadataUtils(outdir=os.path.abspath(outdir), meta=meta,
                              workflow=workflow, comparison_file=comparison)
    (samples_info_dict, sample_pairs, group_pairs, raw_fastq_dir,
     raw_files, cellranger_input_dict) = meta_util.run()
    _smk, raw_json_path = WORKFLOW_DISPATCH[workflow](
        copy.deepcopy(template), samples_info_dict, sample_pairs, group_pairs,
        raw_fastq_dir, os.path.abspath(outdir), meta, raw_files,
        cellranger_input_dict)
    with open(raw_json_path) as fh:
        datajson = json.load(fh)
    return datajson, raw_json_path



def gate_ok(cls, cfg, n_comparisons, gate_any=None):
    """Evaluate a manifest class's gates against the (instantiated) config.

    Args:
        cls: ClassSpec dict.
        cfg: config dict (instantiated json).
        n_comparisons: number of ctr_/exp_ pairs in the run.
        gate_any: any-of gate conditions (one satisfied is enough).

    Returns:
        True when the class is produced under this config.
    """
    for cond in cls.get("gate", []):
        path, expected, default = cond
        if cfg_get(cfg, path, default) != expected:
            return False
    if gate_any:
        if not any(cfg_get(cfg, p, d) == e for p, e, d in gate_any):
            return False
    if cls.get("need_comparisons") and n_comparisons <= 0:
        return False
    for req in cls.get("require", []):
        if req in ("genome_repeat_bed", "genome_subtel_ref"):
            key = "repeat_bed" if req == "genome_repeat_bed" else "subtel_ref"
            g = cfg.get("genome") or {}
            ref = (g.get("references") or {}).get(g.get("default")) or {}
            if not ref.get(key):
                return False
    return True


def is_temp(cls, cfg):
    """True when the class is temp under this config (temp_when satisfied)."""
    tw = cls.get("temp_when")
    if not tw:
        return False
    path, expected, default = tw
    return cfg_get(cfg, path, default) == expected


def load_history(path):
    """Load the calibration history -> {(workflow, class_id): [ratio, ...]}."""
    hist = {}
    if path and os.path.exists(path):
        with open(path, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                try:
                    hist.setdefault((rec["workflow"], rec["class_id"]), []).append(
                        float(rec["ratio"]))
                except (KeyError, ValueError):
                    continue
    return hist


def calibrated_ratio(hist, workflow, cls_id, manifest_lo, manifest_mid, manifest_hi):
    """Return (lo, mid, hi, n_obs, used) for a class.

    With >= MIN_POINTS real observations use p50 with p10/p90 as the range,
    otherwise the empirical manifest triple.
    """
    obs = hist.get((workflow, cls_id)) or []
    if len(obs) < MIN_POINTS:
        return manifest_lo, manifest_mid, manifest_hi, len(obs), False
    xs = sorted(obs)

    def pct(q):
        k = max(0, min(len(xs) - 1, int(round(q * (len(xs) - 1)))))
        return xs[k]
    return pct(0.10), pct(0.50), pct(0.90), len(obs), True



def generic_fill_classes(datajson, covered_res, per_sample, basis_all, samples):
    """Size every registered outfile not covered by a curated manifest.

    Uses FALLBACK_MODEL's per-file ratios (each file sized against the basis
    of the sample owning it, run totals otherwise). Files are matched to
    FALLBACK_PATTERNS classes so measure/calibrate/estimate share class ids.

    Args:
        datajson: instantiated run json.
        covered_res: list of compiled regexes of included manifest patterns.
        per_sample: per-sample input basis (instance_inputs()).
        basis_all: run-wide input basis.
        samples: sample ids.

    Returns:
        List of class dicts in the estimate's internal shape:
        {id, desc, scope, unit, low, mid, hi, temp, calibrated, n_obs, n_files}
    """
    outdir = cfg_get(datajson, "Params.workflow.outdir") or "."
    outfiles = datajson.get("outfiles") or []
    fallback = [(glob_to_re(pat), cid) for cid, pats in FALLBACK_PATTERNS
                for pat in pats]
    agg = {}
    for out in sorted(set(outfiles)):
        try:
            rel = os.path.relpath(out, outdir)
        except ValueError:
            rel = os.path.basename(out)
        if any(rx.match(rel) for rx in covered_res):
            continue
        cid = "fb_other"
        for rx, cand in fallback:
            if rx.match(rel):
                cid = cand
                break
        basis_key, ratio = FALLBACK_MODEL[cid]
        if out.lower().endswith(INDEX_SUFFIXES):
            lo, mid, hi = INDEX_FIXED
        else:
            sid = attribute_to_sample(out, samples) if samples else None
            if basis_key == "fixed":
                lo, mid, hi = ratio
            else:
                unit = 0.0
                if sid is not None:
                    unit = per_sample.get(sid, {}).get(basis_key, 0) or \
                        per_sample.get(sid, {}).get("input_total", 0)
                else:
                    unit = basis_all.get(basis_key, 0) or basis_all.get("input_total", 0)
                lo, mid, hi = (ratio[0] * unit, ratio[1] * unit, ratio[2] * unit)
        d = agg.setdefault(cid, {"lo": 0.0, "mid": 0.0, "hi": 0.0, "files": 0})
        d["lo"] += lo
        d["mid"] += mid
        d["hi"] += hi
        d["files"] += 1
    rows = []
    for cid, d in sorted(agg.items()):
        rows.append(dict(
            id=cid, desc=f"registered outfiles matched by suffix ({cid})",
            scope="run", unit=f"{d['files']} outfiles x per-file model",
            low=d["lo"], mid=d["mid"], hi=d["hi"], temp=False,
            calibrated=False, n_obs=0, n_files=d["files"]))
    return rows


def do_estimate(args):
    """Run the pre-run storage estimate.

    Reads an instantiated json (--run-json) directly, or builds one from
    --meta + --workflow via instantiate_from_meta(); then sizes curated
    manifest classes and the remaining registered outfiles.

    Returns:
        Exit code.
    """
    if args.run_json:
        with open(args.run_json) as fh:
            datajson = json.load(fh)
        if not is_instantiated(datajson):
            raise SystemExit(f"{args.run_json} is not an instantiated run json "
                             "(no outfiles/raw_files); pass it to --config with "
                             "--meta/--workflow instead")
        if args.meta:
            print(f"note: --meta ignored, inputs come from {args.run_json}")
        workflow = args.workflow or detect_workflow(datajson)
        instance_path = args.run_json
    else:
        if not args.meta or not args.workflow:
            raise SystemExit("need either --run-json (instantiated json) or "
                             "--meta + --workflow (template is read from "
                             "config/<workflow>.json and instantiated)")
        workflow = args.workflow
        config_path = args.config or os.path.join(_PROJECT_ROOT, "config",
                                                  f"{workflow}.json")
        if not os.path.exists(config_path):
            raise SystemExit(f"config template not found: {config_path}")
        if args.outdir:
            outdir = os.path.abspath(args.outdir)
        else:
            outdir = tempfile.mkdtemp(prefix="estimate_storage_")
        datajson, instance_path = instantiate_from_meta(
            args.meta, workflow, config_path, outdir, comparison=args.comparison)
        print(f"instance: generated {instance_path}"
              + ("" if args.outdir else " (scratch dir; use --outdir to keep it "
                                       "as the real run dir)"))
    if args.save_json:
        shutil.copyfile(instance_path, args.save_json)
        print(f"instance copy: {args.save_json}")

    if args.manifest:
        with open(args.manifest) as fh:
            MANIFESTS.update(json.load(fh))

    samples = instance_samples(datajson)
    per_sample, basis_all = instance_inputs(datajson, samples,
                                            allow_missing=args.allow_missing)
    if not samples:
        samples = sorted(per_sample)
    if not samples:
        print("note: no sample ids found; sample-scoped classes estimated "
              "as a single unit")
        samples = ["__run__"]
        per_sample = {"__run__": dict(basis_all, n_files=len(instance_raw_files(datajson)))}
    n_cmp = instance_comparisons(datajson)
    gsize = genome_size(datajson)
    hist = load_history(args.calibration or
                        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     HISTORY_NAME))

    curated = MANIFESTS.get(workflow) if workflow else None
    gate_any = GATE_ANY.get(workflow) if workflow else None
    manifest_ids = set()
    table = []
    final = [0.0, 0.0, 0.0]
    transient = [0.0, 0.0, 0.0]
    timeline = []
    skipped = []
    covered_res = []

    def fold(lo, mid, hi, temp):
        for i, v in enumerate((lo, mid, hi)):
            (transient if temp else final)[i] += v

    for cls in (curated or []):
        if not gate_ok(cls, datajson, n_cmp,
                       gate_any=gate_any.get(cls["id"]) if gate_any else None):
            skipped.append((cls["id"], cls["desc"]))
            continue
        r_lo, r_mid, r_hi = cls["ratio"]
        r_lo, r_mid, r_hi, n_obs, used = calibrated_ratio(
            hist, workflow, cls["id"], r_lo, r_mid, r_hi)
        basis_key = cls["basis"] if cls["basis"] != "fixed" else None

        def sid_basis(sid, basis_key=basis_key):
            if not basis_key:
                return 0.0
            if basis_key == "genome_fasta":
                return gsize
            return per_sample.get(sid, {}).get(basis_key, 0)

        if cls["scope"] == "sample":
            if basis_key:
                lo = sum(r_lo * sid_basis(s) for s in samples)
                mid = sum(r_mid * sid_basis(s) for s in samples)
                hi = sum(r_hi * sid_basis(s) for s in samples)
                unit_desc = f"{len(samples)} samples x input ratio"
            else:
                lo, mid, hi = (r_lo * len(samples), r_mid * len(samples),
                               r_hi * len(samples))
                unit_desc = f"{len(samples)} samples x fixed"
        else:
            if basis_key:
                b = gsize if basis_key == "genome_fasta" else basis_all.get(basis_key, 0)
                lo, mid, hi = r_lo * b, r_mid * b, r_hi * b
                unit_desc = "run x input ratio"
            else:
                lo, mid, hi = r_lo, r_mid, r_hi
                unit_desc = "run fixed"
        temp = is_temp(cls, datajson)
        fold(lo, mid, hi, temp)
        born = cls.get("born", 1)
        die = cls.get("die") if temp else None
        timeline.append((born, die, lo, mid, hi))
        manifest_ids.add(cls["id"])
        covered_res.extend(glob_to_re(p) for p in cls["patterns"])
        table.append(dict(
            id=cls["id"], desc=cls["desc"], scope=cls["scope"], unit=unit_desc,
            low=lo, mid=mid, hi=hi, temp=temp,
            calibrated=used, n_obs=n_obs))

    generic_rows = generic_fill_classes(datajson, covered_res, per_sample,
                                        basis_all, samples)
    for row in generic_rows:
        obs = hist.get((workflow, row["id"])) or [] if workflow else []
        if len(obs) >= MIN_POINTS:
            xs = sorted(obs)
            pick = lambda q: xs[max(0, min(len(xs) - 1, int(round(q * (len(xs) - 1)))))]
            if FALLBACK_BASIS.get(row["id"], "fixed") == "fixed":
                row["low"], row["mid"], row["hi"] = pick(0.10), pick(0.50), pick(0.90)
            else:
                b = basis_all.get(FALLBACK_BASIS[row["id"]], 0) or basis_all.get("input_total", 0)
                row["low"], row["mid"], row["hi"] = pick(0.10) * b, pick(0.50) * b, pick(0.90) * b
            row["calibrated"], row["n_obs"] = True, len(obs)
        fold(row["low"], row["mid"], row["hi"], row["temp"])
        timeline.append((GENERIC_BORN, None, row["low"], row["mid"], row["hi"]))
        table.append(row)

    if not curated:
        r_lo, r_mid, r_hi, n_obs, used = calibrated_ratio(
            hist, workflow or "", "run_overhead", 5e8, 2e9, 1e10)
        fold(r_lo, r_mid, r_hi, False)
        timeline.append((1, None, r_lo, r_mid, r_hi))
        table.append(dict(id="run_overhead",
                          desc="raw.json/log/tmp/.snakemake/input-symlink dirs fixed overhead",
                          scope="run", unit="run fixed",
                          low=r_lo, mid=r_mid, hi=r_hi, temp=False,
                          calibrated=used, n_obs=n_obs))

    checkpoints = set()
    for b, d, _lo, _mid, _hi in timeline:
        checkpoints.add(b)
        if d is not None:
            checkpoints.add(d)
            checkpoints.add(d + 1)
    peak_stage, peak = None, [0.0, 0.0, 0.0]
    for s in sorted(checkpoints):
        acc = [0.0, 0.0, 0.0]
        for b, d, lo_i, mid_i, hi_i in timeline:
            if b <= s and (d is None or s <= d):
                acc[0] += lo_i
                acc[1] += mid_i
                acc[2] += hi_i
        if acc[1] > peak[1]:
            peak_stage, peak = s, acc

    print(f"workflow: {workflow or '(not detected)'}   samples: {len(samples)}   "
          f"comparisons: {n_cmp}   genome fasta: {human(gsize)}")
    print(f"inputs: bam {human(basis_all['input_bam'])}, "
          f"fastq {human(basis_all['input_fastq'])}, "
          f"total {human(basis_all['input_total'])}")
    print()
    hdr = f"{'class':<22}{'scope':<8}{'kept/temp':<10}{'low':>10}{'mid':>10}{'high':>10}  note"
    print(hdr)
    print("-" * len(hdr))
    for r in table:
        note = (f"calibrated n={r['n_obs']}" if r["calibrated"]
                else f"empirical ratio n={r['n_obs']}")
        if r.get("n_files"):
            note += f", {r['n_files']} files"
        print(f"{r['id']:<22}{r['scope']:<8}{'temp' if r['temp'] else 'kept':<10}"
              f"{human(r['low']):>10}{human(r['mid']):>10}{human(r['hi']):>10}  {note}")
    if skipped:
        print()
        print("not produced (disabled in config or missing prerequisite):")
        for cid, desc in skipped:
            print(f"  - {cid}: {desc}")
    print()
    print("=" * 60)
    print(f"final kept (long-term after the run): {human(final[1])}  "
          f"[{human(final[0])} ~ {human(final[2])}]")
    print(f"transient (temp, auto-deleted later): {human(transient[1])}  "
          f"[{human(transient[0])} ~ {human(transient[2])}]")
    print(f"peak (stage {peak_stage}, reserve this much): {human(peak[1])}  "
          f"[{human(peak[0])} ~ {human(peak[2])}]")
    print(f"                    (decimal vs binary: {final[1]/1e9:.1f} GB = "
          f"{final[1]/2**30:.1f} GiB)")
    print("=" * 60)
    if not curated:
        print("note: no curated manifest for this workflow - temp "
              "intermediates are not modeled, so the real peak can exceed "
              "the estimate (registered outfiles + overhead only).")
    print("not modeled: GATK java tmp-dir churn (out/tmp, cleaned after the "
          "run) and filesystem block overhead.")

    if args.report_json:
        report = dict(
            workflow=workflow, generated=datetime.datetime.now().isoformat(
                timespec="seconds"),
            instance_json=instance_path, meta=args.meta, run_json=args.run_json,
            curated_manifest=bool(curated),
            n_samples=len(samples), n_comparisons=n_cmp,
            input_bytes=basis_all, genome_fasta_bytes=gsize,
            classes=table, skipped=[c for c, _ in skipped],
            final_bytes=dict(zip(("low", "mid", "high"), final)),
            transient_bytes=dict(zip(("low", "mid", "high"), transient)),
            peak_bytes=dict(zip(("low", "mid", "high"), peak)),
            peak_stage=peak_stage,
            unmodeled=["gatk_java_tmpdir", "fs_block_overhead"])
        with open(args.report_json, "w") as fh:
            json.dump(report, fh, indent=2, ensure_ascii=False)
        print(f"report json: {args.report_json}")
    return 0



def classify_matcher(workflow=None):
    """Build (relpath) -> class_id: manifest patterns first, suffix fallback.

    Args:
        workflow: optional workflow name for manifest-precision classes.

    Returns:
        Matcher function.
    """
    entries = []
    if workflow and workflow in MANIFESTS:
        for cls in MANIFESTS[workflow]:
            for pat in cls["patterns"]:
                entries.append((glob_to_re(pat), cls["id"]))
    for cid, pats in FALLBACK_PATTERNS:
        for pat in pats:
            entries.append((glob_to_re(pat), cid))

    def match(relpath):
        for rx, cid in entries:
            if rx.match(relpath):
                return cid
        return "fb_other"
    return match


def walk_root(root, matcher):
    """Scan root -> {class_id: {"bytes": n, "files": n}} + symlink count.

    Symlinks (input data links) are counted but not sized, so real data is
    not double-counted.
    """
    agg = {}
    n_links = 0
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            if os.path.islink(p):
                n_links += 1
                continue
            try:
                size = os.path.getsize(p)
            except OSError:
                continue
            rel = os.path.relpath(p, root)
            cid = matcher(rel)
            d = agg.setdefault(cid, {"bytes": 0, "files": 0})
            d["bytes"] += size
            d["files"] += 1
    return agg, n_links


def do_measure(args):
    """Measure the actual disk usage of a run dir, classified by product class."""
    matcher = classify_matcher(args.workflow)
    agg, n_links = walk_root(args.root, matcher)
    total = sum(d["bytes"] for d in agg.values())
    print(f"root: {args.root}")
    hdr = f"{'class':<22}{'files':>8}{'bytes':>14}{'size':>10}{'pct':>7}"
    print(hdr)
    print("-" * len(hdr))
    for cid, d in sorted(agg.items(), key=lambda kv: -kv[1]["bytes"]):
        pct = 100.0 * d["bytes"] / total if total else 0
        print(f"{cid:<22}{d['files']:>8}{d['bytes']:>14}{human(d['bytes']):>10}{pct:>6.1f}%")
    print("-" * len(hdr))
    print(f"{'TOTAL':<22}{sum(d['files'] for d in agg.values()):>8}{total:>14}"
          f"{human(total):>10}{100.0:>6.1f}%")
    if n_links:
        print(f"({n_links} symlinks skipped - input data links are not "
              "double-counted)")
    if args.report_json:
        with open(args.report_json, "w") as fh:
            json.dump(dict(root=args.root, generated=datetime.datetime.now().isoformat(
                timespec="seconds"), total_bytes=total, symlinks_skipped=n_links,
                classes=agg), fh, indent=2, ensure_ascii=False)
        print(f"report json: {args.report_json}")
    return 0


def parse_meta_rows(path):
    """Parse a meta tsv/csv into row dicts (delimiter by extension, sniffed).

    Kept as the calibrate() input path; estimate() goes through
    MetadataUtils via instantiate_from_meta() instead.

    Returns:
        List of row dicts (needs a sample_id column).
    """
    delim = "," if path.lower().endswith(".csv") else "\t"
    with open(path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        first_line = sample.splitlines()[0] if sample else ""
        if delim == "\t" and "\t" not in first_line and "," in first_line:
            delim = ","
        rows = list(csv.DictReader(fh, delimiter=delim))
    if not rows or "sample_id" not in rows[0]:
        raise SystemExit(f"meta parse failed (missing sample_id column or empty): {path}")
    return rows


INPUT_FILE_COLS = ("fastq_1", "fastq_2", "bam", "pbi", "bai", "ms_file")


def sample_inputs_from_meta(rows, allow_missing=False):
    """meta rows -> {sample_id: basis bytes}, real stat on input files.

    Multiple rows per sample (data_id one-to-many) are summed. Missing files
    abort unless allow_missing (sizes are never guessed).

    Returns:
        (per_sample, basis_all) with keys input_bam/input_fastq/input_total.
    """
    per = {}
    missing = []
    for row in rows:
        sid = (row.get("sample_id") or "").strip()
        if not sid:
            continue
        d = per.setdefault(sid, {"input_bam": 0, "input_fastq": 0,
                                 "input_total": 0, "n_files": 0})
        for col in INPUT_FILE_COLS:
            p = (row.get(col) or "").strip()
            if not p:
                continue
            if not os.path.exists(p):
                missing.append(p)
                continue
            size = os.path.getsize(p)
            d["input_total"] += size
            kind = input_kind(p)
            if kind:
                d[kind] += size
            d["n_files"] += 1
    if missing and not allow_missing:
        raise SystemExit("meta input files missing (sizes are never guessed), "
                         f"{len(missing)} total, e.g. {missing[0]}")
    basis_all = {k: sum(d[k] for d in per.values())
                 for k in ("input_bam", "input_fastq", "input_total")}
    return per, basis_all


def do_calibrate(args):
    """Measure a finished run and append per-class observed/basis ratios to
    the append-only history so estimate() gets more accurate over time.

    Inputs come from --run-json (or <root>/raw.json when present) or from
    --meta (+ optional --config for the genome fasta basis).

    Returns:
        Exit code.
    """
    run_json = args.run_json
    if not run_json and os.path.isfile(os.path.join(args.root, "raw.json")):
        run_json = os.path.join(args.root, "raw.json")

    workflow = args.workflow
    gsize = 0
    if run_json:
        with open(run_json) as fh:
            datajson = json.load(fh)
        workflow = workflow or detect_workflow(datajson)
        samples = instance_samples(datajson)
        per_sample, basis_all = instance_inputs(datajson, samples,
                                                allow_missing=args.allow_missing)
        if not samples:
            samples = sorted(per_sample)
        gsize = genome_size(datajson)
    else:
        if not args.meta:
            raise SystemExit("need --run-json (or <root>/raw.json) or --meta "
                             "to derive per-sample input baselines")
        rows = parse_meta_rows(args.meta)
        per_sample, basis_all = sample_inputs_from_meta(
            rows, allow_missing=args.allow_missing)
        samples = sorted(per_sample)
        if args.config:
            with open(args.config) as fh:
                gsize = genome_size(json.load(fh))
    if not workflow:
        raise SystemExit("workflow name unknown; pass --workflow (needed as "
                         "the calibration history key)")
    if workflow not in MANIFESTS and not args.manifest:
        raise SystemExit(f"no product manifest for {workflow!r}; cannot map "
                         "files to calibration classes (use --manifest)")
    if args.manifest:
        with open(args.manifest) as fh:
            MANIFESTS.update(json.load(fh))

    matcher = classify_matcher(workflow)
    cls_by_id = {c["id"]: c for c in MANIFESTS[workflow]}
    history = args.history or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), HISTORY_NAME)

    per_class_sample = {}
    per_class_run = {}
    n_links = 0
    for dirpath, _dirnames, filenames in os.walk(args.root):
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            if os.path.islink(p):
                n_links += 1
                continue
            try:
                size = os.path.getsize(p)
            except OSError:
                continue
            rel = os.path.relpath(p, args.root)
            cid = matcher(rel)
            sid = None
            for part in rel.split(os.sep):
                if part in per_sample:
                    sid = part
                    break
            cls = cls_by_id.get(cid)
            if cls is not None and cls["scope"] == "sample" and sid:
                per_class_sample[(cid, sid)] = per_class_sample.get((cid, sid), 0) + size
            else:
                per_class_run[cid] = per_class_run.get(cid, 0) + size

    run_id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    ts = datetime.datetime.now().isoformat(timespec="seconds")
    new_rows = 0
    exists = os.path.exists(history)
    with open(history, "a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        if not exists:
            w.writerow(["run_id", "ts", "workflow", "class_id", "scope",
                        "sample_id", "basis_bytes", "observed_bytes", "ratio", "source"])

        def emit(cid, scope, sid, basis_b, obs_b, cls):
            nonlocal new_rows
            if obs_b <= 0:
                return
            fixed = (cls is not None and cls["basis"] == "fixed") or \
                (cls is None and FALLBACK_BASIS.get(cid) == "fixed")
            if fixed:
                ratio = obs_b
                basis_b = ""
            else:
                if basis_b <= 0:
                    return
                ratio = obs_b / basis_b
            w.writerow([run_id, ts, workflow, cid, scope, sid or "",
                        f"{basis_b:.0f}" if basis_b != "" else "",
                        f"{obs_b:.0f}", f"{ratio:.6g}", args.root])
            new_rows += 1

        for (cid, sid), obs in sorted(per_class_sample.items()):
            cls = cls_by_id.get(cid)
            basis_key = cls["basis"] if cls and cls["basis"] != "fixed" else None
            basis_b = per_sample[sid].get(basis_key, 0) if basis_key else 0
            emit(cid, "sample", sid, basis_b, obs, cls)
        for cid, obs in sorted(per_class_run.items()):
            cls = cls_by_id.get(cid)
            if cls is not None and cls["scope"] == "sample":
                continue
            basis_key = None
            if cls is not None and cls["basis"] != "fixed":
                basis_key = cls["basis"]
            elif cls is None:
                fb = FALLBACK_BASIS.get(cid)
                if fb and fb != "fixed":
                    basis_key = fb
            basis_b = basis_all.get(basis_key, 0) if basis_key else 0
            emit(cid, "run", None, basis_b, obs, cls)
    print(f"calibrate: appended {new_rows} rows -> {history}")
    print(f"(estimate reads this table once a class has >= {MIN_POINTS} "
          "points; until then the empirical ratios apply)")
    return 0



def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("estimate", help="estimate disk usage before running")
    e.add_argument("-j", "--run-json", default=None, metavar="JSON",
                   help="instantiated run json (raw.json): inputs, samples and "
                        "module switches all come from it; meta not needed")
    e.add_argument("-m", "--meta", default=None, metavar="TSV",
                   help="meta tsv/csv (uninstantiated mode; combined with "
                        "--workflow via config/<workflow>.json)")
    e.add_argument("-w", "--workflow", default=None, metavar="NAME",
                   help="workflow name (required with --meta; optional "
                        "override with --run-json when auto-detection is wrong)")
    e.add_argument("-c", "--config", default=None, metavar="JSON",
                   help="uninstantiated config template (default "
                        "config/<workflow>.json)")
    e.add_argument("-o", "--outdir", default=None, metavar="DIR",
                   help="run output dir to instantiate into (default: a "
                        "scratch dir; pass the real dir to keep raw.json)")
    e.add_argument("-s", "--save-json", default=None, metavar="JSON",
                   help="write a copy of the instantiated json here")
    e.add_argument("-p", "--comparison", default=None, metavar="FILE",
                   help="comparison file (run.py --comparison), used during "
                        "instantiation")
    e.add_argument("-x", "--manifest", default=None, metavar="JSON",
                   help="extra product-class manifest merged into MANIFESTS")
    e.add_argument("-t", "--calibration", default=None, metavar="TSV",
                   help=f"calibration history (default script dir {HISTORY_NAME})")
    e.add_argument("-r", "--report-json", default=None, metavar="JSON",
                   help="also write the report as JSON")
    e.add_argument("-a", "--allow-missing", action="store_true",
                   help="do not fail on missing input files (counted as 0, "
                        "underestimates)")
    e.set_defaults(func=do_estimate)

    m = sub.add_parser("measure", help="measure an existing run dir (any workflow)")
    m.add_argument("-R", "--root", default=".", metavar="DIR",
                   help="run output root (default: .)")
    m.add_argument("-w", "--workflow", default=None, metavar="NAME",
                   help="classify with this workflow's manifest (optional)")
    m.add_argument("-r", "--report-json", default=None, metavar="JSON",
                   help="also write the report as JSON")
    m.set_defaults(func=do_measure)

    c = sub.add_parser("calibrate",
                       help="measure a real run and append ratios to the history")
    c.add_argument("-R", "--root", required=True, metavar="DIR",
                   help="run output root")
    c.add_argument("-j", "--run-json", default=None, metavar="JSON",
                   help="instantiated run json for input baselines "
                        "(default: <root>/raw.json when present)")
    c.add_argument("-m", "--meta", default=None, metavar="TSV",
                   help="meta tsv of this run (alternative to --run-json)")
    c.add_argument("-w", "--workflow", default=None, metavar="NAME",
                   help="workflow name (history key; auto-detected from "
                        "--run-json when omitted)")
    c.add_argument("-c", "--config", default=None, metavar="JSON",
                   help="config json (genome fasta baseline, --meta mode)")
    c.add_argument("-x", "--manifest", default=None, metavar="JSON",
                   help="extra product-class manifest")
    c.add_argument("-H", "--history", default=None, metavar="TSV",
                   help=f"calibration history (default script dir {HISTORY_NAME}, "
                        "append-only)")
    c.add_argument("-a", "--allow-missing", action="store_true",
                   help="do not fail on missing input files")
    c.set_defaults(func=do_calibrate)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
