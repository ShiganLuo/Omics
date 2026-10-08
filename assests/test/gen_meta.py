#!/usr/bin/env python3
"""Generate meta_*.tsv for a subworkflow by scanning a directory for sequencing files.

Usage:
    python gen_meta.py -w RNAseq -d /path/to/data_dir [-o meta_RNAseq.tsv]

Column patterns follow assests/test/meta_*.tsv; input file kinds follow
src/common/util/MetaUtil.py (fastq / bam+pbi / ms_file):
  - fastq   : sample_id/data_id/design/fastq_1/fastq_2 (fastq_2 empty for SE)
  - bam_pbi : sample_id/data_id/design/bam/pbi (Fiberseq, PacBio)
  - ms      : sample_id/data_id/design/ms_file (QuantMS, .raw/.mzML/.mgf)

Rules:
  - data_id identifies one unique sequencing file (read pair / single file / bam / ms file)
  - sample_id defaults to data_id (override via --sample-map to merge replicates)
  - only relevant columns are emitted; optional fields stay empty unless given
  - paths are written absolute (meta_tsv_filling.md pitfall #4)
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
try:
    from src.common.util.LogUtil import setup_logger
    logger = setup_logger(__name__)
except Exception:
    import logging
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    logger = logging.getLogger(__name__)

FASTQ_COLS = ["sample_id", "data_id", "design", "fastq_1", "fastq_2"]

# kind maps to MetaUtil.prepare_fastq_meta / prepare_pacbio_meta / prepare_ms_meta;
# columns match assests/test/meta_*.tsv
WORKFLOWS = {
    "CoCulture":       {"kind": "fastq", "columns": FASTQ_COLS},
    "KARRseq":         {"kind": "fastq", "columns": FASTQ_COLS},
    "MERIP":           {"kind": "fastq", "columns": FASTQ_COLS},
    "Mutation":        {"kind": "fastq", "columns": FASTQ_COLS},
    "ncRNAseq":        {"kind": "fastq", "columns": FASTQ_COLS},
    "PeakCalling":     {"kind": "fastq", "columns": FASTQ_COLS},
    "tRNAseq":         {"kind": "fastq", "columns": FASTQ_COLS},
    "CLIP":            {"kind": "fastq", "columns": FASTQ_COLS},
    "PacVar":          {"kind": "fastq", "columns": FASTQ_COLS},
    "scRNAseq":        {"kind": "fastq", "columns": FASTQ_COLS},
    "LRtranscriptome": {"kind": "fastq", "columns": FASTQ_COLS + ["organism"]},
    "RNAseq":          {"kind": "fastq", "columns": ["sample_id", "organism", "group", "data_id",
                                                     "design", "fastq_1", "fastq_2", "contaminated_organism"]},
    "Fiberseq":        {"kind": "bam_pbi", "columns": ["sample_id", "data_id", "design", "bam", "pbi"]},
    "QuantMS":         {"kind": "ms", "columns": ["sample_id", "data_id", "design", "ms_file"]},
}

FQ_PE_RE = re.compile(r"^(?P<base>.+?)(?:_R?(?P<read>[12]))[^/]*\.(?:f(?:ast)?q)(?:\.gz)?$", re.IGNORECASE)
FQ_EXT = (".fq.gz", ".fastq.gz", ".fq", ".fastq")
MS_EXT = (".raw", ".mzml", ".mgf", ".wiff", ".d")
REP_SUFFIX_RE = re.compile(r"(?:[-_](?:Rep)?\d+)$", re.IGNORECASE)


def die(msg):
    logger.error(msg)
    sys.exit(1)


def find_fastq_units(search_dir):
    """Collect (data_id, fastq_1, fastq_2) units; pair _1/_2 or _R1/_R2, else SE."""
    groups = defaultdict(lambda: {"1": [], "2": [], "se": []})
    for f in sorted(search_dir.rglob("*")):
        if not f.is_file() or not f.name.lower().endswith(FQ_EXT):
            continue
        m = FQ_PE_RE.match(f.name)
        if m:
            groups[m.group("base")][m.group("read")].append(f)
        else:
            stem = re.sub(r"\.(?:f(?:ast)?q)(?:\.gz)?$", "", f.name, flags=re.IGNORECASE)
            groups[stem]["se"].append(f)

    units = []
    for base in sorted(groups):
        g = groups[base]
        r1s, r2s, ses = sorted(g["1"]), sorted(g["2"]), sorted(g["se"])
        if r1s and r2s:
            if len(r1s) != len(r2s):
                die(f"{base}: R1 count ({len(r1s)}) != R2 count ({len(r2s)}), check file naming")
            multi = len(r1s) > 1
            for i, (r1, r2) in enumerate(zip(r1s, r2s), 1):
                units.append((f"{base}_L{i:03d}" if multi else base, r1, r2))
        elif r1s or r2s:
            if not r1s:
                logger.warning(f"{base}: only R2 files found, treating as single-end")
            files = sorted(r1s or r2s)
            for i, fq in enumerate(files, 1):
                units.append((f"{base}_L{i:03d}" if len(files) > 1 else base, fq, None))
        else:
            for i, se in enumerate(ses, 1):
                units.append((f"{base}_L{i:03d}" if len(ses) > 1 else base, se, None))
    return units


def find_bam_pbi_units(search_dir):
    """Collect (data_id, bam, pbi) units; pbi is required by MetaUtil pacbio handling."""
    units = []
    for bam in sorted(search_dir.rglob("*.bam")):
        if not bam.is_file():
            continue
        pbi = next((c for c in (Path(str(bam) + ".pbi"), bam.with_suffix(".pbi")) if c.is_file()), None)
        if pbi is None:
            logger.warning(f"{bam}: no .pbi index found, skipping (required for Fiberseq)")
            continue
        units.append((bam.name[:-4], bam, pbi))
    return units


def find_ms_units(search_dir):
    """Collect (data_id, ms_file) units from MS raw files."""
    return [(f.stem, f, None)
            for f in sorted(search_dir.rglob("*"))
            if f.is_file() and f.name.lower().endswith(MS_EXT)]


def load_map(path):
    """Load a two-column TSV (key<TAB>value) into a dict."""
    m = {}
    if path:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line and not line.startswith("#"):
                    k, _, v = line.partition("\t")
                    m[k.strip()] = v.strip()
    return m


def main():
    ap = argparse.ArgumentParser(description="Generate metadata TSV for a subworkflow")
    ap.add_argument("-w", "--workflow", required=True, choices=sorted(WORKFLOWS))
    ap.add_argument("-d", "--search-dir", required=True, nargs="+",
                    help="one or more directories to scan recursively")
    ap.add_argument("-o", "--out", default=None, help="output TSV (default meta_<workflow>.tsv)")
    ap.add_argument("--kind", choices=["fastq", "bam_pbi", "ms"], default=None,
                    help="override default input kind (e.g. --kind fastq for QuantMS test data)")
    ap.add_argument("--layout", choices=["any", "pe", "se"], default="any",
                    help="fastq layout filter (any = auto-detect)")
    ap.add_argument("--organism", default="", help="fill organism column (workflows that have it)")
    ap.add_argument("--contaminated-organism", default="", help="fill contaminated_organism (RNAseq)")
    ap.add_argument("--group", default="", help="fill group column (RNAseq)")
    ap.add_argument("--sample-map", default=None, help="TSV: data_id<TAB>sample_id (merge replicates)")
    ap.add_argument("--design-map", default=None, help="TSV: data_id or sample_id<TAB>design")
    ap.add_argument("--design-strip-rep", action="store_true",
                    help="fill design from sample_id with replicate suffix stripped")
    args = ap.parse_args()

    spec = WORKFLOWS[args.workflow]
    kind = args.kind or spec["kind"]
    file_cols = {"fastq": ["fastq_1", "fastq_2"], "bam_pbi": ["bam", "pbi"], "ms": ["ms_file"]}[kind]
    columns = spec["columns"]
    if not set(file_cols) <= set(columns):
        # kind overridden: rebuild columns around the file columns of that kind
        base = ["sample_id", "data_id", "design"]
        extras = [c for c in columns if c not in base + ["fastq_1", "fastq_2", "bam", "pbi", "ms_file"]]
        columns = base + file_cols + extras

    units = []
    find = {"fastq": find_fastq_units, "bam_pbi": find_bam_pbi_units, "ms": find_ms_units}[kind]
    for d in args.search_dir:
        search_dir = Path(d).resolve()
        if not search_dir.is_dir():
            die(f"search directory not found: {search_dir}")
        units.extend(find(search_dir))
    if kind == "fastq" and args.layout != "any":
        want_pe = args.layout == "pe"
        units = [u for u in units if (u[2] is not None) == want_pe]
    if not units:
        die(f"no {kind} sequencing files found under {', '.join(args.search_dir)}")

    seen = set()
    for did, _, _ in units:
        if did in seen:
            die(f"duplicate data_id: {did}")
        seen.add(did)

    sample_map = load_map(args.sample_map)
    design_map = load_map(args.design_map)

    out = Path(args.out) if args.out else Path(f"meta_{args.workflow}.tsv")
    n = 0
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for did, f1, f2 in units:
            sid = sample_map.get(did, did)  # sample_id defaults to data_id
            if re.search(r"\s", sid):
                die(f"sample_id '{sid}' contains whitespace (MetaUtil rejects it)")
            design = design_map.get(did) or design_map.get(sid) or ""
            if not design and args.design_strip_rep:
                design = REP_SUFFIX_RE.sub("", sid)
            row = {"sample_id": sid, "data_id": did, "design": design}
            if kind == "fastq":
                row["fastq_1"] = os.path.abspath(f1)
                row["fastq_2"] = os.path.abspath(f2) if f2 else ""  # SE: leave empty, never "NA"
            elif kind == "bam_pbi":
                row["bam"], row["pbi"] = os.path.abspath(f1), os.path.abspath(f2)
            else:
                row["ms_file"] = os.path.abspath(f1)
            row["organism"] = args.organism
            row["group"] = args.group
            row["contaminated_organism"] = args.contaminated_organism
            w.writerow([row.get(c, "") for c in columns])
            n += 1

    logger.info(f"workflow={args.workflow} kind={kind} -> {out} ({n} samples, columns: {','.join(columns)})")


if __name__ == "__main__":
    main()
