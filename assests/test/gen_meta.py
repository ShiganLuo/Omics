#!/usr/bin/env python3
"""Generate meta_*.tsv for a subworkflow by scanning a directory for sequencing files.

Usage:
    python gen_meta.py -w RNAseq -d /path/to/data_dir [-o meta_RNAseq.tsv]

Column patterns follow assests/test/meta_*.tsv; input file kinds follow
src/common/util/MetaUtil.py (fastq / bam+pbi / ms_file):
  - fastq   : sample_id/data_id/design/fastq_1/fastq_2 (fastq_2 empty for SE)
  - bam_pbi : sample_id/data_id/design/bam/pbi (Fiberseq, PacBio)
  - ms      : sample_id/data_id/design/ms_file (QuantMS, .raw/.mzML/.mgf)
  - data_type column (DR): DNA/RNA, filled via --data-type / --data-type-map
    and routing raw staging into {outdir}/{DNA|RNA}/common/1_raw_data
  - organism column: genome version / species alias (e.g. GRCm39, mouse);
    PacVar/Fiberseq default to the mouse T2T assembly T2T_mhaESC_v1.5 when
    --organism is empty (WORKFLOW_DEFAULT_ORGANISM)
  - assay column (Mutation, DR): WGS / Exome / Panel — node.py derives CNVkit
    --seq-method and the exome capture targets BED from it; DR fills it on DNA
    rows only (RNA rows empty) via --assay / --assay-map

Rules:
  - data_id identifies one unique sequencing file (read pair / single file / bam / ms file)
  - sample_id defaults to data_id (override via --sample-map to merge replicates)
  - -d may be repeated; every directory is scanned (-d A -d B == -d A B);
    a search dir may be tagged with its data_type: -d DIR=DNA / -d DIR=RNA
  - duplicate-name detection is scoped per data_type: the same name under two
    data_types is ONE sample assayed twice (a DRclass pairing) and passes
    through untouched; within one data_type every colliding member gets a
    numbered suffix in scan order, first -1 then -2/-3/...
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

WORKFLOWS = {
    "CoCulture":       {"kind": "fastq", "columns": FASTQ_COLS},
    "DR":              {"kind": "fastq", "columns": ["sample_id", "organism", "group", "data_id",
                                                     "design", "data_type", "assay", "fastq_1", "fastq_2",
                                                     "contaminated_organism", "strain"]},
    "KARRseq":         {"kind": "fastq", "columns": FASTQ_COLS},
    "MERIP":           {"kind": "fastq", "columns": FASTQ_COLS},
    "Mutation":        {"kind": "fastq", "columns": ["sample_id", "organism", "data_id", "design",
                                                     "group", "assay", "strain", "fastq_1", "fastq_2"]},
    "ncRNAseq":        {"kind": "fastq", "columns": FASTQ_COLS},
    "PeakCalling":     {"kind": "fastq", "columns": FASTQ_COLS},
    "tRNAseq":         {"kind": "fastq", "columns": FASTQ_COLS},
    "CLIP":            {"kind": "fastq", "columns": FASTQ_COLS},
    "PacVar":          {"kind": "bam_pbi", "columns": ["sample_id", "organism", "data_id", "design", "group", "bam", "pbi"]},
    "scRNAseq":        {"kind": "fastq", "columns": FASTQ_COLS},
    "LRtranscriptome": {"kind": "fastq", "columns": FASTQ_COLS + ["organism"]},
    "RNAseq":          {"kind": "fastq", "columns": ["sample_id", "organism", "group", "data_id",
                                                     "design", "fastq_1", "fastq_2", "contaminated_organism"]},
    "Fiberseq":        {"kind": "bam_pbi", "columns": ["sample_id", "organism", "data_id", "design", "group", "bam", "pbi"]},
    "QuantMS":         {"kind": "ms", "columns": ["sample_id", "data_id", "design", "ms_file"]},
}

WORKFLOW_DEFAULT_ORGANISM = {
    "PacVar": "T2T_mhaESC_v1.5",
    "Fiberseq": "T2T_mhaESC_v1.5",
}

FQ_PE_RE = re.compile(
    r"^(?P<base>.+?)[._-](?:R(?P<read_r>[12])|(?P<read_d>[12]))"
    r"(?:(?:[_-](?:00)?\d+)|(?:\.raw)*)*"
    r"\.(?:f(?:ast)?q)(?:\.gz)?$",
    re.IGNORECASE,
)
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
            groups[m.group("base")][m.group("read_r") or m.group("read_d")].append(f)
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
        stem = re.sub(r"\.(hifi_reads|subreads|ccs|polymerase|scraps)$", "", bam.name[:-4])
        units.append((stem, bam, pbi))
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
    ap.add_argument("-d", "--search-dir", required=True, nargs="+", action="extend",
                    help="one or more directories to scan recursively; repeatable "
                         "(-d A -d B == -d A B); tag with data_type as DIR=DNA / DIR=RNA")
    ap.add_argument("-o", "--out", default=None, help="output TSV (default meta_<workflow>.tsv)")
    ap.add_argument("--kind", choices=["fastq", "bam_pbi", "ms"], default=None,
                    help="override default input kind (e.g. --kind fastq for QuantMS test data)")
    ap.add_argument("--layout", choices=["any", "pe", "se"], default="any",
                    help="fastq layout filter (any = auto-detect)")
    ap.add_argument("--organism", default="", help="fill organism column (workflows that have it); "
                    "when empty, defaults per workflow (PacVar/Fiberseq -> T2T_mhaESC_v1.5)")
    ap.add_argument("--contaminated-organism", default="", help="fill contaminated_organism (RNAseq)")
    ap.add_argument("--group", default="", help="fill group column (RNAseq, Mutation)")
    ap.add_argument("--assay", default="", help="fill assay column, e.g. WGS / Exome / Panel "
                                               "(Mutation, DR DNA rows; drives CNVkit --seq-method "
                                               "and exome capture targets)")
    ap.add_argument("--assay-map", default=None, help="TSV: data_id or sample_id<TAB>assay "
                                                      "(DR: fill DNA rows only, RNA rows stay empty)")
    ap.add_argument("--strain", default="", help="fill strain column, e.g. 129P2_OlaHsd (Mutation, DR; "
                                                  "enables strain-aware germline filtering, empty = all-strain resource)")
    ap.add_argument("--data-type", default="", help="fill data_type column: DNA or RNA (DR; "
                                                    "per-sample override via --data-type-map)")
    ap.add_argument("--data-type-map", default=None, help="TSV: data_id or sample_id<TAB>DNA|RNA "
                                                          "(DR metas mix DNA and RNA rows)")
    ap.add_argument("--sample-map", default=None, help="TSV: data_id<TAB>sample_id (merge replicates)")
    ap.add_argument("--design-map", default=None, help="TSV: data_id or sample_id<TAB>design")
    ap.add_argument("--group-map", default=None, help="TSV: data_id or sample_id<TAB>group")
    ap.add_argument("--design-strip-rep", action="store_true",
                    help="fill design from sample_id with replicate suffix stripped")
    ap.add_argument("--group-strip-rep", action="store_true",
                    help="fill group from sample_id with replicate suffix stripped (e.g. E14-1 -> E14)")
    args = ap.parse_args()

    spec = WORKFLOWS[args.workflow]
    kind = args.kind or spec["kind"]
    file_cols = {"fastq": ["fastq_1", "fastq_2"], "bam_pbi": ["bam", "pbi"], "ms": ["ms_file"]}[kind]
    columns = spec["columns"]
    if not set(file_cols) <= set(columns):
        base = ["sample_id", "data_id", "design"]
        extras = [c for c in columns if c not in base + ["fastq_1", "fastq_2", "bam", "pbi", "ms_file"]]
        columns = base + file_cols + extras

    search_specs = []
    for d in args.search_dir:
        base, eq, tag = d.rpartition("=")
        if eq and tag.strip().upper() in ("DNA", "RNA"):
            search_specs.append((base, tag.strip().upper()))
        else:
            search_specs.append((d, ""))

    def scan_with(k):
        found = []
        seen_dirs = set()
        for base, dt in search_specs:
            sd = Path(base).resolve()
            if not sd.is_dir():
                die(f"search directory not found: {sd}")
            if sd in seen_dirs:
                logger.warning(f"search directory listed twice, ignoring repeat: {sd}")
                continue
            seen_dirs.add(sd)
            for unit in {"fastq": find_fastq_units, "bam_pbi": find_bam_pbi_units, "ms": find_ms_units}[k](sd):
                found.append(tuple(unit) + (dt,))
        return found

    units = scan_with(kind)
    if not units and not args.kind:
        for alt in ("fastq", "bam_pbi", "ms"):
            if alt == kind:
                continue
            alt_units = scan_with(alt)
            if alt_units:
                logger.warning(f"no {kind} files found, auto-detected {alt} input ({len(alt_units)} units)")
                kind = alt
                file_cols = {"fastq": ["fastq_1", "fastq_2"], "bam_pbi": ["bam", "pbi"], "ms": ["ms_file"]}[kind]
                if not set(file_cols) <= set(columns):
                    base = ["sample_id", "data_id", "design"]
                    extras = [c for c in columns if c not in base + ["fastq_1", "fastq_2", "bam", "pbi", "ms_file"]]
                    columns = base + file_cols + extras
                units = alt_units
                break
    if kind == "fastq" and args.layout != "any":
        want_pe = args.layout == "pe"
        units = [u for u in units if (u[2] is not None) == want_pe]
    if not units:
        die(f"no {kind} sequencing files found under {', '.join(args.search_dir)}")

    sample_map = load_map(args.sample_map)
    design_map = load_map(args.design_map)
    group_map = load_map(args.group_map)
    data_type_map = load_map(args.data_type_map)
    assay_map = load_map(args.assay_map)

    def norm_data_type(value, what):
        v = (value or "").strip().upper()
        if v and v not in ("DNA", "RNA"):
            die(f"data_type '{value}' ({what}): expected DNA or RNA, or empty")
        return v

    default_data_type = norm_data_type(args.data_type, "--data-type")

    typed = []
    for did, f1, f2, dir_dt in units:
        dt = norm_data_type(data_type_map.get(did) or dir_dt or default_data_type,
                            f"data_type for {did}")
        typed.append((did, f1, f2, dt))

    deduped = []
    seen_files = set()
    for did, f1, f2, dt in typed:
        key = (str(f1), str(f2) if f2 else "")
        if key in seen_files:
            logger.info(f"{did}: identical files already collected, skipping repeat scan hit")
            continue
        seen_files.add(key)
        deduped.append((did, f1, f2, dt))

    total = {}
    for did, _, _, dt in deduped:
        total[(dt, did)] = total.get((dt, did), 0) + 1
    reserved = set(total)
    emitted = set()
    next_no = {}
    units = []
    for did, f1, f2, dt in deduped:
        final = did
        scope = (dt, did)
        if total[scope] > 1:
            n = next_no.get(scope, 0) + 1
            final = f"{did}-{n}"
            while (dt, final) in reserved or (dt, final) in emitted:
                n += 1
                final = f"{did}-{n}"
            next_no[scope] = n
            logger.warning(f"duplicate data_id '{did}' under data_type '{dt or '-'}' "
                           f"({os.path.dirname(os.path.abspath(f1))}) -> renamed to '{final}'")
        emitted.add((dt, final))
        units.append((did, final, dt, f1, f2))

    out = Path(args.out) if args.out else Path(f"meta_{args.workflow}.tsv")
    n = 0
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for did, final_did, dt, f1, f2 in units:
            sid = sample_map.get(did, final_did)
            if re.search(r"\s", sid):
                die(f"sample_id '{sid}' contains whitespace (MetaUtil rejects it)")
            design = design_map.get(did) or design_map.get(sid) or ""
            if not design and args.design_strip_rep:
                design = REP_SUFFIX_RE.sub("", sid)
            group = group_map.get(did) or group_map.get(sid) or ""
            if not group and args.group_strip_rep:
                group = REP_SUFFIX_RE.sub("", sid)
            row = {"sample_id": sid, "data_id": final_did, "design": design, "group": group}
            if kind == "fastq":
                row["fastq_1"] = os.path.abspath(f1)
                row["fastq_2"] = os.path.abspath(f2) if f2 else ""
            elif kind == "bam_pbi":
                row["bam"], row["pbi"] = os.path.abspath(f1), os.path.abspath(f2)
            else:
                row["ms_file"] = os.path.abspath(f1)
            row["organism"] = args.organism or WORKFLOW_DEFAULT_ORGANISM.get(args.workflow, "")
            row["assay"] = "" if dt == "RNA" else (assay_map.get(did) or assay_map.get(sid) or args.assay)
            row["strain"] = args.strain
            row["data_type"] = dt
            if not row.get("group"):
                row["group"] = args.group
            row["contaminated_organism"] = args.contaminated_organism
            w.writerow([row.get(c, "") for c in columns])
            n += 1

    logger.info(f"workflow={args.workflow} kind={kind} -> {out} ({n} samples, columns: {','.join(columns)})")


if __name__ == "__main__":
    main()
