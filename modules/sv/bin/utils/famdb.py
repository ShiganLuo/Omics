#!/usr/bin/env python3
"""famdb.py — extract a species-specific TE consensus library from a Dfam
FamDB HDF5 file (e.g. dfam39_full.7.h5), without the full famdb toolchain.

Produces a RepeatMasker ``-lib``-ready FASTA of consensus sequences for
families that are active in the requested species, mirroring
``famdb.py -i dfam.h5 curate --taxon <species>`` (stage + threshold
filtered).

Selection model (from the FamDB format):
  - Every family carries ``clades`` (NCBI taxon ids where its HMM is
    searched) and ``taxa_thresholds`` (per-taxon curated score rows
    ``taxon, threshold, ...``). A family belongs to a species when the
    species taxid (or one of its ancestor taxa) appears in either.
  - A ``stages`` mask optionally restricts searches to specific
    RepeatMasker stages; families whose mask excludes all default
    stages are skipped by curate, so we honour it.

Header format follows FamDB's curated-library convention so
RepeatMasker prints species-specific cutoffs:
  >NAME#TYPE/SUBTYPE @Clade [S:stage LL,UL]

Performance notes (dfam39_full.7.h5 is ~60 GB on spinning rust):
  - family accessions are enumerated through GROUP links only (datasets
    are never opened during enumeration);
  - selection attrs (name/class/clades/thresholds/stages) are read before
    the consensus, and the consensus attr is only fetched for families
    that pass selection;
  - ``--threads N`` scans the family tree in parallel worker processes
    (threads are useless here: the HDF5 C library serialises access
    behind a global lock).

Compressed input: ``-i dfam.h5.gz`` is accepted. When the uncompressed
sibling exists it is used directly; otherwise the archive is decompressed
with pigz (fallback: gzip) next to the archive (or into ``--tmp-dir``),
after a free-space check. ``-o out.fa.gz`` writes gzip-compressed output.

Usage:
  python famdb.py -i dfam39_full.7.h5 -s mouse -o mouse_te_lib.fa
  python famdb.py -i dfam39_full.7.h5.gz -s mouse -o mouse_te_lib.fa.gz
  python famdb.py -i dfam39_full.7.h5 --taxon 10090 -o mouse_te_lib.fa
  python famdb.py -i dfam39_full.7.h5 -s human -o human_te_lib.fa \\
      --threads 8
"""
import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import h5py

logger = logging.getLogger("famdb")

# Common species aliases -> NCBI taxid (extend as needed)
SPECIES_TAXID = {
    "human": 9606,
    "homo sapiens": 9606,
    "mouse": 10090,
    "mus musculus": 10090,
    "rat": 10116,
    "rattus norvegicus": 10116,
    "zebrafish": 7955,
    "danio rerio": 7955,
    "drosophila": 7227,
    "worm": 6239,
    "c. elegans": 6239,
    "yeast": 4932,
    "arabidopsis": 3702,
    "rice": 39947,
    "maize": 4577,
    "rhesus": 9544,
    "macaca mulatta": 9544,
    "macaque": 9544,
    "chimp": 9598,
    "dog": 9615,
    "cow": 9913,
    "chicken": 9031,
    "frog": 8355,
    "mulatta": 9544,
}

DEFAULT_STAGES = (1, 40, 60, 65, 80)


def setup_logger(log_file: Optional[str] = None) -> None:
    handlers: List[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
    )


def resolve_taxon(species: Optional[str], taxon: Optional[int]) -> int:
    if taxon:
        return taxon
    if not species:
        raise SystemExit("either --species or --taxon is required")
    key = species.strip().lower()
    if key in SPECIES_TAXID:
        return SPECIES_TAXID[key]
    raise SystemExit(
        f"unknown species '{species}' — pass its NCBI taxid via --taxon "
        f"(known aliases: {', '.join(sorted(set(SPECIES_TAXID)))})"
    )

def parse_taxa_thresholds(raw: str) -> Dict[int, Tuple[float, float, float, float]]:
    """Parse the per-taxon curated rows: taxon, th, th_anchored, th_raw, div."""
    out: Dict[int, Tuple[float, float, float, float]] = {}
    for line in (raw or "").splitlines():
        f = line.strip().split(",")
        if len(f) >= 5:
            try:
                out[int(f[0])] = tuple(float(x) for x in f[1:5])  # type: ignore[assignment]
            except ValueError:
                continue
    return out

def parse_stages(attr) -> Optional[Set[int]]:
    """buffer_stages attr -> set of stages this family participates in."""
    if attr is None:
        return None
    stages: Set[int] = set()
    for m in re.finditer(r"(\d+)\[", str(attr)):
        stages.add(int(m.group(1)))
    # single-number form means "that stage only"
    if not stages:
        m = re.fullmatch(r"\s*(\d+)\s*", str(attr))
        if m:
            stages.add(int(m.group(1)))
    return stages or None


# ----------------------------------------------------------- gz input/output
def _free_bytes(path: Path) -> int:
    return shutil.disk_usage(str(path)).free


def resolve_h5_input(path: str, tmp_dir: Optional[str]) -> Tuple[str, Optional[Path]]:
    """Return (h5_path, cleanup_path).

    ``.gz``/``.bgz`` inputs are accepted: when the uncompressed sibling
    exists it is used directly, otherwise the archive is decompressed with
    pigz (fallback: the gzip module) next to the archive, or into
    ``tmp_dir`` when given. The decompressed file is KEPT (it is the .h5
    the caller wants anyway, and future runs hit the sibling fast path);
    free space is checked first.
    """
    p = Path(path)
    if not p.name.endswith((".gz", ".bgz")):
        return str(p), None

    stem = p.name[: -len(".gz")] if p.name.endswith(".gz") else p.name[: -len(".bgz")]
    sibling = p.with_name(stem)
    if sibling.exists() and sibling.stat().st_size > 0:
        logger.info("gz input: using uncompressed sibling %s", sibling)
        return str(sibling), None

    target_dir = Path(tmp_dir) if tmp_dir else p.parent
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / stem

    # Dfam FamDB h5 compresses ~7x; require 10x + 1 GB margin
    need = p.stat().st_size * 10 + (1 << 30)
    have = _free_bytes(target_dir)
    if have < need:
        raise SystemExit(
            f"not enough space to decompress {p.name}: need ~{need / 1e9:.0f} GB "
            f"free in {target_dir}, have {have / 1e9:.0f} GB — free space or "
            f"point --tmp-dir elsewhere"
        )

    logger.info("decompressing %s -> %s (this takes a few minutes)", p, target)
    pigz = shutil.which("pigz")
    if pigz:
        with open(target, "wb") as fh:
            r = subprocess.run([pigz, "-dc", str(p)], stdout=fh,
                               stderr=subprocess.PIPE)
        if r.returncode != 0:
            target.unlink(missing_ok=True)
            raise SystemExit(f"pigz failed: {r.stderr.decode()[:500]}")
    else:
        try:
            with gzip.open(p, "rb") as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst, length=64 << 20)
        except Exception as e:
            target.unlink(missing_ok=True)
            raise SystemExit(f"gzip decompression failed: {e}")
    logger.info("decompressed %s (%.1f GB)", target, target.stat().st_size / 1e9)
    return str(target), None


def open_output(path: str):
    """Open output FASTA, transparently gzipped when the name ends in .gz."""
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "wt", encoding="utf-8", compresslevel=6)
    return open(path, "w", encoding="utf-8")


# --------------------------------------------------------------- family scan
def enumerate_family_paths(h5: h5py.File) -> List[str]:
    """Accession -> object path for every family in Families/DF/**.

    Walks GROUP links only: family datasets are never opened here, so this
    is cheap even though the file is tens of GB (dataset object headers
    are the expensive random reads; groups are few).
    """
    out: List[str] = []

    def _walk(group, prefix: str) -> None:
        for key in group:
            if key.startswith("DF"):
                out.append(f"{prefix}/{key}")
            else:  # numeric bucket group
                _walk(group[key], f"{prefix}/{key}")

    _walk(h5["Families/DF"], "Families/DF")
    return out


def _to_str(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    return str(v)


def scan_chunk(h5path: str, chunk: List[Tuple[str, str]], taxid: int,
               min_len: int, use_stages: bool) -> List[Tuple[str, str, str]]:
    """Process [(accession, obj_path)] and return [(accession, header, seq)].

    One HDF5 handle per call; selection attrs are read first and the
    consensus attr only for families that pass selection.
    """
    kept: List[Tuple[str, str, str]] = []
    with h5py.File(h5path, "r") as h5:
        for acc, path in chunk:
            node = h5[path]
            attrs = node.attrs

            clades = attrs.get("clades")
            clade_ids = set(int(x) for x in clades) if clades is not None else set()
            thresholds = parse_taxa_thresholds(_to_str(attrs.get("taxa_thresholds", "")))
            if taxid not in clade_ids and taxid not in thresholds:
                continue

            stages = parse_stages(attrs.get("buffer_stages"))
            if use_stages and stages is not None and not (stages & set(DEFAULT_STAGES)):
                continue

            consensus = attrs.get("consensus")
            if consensus is None or len(consensus) == 0:
                continue
            consensus = _to_str(consensus)
            if len(consensus) < min_len:
                continue

            name = _to_str(attrs.get("name")) or acc
            rtype = _to_str(attrs.get("repeat_type")) or "Unknown"
            rsub = _to_str(attrs.get("repeat_subtype")) or ""
            type_field = f"{rtype}/{rsub}" if rsub else rtype

            stage_info = ""
            if taxid in thresholds:
                th = thresholds[taxid][0]
                th_up = th + 5.0
                active = sorted((stages or set(DEFAULT_STAGES)) & set(DEFAULT_STAGES))
                stage_info = f" [S:{','.join(str(s) for s in active)} {th:.0f},{th_up:.0f}]"
            header = f"{name}#{type_field} @{clade_name(taxid)}{stage_info}"
            kept.append((acc, header, consensus))
    return kept


def curate_library(h5path: str, taxid: int, min_len: int, use_stages: bool,
                   threads: int = 1) -> List[Tuple[str, str]]:
    """Return [(header, consensus)] for families active in the taxon."""
    with h5py.File(h5path, "r") as h5:
        paths = enumerate_family_paths(h5)
    logger.info("enumerated %d family accessions", len(paths))
    n_total = len(paths)

    accs = [(p.rsplit("/", 1)[-1], p) for p in paths]
    accs.sort(key=lambda x: x[0])  # deterministic: first-wins dedupe by accession

    threads = max(1, min(threads, 16))
    chunks: List[List[Tuple[str, str]]] = []
    size = max(64, len(accs) // (threads * 4) or 1)
    for i in range(0, len(accs), size):
        chunks.append(accs[i : i + size])

    results: List[Tuple[str, str, str]] = []
    if threads == 1 or len(chunks) == 1:
        for ch in chunks:
            results.extend(scan_chunk(h5path, ch, taxid, min_len, use_stages))
    else:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futs = [pool.submit(scan_chunk, h5path, ch, taxid, min_len, use_stages)
                    for ch in chunks]
            for fut in futs:
                results.extend(fut.result())

    results.sort(key=lambda x: x[0])
    out: List[Tuple[str, str]] = []
    seen: Set[str] = set()
    for _acc, header, seq in results:
        if header in seen:
            continue
        seen.add(header)
        out.append((header, seq))

    logger.info("curated %d/%d families for taxon %d (min_len=%d)",
                len(out), n_total, taxid, min_len)
    return out

def clade_name(taxid: int) -> str:
    """Clade label: taxonomy names are not stored per-family in FamDB; use
    the species alias when known, else the taxid."""
    for alias, t in SPECIES_TAXID.items():
        if t == taxid:
            return alias.replace(" ", "_").capitalize()
    return str(taxid)


# Sentinel subfamilies per species alias: a curated species library MUST
# contain at least one family from each of these (active-clade marker
# names). A library lacking them was built against the wrong lineage
# branch (observed 2026-10-29: mouse_te_lib.fa carried only
# Amniota/Vertebrata-level families and silently classified every mouse
# insertion as noise). Checked at build time (refuse to write) and usable
# by consumers via validate_library().
SENTINEL_SUBFAMILIES = {
    "mouse": [["L1Md", "L1M", "L1Mus"], ["B1_Mm", "B1_Mus", "B1"],
              ["B2_Mm", "B2"], ["IAP"], ["ETn", "MERVL"]],
    "human": [["L1PA", "L1HS", "L1"], ["AluY", "AluJ", "AluS", "Alu"],
              ["SVA"], ["HERV", "LTR"]],
    "rhesus": [["L1PA", "L1"], ["AluY", "AluJ", "AluS", "Alu"],
               ["SVA"], ["HERV", "LTR"]],
}

def validate_library(records, species_alias: str) -> List[str]:
    """Return list of missing sentinel groups ([] = library OK).

    records: [(header, seq)] as produced by curate_library. A sentinel
    group passes when ANY of its names appears as a word-start in some
    header family name (case-insensitive). Deterministic: catches
    wrong-branch/wrong-species curation and empty libraries; it is not a
    checksum for sequence-level corruption.
    """
    groups = SENTINEL_SUBFAMILIES.get(species_alias)
    if not groups:  # unmapped species: no sentinels known, nothing to check
        return []
    headers = [h.lower() for h, _ in records]
    missing = []
    for group in groups:
        if not any(any(g.lower() in h for h in headers) for g in group):
            missing.append("/".join(group))
    return missing

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Extract a species-specific TE consensus library from a Dfam FamDB HDF5 file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Dfam FamDB .h5 file (or .h5.gz — decompressed on demand)")
    ap.add_argument("-s", "--species", default=None, help="species name/alias (e.g. mouse, human)")
    ap.add_argument("--taxon", type=int, default=None, help="NCBI taxid (alternative to --species)")
    ap.add_argument("-o", "--output", required=True,
                    help="output consensus FASTA (RepeatMasker -lib ready; .gz name -> gzip output)")
    ap.add_argument("--min-length", type=int, default=20, help="minimum consensus length")
    ap.add_argument("--line-width", type=int, default=60, help="FASTA sequence line width")
    ap.add_argument("--no-stage-filter", action="store_true",
                    help="keep families regardless of their buffer_stages mask")
    ap.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1),
                    help="parallel scan workers (processes; one HDF5 handle each)")
    ap.add_argument("--tmp-dir", default=None,
                    help="directory for on-demand gz decompression (default: next to the archive)")
    ap.add_argument("-l", "--log", default=None, help="log file path")
    args = ap.parse_args(argv)

    setup_logger(args.log)
    taxid = resolve_taxon(args.species, args.taxon)
    h5path, cleanup = resolve_h5_input(args.input, args.tmp_dir)
    logger.info("Reading %s for taxon %d (threads=%d)", h5path, taxid, args.threads)

    try:
        records = curate_library(h5path, taxid, args.min_length,
                                 use_stages=not args.no_stage_filter,
                                 threads=args.threads)
    finally:
        if cleanup is not None:
            cleanup.unlink(missing_ok=True)

    if not records:
        logger.error("no families matched taxon %d — wrong taxid or unsupported file?", taxid)
        return 1

    # refuse to write a wrong-branch library: every sentinel group of the
    # species must be present (build-time half of the integrity check)
    alias = next((a for a, t in SPECIES_TAXID.items()
                  if t == taxid and " " not in a), str(taxid))
    missing = validate_library(records, alias)
    if missing:
        logger.error(
            "curated library for '%s' is missing sentinel families: %s — "
            "the curation likely selected the wrong lineage branch; "
            "NOT writing %s", alias, ", ".join(missing), args.output)
        return 1

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open_output(str(out)) as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), args.line_width):
                fh.write(seq[i : i + args.line_width] + "\n")
    logger.info("Wrote %d consensus sequences to %s", len(records), out)
    return 0

if __name__ == "__main__":
    sys.exit(main())
