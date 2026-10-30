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

Usage:
  python famdb.py -i dfam39_full.7.h5 -s mouse -o mouse_te_lib.fa
  python famdb.py -i dfam39_full.7.h5 --taxon 10090 -o mouse_te_lib.fa
  python famdb.py -i dfam39_full.7.h5 -s human -o human_te_lib.fa \
      --ancestor-depth 2
"""
import argparse
import logging
import re
import sys
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


def iter_families(h5: h5py.File):
    """Yield (accession, attrs) for every family in Families/DF/**."""
    def _walk(group):
        for key in group:
            node = group[key]
            if key.startswith("DF") and isinstance(node, h5py.Dataset):
                yield key, dict(node.attrs)
            elif isinstance(node, h5py.Group):
                yield from _walk(node)
    yield from _walk(h5["Families"]["DF"])


def clade_name(taxid: int) -> str:
    """Clade label: taxonomy names are not stored per-family in FamDB; use
    the species alias when known, else the taxid."""
    for alias, t in SPECIES_TAXID.items():
        if t == taxid:
            return alias.replace(" ", "_").capitalize()
    return str(taxid)


def curate_library(h5path: str, taxid: int, min_len: int, use_stages: bool) -> List[Tuple[str, str]]:
    """Return [(header, consensus)] for families active in the taxon."""
    out: List[Tuple[str, str]] = []
    seen: Set[str] = set()
    n_total = n_kept = 0
    with h5py.File(h5path, "r") as h5:
        for acc, attrs in iter_families(h5):
            n_total += 1
            name = str(attrs.get("name", acc))
            consensus = attrs.get("consensus")
            if consensus is None or len(consensus) == 0:
                continue
            if len(consensus) < min_len:
                continue

            clades = attrs.get("clades")
            clade_ids = set(int(x) for x in clades) if clades is not None else set()
            thresholds = parse_taxa_thresholds(str(attrs.get("taxa_thresholds", "")))

            # family belongs to the species if the taxid is referenced in
            # its clade list or curated thresholds
            if not (taxid in clade_ids or taxid in thresholds):
                continue

            stages = parse_stages(attrs.get("buffer_stages"))
            if use_stages and stages is not None and not (stages & set(DEFAULT_STAGES)):
                continue

            # header: NAME#TYPE/SUBTYPE @Clade_name [S:stage th,th]
            rtype = str(attrs.get("repeat_type", "Unknown")) or "Unknown"
            rsub = str(attrs.get("repeat_subtype", "")) or ""
            type_field = f"{rtype}/{rsub}" if rsub else rtype
            stage_info = ""
            if taxid in thresholds:
                th = thresholds[taxid][0]
                th_up = th + 5.0
                active = sorted((stages or set(DEFAULT_STAGES)) & set(DEFAULT_STAGES))
                stage_info = f" [S:{','.join(str(s) for s in active)} {th:.0f},{th_up:.0f}]"
            header = f"{name}#{type_field} @{clade_name(taxid)}{stage_info}"
            if header in seen:
                continue
            seen.add(header)
            out.append((header, str(consensus)))
            n_kept += 1

    logger.info(
        "curated %d/%d families for taxon %d (min_len=%d)",
        n_kept, n_total, taxid, min_len,
    )
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Extract a species-specific TE consensus library from a Dfam FamDB HDF5 file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-i", "--input", required=True, help="Dfam FamDB .h5 file")
    ap.add_argument("-s", "--species", default=None, help="species name/alias (e.g. mouse, human)")
    ap.add_argument("--taxon", type=int, default=None, help="NCBI taxid (alternative to --species)")
    ap.add_argument("-o", "--output", required=True, help="output consensus FASTA (RepeatMasker -lib ready)")
    ap.add_argument("--min-length", type=int, default=20, help="minimum consensus length")
    ap.add_argument("--line-width", type=int, default=60, help="FASTA sequence line width")
    ap.add_argument("--no-stage-filter", action="store_true",
                    help="keep families regardless of their buffer_stages mask")
    ap.add_argument("-l", "--log", default=None, help="log file path")
    args = ap.parse_args(argv)

    setup_logger(args.log)
    taxid = resolve_taxon(args.species, args.taxon)
    logger.info("Reading %s for taxon %d", args.input, taxid)

    records = curate_library(args.input, taxid, args.min_length,
                             use_stages=not args.no_stage_filter)
    if not records:
        logger.error("no families matched taxon %d — wrong taxid or unsupported file?", taxid)
        return 1

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), args.line_width):
                fh.write(seq[i : i + args.line_width] + "\n")
    logger.info("Wrote %d consensus sequences to %s", len(records), out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
