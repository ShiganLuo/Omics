#!/usr/bin/env python3
"""Cross-match Arriba fusion breakpoints against decoil ecDNA breakpoints.

Evidence level 1 in METHOD.md: an Arriba fusion whose breakpoint falls near an
ecDNA fragment boundary suggests the ecDNA is transcribed as a fusion
transcript. This is indirect evidence — the fusion may originate from a linear
rearrangement — and must be confirmed with junction-spanning reads on the
reconstructed ecDNA reference (see align_ecdna_junction.py).

Reads:
    decoil  reconstruct.ecDNA.filtered.bed   (#chr start end circ_id ...)
    decoil  reconstruct.links.txt (optional, unused for matching but validated)
    arriba  <sample>_passed_fusions.tsv and/or <sample>_discarded_fusions.tsv
            (Arriba TSV: breakpoint1 = "chrN:pos", breakpoint2 = "chrN:pos")

Writes:
    <out>/arriba_decoil_matches.tsv
        One row per (fusion, ecDNA, breakpoint-side) match with distance and
        the ecDNA fragment the breakpoint falls on/near.
    <out>/arriba_decoil_summary.tsv
        One row per ecDNA circ_id: matched fusion genes, min distance, evidence.

References (METHOD.md):
    Zhu et al., Nat Genet 2023, doi:10.1038/s41588-023-01386-y
    Turner et al., Nature 2017, doi:10.1038/nature21356
"""
from __future__ import annotations

import argparse
import logging
import os
import sys

LogUtil = None  # set up below; scripts in this repo use src.common.util.LogUtil

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
)
logger = logging.getLogger("match_arriba_decoil")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--decoil-dir", required=True,
                        help="decoil output dir containing reconstruct.ecDNA.filtered.bed")
    parser.add_argument("--arriba-tsv", required=True,
                        help="Arriba <sample>_passed_fusions.tsv")
    parser.add_argument("--arriba-discarded", default=None,
                        help="Arriba <sample>_discarded_fusions.tsv (low-confidence fusions still carry breakpoints)")
    parser.add_argument("--window", type=int, default=1000,
                        help="Max distance (bp) between Arriba breakpoint and ecDNA fragment boundary (default: 1000)")
    parser.add_argument("--out", required=True, help="Output directory")
    return parser.parse_args(argv)


def load_ecdna_fragments(bed_path: str) -> list[dict]:
    """Parse decoil reconstruct.ecDNA.filtered.bed into fragment records."""
    fragments = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                logger.warning("Skip malformed bed line: %s", line)
                continue
            fragments.append({
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "circ_id": fields[3],
                "fragment_id": fields[4],
                "strand": fields[5],
                "coverage": float(fields[6]) if len(fields) > 6 and fields[6] else None,
            })
    logger.info("Loaded %d ecDNA fragments from %s", len(fragments), bed_path)
    return fragments


def parse_breakpoint(bp: str) -> tuple[str, int] | None:
    """Parse an Arriba breakpoint 'chrN:pos' into (chrom, pos)."""
    if not bp or bp == "." or ":" not in bp:
        return None
    chrom, _, pos = bp.rpartition(":")
    try:
        return chrom, int(pos)
    except ValueError:
        logger.warning("Unparseable breakpoint: %r", bp)
        return None


def load_arriba_fusions(paths: list[str]) -> list[dict]:
    """Load Arriba fusion TSVs (passed + discarded) with a source column."""
    fusions = []
    for path, source in zip(paths, ["passed", "discarded"]):
        if not path or not os.path.isfile(path):
            if path:
                logger.warning("Arriba TSV not found, skip: %s", path)
            continue
        with open(path) as fh:
            header = None
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    header = line.lstrip("#").split("\t")
                    continue
                if header is None:
                    header = line.split("\t")
                    continue
                fields = line.split("\t")
                row = dict(zip(header, fields))
                bp1 = parse_breakpoint(row.get("breakpoint1", ""))
                bp2 = parse_breakpoint(row.get("breakpoint2", ""))
                if bp1 is None or bp2 is None:
                    continue
                fusions.append({
                    "source": source,
                    "gene1": row.get("gene1", ""),
                    "gene2": row.get("gene2", ""),
                    "bp1": bp1,
                    "bp2": bp2,
                    "type": row.get("type", ""),
                    "split_reads": row.get("split_reads1", ""),
                    "confidence": row.get("confidence", ""),
                })
        logger.info("Loaded %d fusions (cumulative) from %s [%s]", len(fusions), path, source)
    return fusions


def match_breakpoint_to_fragments(bp: tuple[str, int],
                                  fragments: list[dict],
                                  window: int) -> list[dict]:
    """Return ecDNA fragments whose boundary lies within *window* of *bp*.

    A breakpoint "falls on" an ecDNA fragment when it is inside the fragment or
    within *window* of either boundary: circularization can interrupt a gene
    anywhere, and RNA splicing may drop a few bases around the junction, so the
    RNA breakpoint need not coincide exactly with the DNA circle boundary.
    """
    chrom, pos = bp
    hits = []
    for frag in fragments:
        if frag["chrom"] != chrom:
            continue
        if frag["start"] - window <= pos <= frag["end"] + window:
            dist = 0 if frag["start"] <= pos <= frag["end"] else min(
                abs(pos - frag["start"]), abs(pos - frag["end"]))
            hits.append({**frag, "distance": dist})
    hits.sort(key=lambda h: h["distance"])
    return hits


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    os.makedirs(args.out, exist_ok=True)

    bed_path = os.path.join(args.decoil_dir, "reconstruct.ecDNA.filtered.bed")
    if not os.path.isfile(bed_path):
        logger.error("Missing decoil bed: %s", bed_path)
        return 1
    fragments = load_ecdna_fragments(bed_path)
    fusions = load_arriba_fusions([args.arriba_tsv, args.arriba_discarded])
    if not fusions:
        logger.error("No fusions loaded from %s", args.arriba_tsv)
        return 1

    match_rows: list[dict] = []
    per_circ: dict[str, dict] = {}

    for fusion in fusions:
        for side, bp in (("breakpoint1", fusion["bp1"]), ("breakpoint2", fusion["bp2"])):
            hits = match_breakpoint_to_fragments(bp, fragments, args.window)
            for hit in hits[:5]:  # cap per-breakpoint reporting at 5 nearest fragments
                row = {
                    "fusion_genes": f"{fusion['gene1']}--{fusion['gene2']}",
                    "fusion_source": fusion["source"],
                    "fusion_type": fusion["type"],
                    "confidence": fusion["confidence"],
                    "breakpoint_side": side,
                    "breakpoint": f"{bp[0]}:{bp[1]}",
                    "circ_id": hit["circ_id"],
                    "fragment_id": hit["fragment_id"],
                    "fragment": f"{hit['chrom']}:{hit['start']}-{hit['end']}",
                    "distance_bp": hit["distance"],
                }
                match_rows.append(row)
                agg = per_circ.setdefault(hit["circ_id"], {
                    "circ_id": hit["circ_id"],
                    "fusions": set(),
                    "min_distance": hit["distance"],
                    "sources": set(),
                })
                agg["fusions"].add(row["fusion_genes"])
                agg["min_distance"] = min(agg["min_distance"], hit["distance"])
                agg["sources"].add(fusion["source"])

    match_tsv = os.path.join(args.out, "arriba_decoil_matches.tsv")
    with open(match_tsv, "w") as fh:
        cols = ["fusion_genes", "fusion_source", "fusion_type", "confidence",
                "breakpoint_side", "breakpoint", "circ_id", "fragment_id",
                "fragment", "distance_bp"]
        fh.write("\t".join(cols) + "\n")
        for row in match_rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    summary_tsv = os.path.join(args.out, "arriba_decoil_summary.tsv")
    with open(summary_tsv, "w") as fh:
        fh.write("circ_id\tn_fusions\tmin_distance_bp\tfusion_genes\tsources\tevidence\n")
        for circ_id in sorted(per_circ, key=lambda c: per_circ[c]["min_distance"]):
            agg = per_circ[circ_id]
            if "passed" in agg["sources"]:
                evidence = "fusion_breakpoint_on_ecDNA"
            else:
                # discarded-only hits come from Arriba's min_support rejects
                # (hundreds of thousands of noise chimeras); treat as
                # noise-prone, do not claim expression from these alone.
                evidence = "discarded_only_noise_prone"
            fh.write(f"{circ_id}\t{len(agg['fusions'])}\t{agg['min_distance']}\t"
                     f"{';'.join(sorted(agg['fusions']))}\t{';'.join(sorted(agg['sources']))}\t"
                     f"{evidence}\n")

    logger.info("Matches: %d rows across %d ecDNAs -> %s",
                len(match_rows), len(per_circ), summary_tsv)
    if not match_rows:
        logger.info("No Arriba breakpoint falls within %d bp of any ecDNA fragment.", args.window)
    return 0


if __name__ == "__main__":
    sys.exit(main())
