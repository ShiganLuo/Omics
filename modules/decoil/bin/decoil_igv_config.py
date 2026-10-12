#!/usr/bin/env python3
"""Build an IGV.js browser config for Decoil ecDNA results.

Generates the config JSON consumed by ``modules/track/bin/track.py`` in
``--mode igv``. The browser shows, per sample:

- coverage bigwig (read depth track)
- ecDNA fragments BED (``reconstruct.ecDNA.filtered.bed``) as an
  annotation track
- reference GTF gene annotation (from the workflow igv config)

The default ``locus`` is the largest filtered ecDNA region of the first
sample, so the browser opens directly on the most relevant view.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List, Optional


def parse_pairs(pairs: Optional[List[str]], flag: str) -> List:
    """Parse repeated ``SAMPLE:PATH`` CLI arguments.

    Parameters
    ----------
    pairs : list of str or None
        Raw argument values.
    flag : str
        Flag name used in error messages.

    Returns
    -------
    list of tuple(str, str)
        (sample_id, path) pairs; path may contain colons.

    Raises
    ------
    SystemExit
        On malformed entries or missing files.
    """
    parsed = []
    for item in pairs or []:
        if ":" not in item:
            sys.exit(f"Invalid {flag} entry '{item}', expected SAMPLE:PATH")
        name, path = item.split(":", 1)
        if not os.path.isfile(path):
            sys.exit(f"File not found for {flag}: {path}")
        parsed.append((name, path))
    return parsed


def largest_ecdna_locus(bed_path: str, default_locus: str) -> str:
    """Return ``chr:start-end`` of the largest region in a BED.

    Parameters
    ----------
    bed_path : str
        BED path (comment lines tolerated).
    default_locus : str
        Fallback when the BED is empty/unreadable.

    Returns
    -------
    str
        Locus string with full coordinates (never abbreviated).
    """
    best = None
    best_len = -1
    try:
        with open(bed_path, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    continue
                if end - start > best_len:
                    best_len = end - start
                    best = f"{fields[0]}:{start}-{end}"
    except OSError as e:
        print(f"WARNING: cannot read {bed_path}: {e}", file=sys.stderr)
    return best or default_locus


def build_config(args: argparse.Namespace) -> Dict:
    """Assemble the IGV browser config dict.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments (sample/bw/bed pairs, igv base config, locus).

    Returns
    -------
    dict
        IGV.js browser config: reference + base tracks (GTF) + per-sample
        bigwig/bed tracks, with ``locus`` set for the initial view.
    """
    with open(args.igv_config, "r", encoding="utf-8") as handle:
        base = json.load(handle)

    sample_pairs = parse_pairs(args.sample, "--sample")
    bw_pairs = parse_pairs(args.bw, "--bw")
    if len(sample_pairs) != len(bw_pairs):
        sys.exit("--sample and --bw must have the same count")

    tracks: List[Dict] = []
    # Reference annotation tracks first (loaded from the workflow igv config)
    for track in base.get("tracks", []):
        tracks.append(track)

    locus = args.locus or ""
    for i, (sid, bed_path) in enumerate(sample_pairs):
        bw_path = bw_pairs[i][1]
        tracks.append({
            "name": f"{sid} coverage",
            "url": bw_path,
            "format": "bigwig",
            "type": "wig",
            "color": "#1f77b4",
            "displayMode": "COLLAPSED",
            "autoscale": True,
            "height": 40,
        })
        tracks.append({
            "name": f"{sid} ecDNA",
            "url": bed_path,
            "format": "bed",
            "type": "annotation",
            "color": "#E74C3C",
            "height": 40,
        })
        if not locus:
            locus = largest_ecdna_locus(bed_path, "")

    config = {
        "id": base.get("id", "decoil"),
        "name": base.get("name", "Decoil ecDNA"),
        "js": base.get("js", "https://cdn.jsdelivr.net/npm/igv/dist/igv.min.js"),
        "publicPathMap": base.get("publicPathMap", []),
        "fastaURL": base.get("fastaURL"),
        "indexURL": base.get("indexURL"),
        "cytobandURL": base.get("cytobandURL"),
        "tracks": tracks,
    }
    if locus:
        config["locus"] = locus
    return config


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build IGV.js config JSON for Decoil ecDNA results.")
    parser.add_argument("--sample", action="append", required=True,
                        metavar="SAMPLE:ECDNA_BED",
                        help="sample id + reconstruct.ecDNA.filtered.bed; "
                             "repeatable")
    parser.add_argument("--bw", action="append", required=True,
                        metavar="SAMPLE:BIGWIG",
                        help="sample id + coverage bigwig; repeatable, same "
                             "order as --sample")
    parser.add_argument("--igv-config", required=True,
                        help="workflow igv config JSON (reference fasta/"
                             "index/cytoband/GTF + publicPathMap)")
    parser.add_argument("--locus", default="",
                        help="initial genomic locus; default: largest "
                             "ecDNA region of the first sample")
    parser.add_argument("-o", "--output", required=True,
                        help="output config JSON path")
    args = parser.parse_args()

    config = build_config(args)
    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2, ensure_ascii=False)
    print(f"IGV config written: {args.output} "
          f"(locus={config.get('locus', 'default')}, "
          f"tracks={len(config['tracks'])})")


if __name__ == "__main__":
    main()
