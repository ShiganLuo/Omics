#!/usr/bin/env python3
"""Extract centromere satellite DNA statistics from RepeatMasker .out file.

Usage:
    python extract_centromere_stats.py --rm_out input.fa.out --output stats.txt
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

CENTROMERE_FAMILIES = {
    "mouse": {
        "MajorSatellite": ["MSAT", "Major_satellite", "GSAT"],
        "MinorSatellite": ["SATMIN", "Minor_satellite"],
        "Pericentromeric": ["SAT", "SATB"],
    },
    "human": {
        "AlphaSatellite": ["ALR", "Alpha"],
        "Satellite": ["SATA", "SATR", "HSAT"],
    },
}


def parse_repeatmasker_out(rm_out: str) -> list[dict]:
    records = []
    with open(rm_out) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score") or \
               line.startswith("Kimura") or line.startswith("("):
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            records.append({
                "div": parts[1],
                "query": parts[4],
                "q_begin": int(parts[5]),
                "q_end": int(parts[6]),
                "repeat": parts[9],
                "class_family": parts[10] if len(parts) > 10 else "",
            })
    return records


def extract_satellite_stats(records: list[dict], species: str) -> dict[str, list[int]]:
    families = CENTROMERE_FAMILIES.get(species, {})
    sat_stats: dict[str, list[int]] = defaultdict(list)
    for rec in records:
        cf = rec["class_family"]
        rn = rec["repeat"]
        span = abs(rec["q_end"] - rec["q_begin"]) + 1
        for label, patterns in families.items():
            for pat in patterns:
                if pat in cf.upper() or pat in rn.upper():
                    sat_stats[label].append(span)
                    break
        if "Satellite" in cf or "satellite" in cf:
            sat_stats["AllSatellite"].append(span)
    return dict(sat_stats)


def find_satellite_blocks(records: list[dict], species: str,
                          max_gap: int = 50_000,
                          min_block_len: int = 100_000) -> list[dict]:
    """Find contiguous satellite DNA blocks on each contig.

    Groups nearby satellite repeat hits into blocks, allowing gaps up to
    max_gap between consecutive hits. Returns blocks >= min_block_len.

    Args:
        records: parsed RepeatMasker records
        species: species name for centromere family lookup
        max_gap: max gap (bp) between hits to merge into one block
        min_block_len: minimum block length to report

    Returns:
        list of dicts with keys: contig, start, end, length, sat_types, n_regions
    """
    families = CENTROMERE_FAMILIES.get(species, {})
    sat_keywords = set()
    for patterns in families.values():
        sat_keywords.update(p.upper() for p in patterns)
    # Also include generic satellite
    sat_keywords.add("SATELLITE")

    # Group satellite hits by contig
    contig_hits: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for rec in records:
        cf = rec["class_family"].upper()
        rn = rec["repeat"].upper()
        is_sat = any(kw in cf or kw in rn for kw in sat_keywords)
        if not is_sat:
            continue
        start = min(rec["q_begin"], rec["q_end"])
        end = max(rec["q_begin"], rec["q_end"])
        contig_hits[rec["query"]].append((start, end, rec["class_family"] + "/" + rec["repeat"]))

    blocks = []
    for contig, hits in contig_hits.items():
        hits.sort(key=lambda x: x[0])
        # Cluster hits into blocks
        current_block = [hits[0]]
        for hit in hits[1:]:
            if hit[0] - current_block[-1][1] <= max_gap:
                current_block.append(hit)
            else:
                block_start = current_block[0][0]
                block_end = current_block[-1][1]
                block_len = block_end - block_start + 1
                if block_len >= min_block_len:
                    sat_types = set()
                    for h in current_block:
                        sat_types.add(h[2])
                    blocks.append({
                        "contig": contig,
                        "start": block_start,
                        "end": block_end,
                        "length": block_len,
                        "sat_types": sorted(sat_types),
                        "n_regions": len(current_block),
                    })
                current_block = [hit]
        # Last block
        if current_block:
            block_start = current_block[0][0]
            block_end = current_block[-1][1]
            block_len = block_end - block_start + 1
            if block_len >= min_block_len:
                sat_types = set()
                for h in current_block:
                    sat_types.add(h[2])
                blocks.append({
                    "contig": contig,
                    "start": block_start,
                    "end": block_end,
                    "length": block_len,
                    "sat_types": sorted(sat_types),
                    "n_regions": len(current_block),
                })

    blocks.sort(key=lambda x: -x["length"])
    return blocks


def write_stats(rm_out: str, output: str, species: str,
                max_gap: int = 50_000,
                min_block_len: int = 100_000) -> None:
    records = parse_repeatmasker_out(rm_out)
    sat_stats = extract_satellite_stats(records, species)
    total_span = sum(abs(r["q_end"] - r["q_begin"]) + 1 for r in records)

    with open(output, "w") as f:
        f.write(f"RepeatMasker output: {rm_out}\n")
        f.write(f"Species: {species}\n")
        f.write(f"Total repeats: {len(records)}\n")
        f.write(f"Total genome span: {total_span:,} bp\n\n")

        f.write("=== Centromere-Associated Satellite DNA ===\n\n")
        for label, spans in sorted(sat_stats.items()):
            if not spans:
                continue
            total = sum(spans)
            n = len(spans)
            f.write(f"  {label}:\n")
            f.write(f"    Regions:   {n}\n")
            f.write(f"    Total bp:  {total:,}\n")
            f.write(f"    Mean span: {total // n:,} bp\n")
            f.write(f"    Max span:  {max(spans):,} bp\n\n")

        class_totals: dict[str, int] = defaultdict(int)
        class_counts: dict[str, int] = defaultdict(int)
        for rec in records:
            cf = rec["class_family"]
            span = abs(rec["q_end"] - rec["q_begin"]) + 1
            class_totals[cf] += span
            class_counts[cf] += 1

        f.write("=== Top Repeat Classes by Total Bp ===\n\n")
        for cf, total in sorted(class_totals.items(), key=lambda x: -x[1])[:15]:
            f.write(f"  {cf:40s}  {class_counts[cf]:>8,} regions  {total:>12,} bp\n")

        # Per-contig satellite blocks
        blocks = find_satellite_blocks(records, species, max_gap, min_block_len)
        f.write(f"\n=== Per-Contig Satellite Blocks (>= {min_block_len:,} bp, gap <= {max_gap:,} bp) ===\n\n")
        if blocks:
            f.write(f"{'Contig':<20} {'Start':>12} {'End':>12} {'Length':>12} {'Regions':>8}  Sat Types\n")
            f.write(f"{'-'*20} {'-'*12} {'-'*12} {'-'*12} {'-'*8}  {'-'*40}\n")
            for b in blocks:
                types_str = ", ".join(b["sat_types"][:5])
                if len(b["sat_types"]) > 5:
                    types_str += f" (+{len(b['sat_types'])-5} more)"
                f.write(f"{b['contig']:<20} {b['start']:>12,} {b['end']:>12,} {b['length']:>12,} {b['n_regions']:>8}  {types_str}\n")

            # Summary
            total_block_bp = sum(b["length"] for b in blocks)
            f.write(f"\n  Total satellite block bp: {total_block_bp:,}\n")
            f.write(f"  Number of blocks: {len(blocks)}\n")
        else:
            f.write(f"  No satellite blocks >= {min_block_len:,} bp found.\n")

    print(f"Wrote: {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rm_out", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--species", default="mouse",
                        choices=list(CENTROMERE_FAMILIES.keys()))
    parser.add_argument("--max_gap", type=int, default=50_000,
                        help="Max gap (bp) between satellite hits to merge into one block (default: 50000)")
    parser.add_argument("--min_block_len", type=int, default=100_000,
                        help="Minimum satellite block length to report (default: 100000)")
    args = parser.parse_args()
    write_stats(args.rm_out, args.output, args.species,
                max_gap=args.max_gap, min_block_len=args.min_block_len)


if __name__ == "__main__":
    main()
