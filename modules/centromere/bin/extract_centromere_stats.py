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


def write_stats(rm_out: str, output: str, species: str) -> None:
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

    print(f"Wrote: {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rm_out", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--species", default="mouse",
                        choices=list(CENTROMERE_FAMILIES.keys()))
    args = parser.parse_args()
    write_stats(args.rm_out, args.output, args.species)


if __name__ == "__main__":
    main()
