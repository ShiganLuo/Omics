#!/usr/bin/env python3
"""Breakpoint-sharing families of ecDNA circles across samples.

Groups ecDNA circles (decoil reconstruct.json, label=ecDNA) into families:
two circles from different samples are in the same family when they share
>= min_shared fragment-boundary breakpoints (same chromosome, within
--tolerance bp). For each family, reports the exact breakpoint coordinates,
which circles carry them, and the per-circle coordinate spread.

Outputs a TSV (per-family breakpoint table) and a summary to stdout.

Parameters
----------
--decoil-dir : directory containing <sample>/decoil/reconstruct.json
--samples : sample names (repeatable)
--tolerance : breakpoint match tolerance in bp (default 200)
--min-shared : minimum shared breakpoints to link two circles (default 2)
--out : output TSV path
"""
from __future__ import annotations

import argparse
import csv
import itertools
import logging
import os
from collections import defaultdict

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("ecdna_breakpoint_families")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--decoil-dir", required=True,
                        help="Directory containing <sample>/decoil/reconstruct.json")
    parser.add_argument("--samples", action="append", required=True,
                        help="Sample name, repeatable")
    parser.add_argument("--tolerance", type=int, default=200,
                        help="Breakpoint match tolerance in bp (default: 200)")
    parser.add_argument("--min-shared", type=int, default=2,
                        help="Min shared breakpoints to link two circles (default: 2)")
    parser.add_argument("--out", required=True, help="Output TSV path")
    return parser.parse_args(argv)


def load_circles(decoil_dir: str, samples: list[str]) -> dict:
    """Load {(sample, circ_id): [(chrom, start, end), ...]} from
    reconstruct.ecDNA.filtered.bed (the filtered, downstream-standard source).

    One bed row = one merged fragment; fragments of the same circ_id combine
    into one circle.
    """
    circles = {}
    for s in samples:
        path = os.path.join(decoil_dir, s, "decoil",
                            "reconstruct.ecDNA.filtered.bed")
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                circles.setdefault((s, f[3]), []).append(
                    (f[0], int(f[1]), int(f[2])))
    return circles


def edge_set(frags: list[tuple]) -> set:
    out = set()
    for chrom, start, end in frags:
        out.add((chrom, start))
        out.add((chrom, end))
    return out


def shared_breakpoints(e1: set, e2: set, tol: int) -> list:
    """Return [(chrom, x, y, dist)] breakpoint matches within tol bp."""
    hits = []
    for ch1, x in e1:
        for ch2, y in e2:
            if ch1 == ch2 and abs(x - y) <= tol:
                hits.append((ch1, x, y, abs(x - y)))
    return hits


class UnionFind:
    def __init__(self, keys):
        self.parent = {k: k for k in keys}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    circles = load_circles(args.decoil_dir, args.samples)
    logger.info("Loaded %d ecDNA circles from %d samples",
                len(circles), len(args.samples))

    uf = UnionFind(list(circles))
    n_links = 0
    for (s1, c1), (s2, c2) in itertools.combinations(circles, 2):
        if s1 == s2:
            continue
        hits = shared_breakpoints(edge_set(circles[(s1, c1)]),
                                  edge_set(circles[(s2, c2)]), args.tolerance)
        # distinct breakpoint positions (collapse mirror pairs of the same edge)
        positions = {h[0] for h in hits}
        if len(hits) >= args.min_shared:
            uf.union((s1, c1), (s2, c2))
            n_links += 1
    logger.info("%d cross-sample circle pairs linked (>= %d shared breakpoints)",
                n_links, args.min_shared)

    families = defaultdict(list)
    for key in circles:
        families[uf.find(key)].append(key)
    ordered = sorted(families.values(),
                     key=lambda m: (-len(m), min(m)))

    rows = []
    for fam_i, members in enumerate(ordered, 1):
        members = sorted(members)
        loci = {(frags[0][0], frags[0][1] // 1_000_000)
                for frags in (circles[key] for key in members)}
        locus_str = ";".join(f"{ch}:{mb}Mb" for ch, mb in sorted(loci))
        # breakpoint columns: positions shared by >= 2 members
        edge_sets = {f"{s}:circ{c}": edge_set(circles[(s, c)])
                     for s, c in members}
        all_edges = set().union(*edge_sets.values())
        # cluster matched edges into distinct breakpoint positions (edges of the
        # same locus within tolerance collapse to one breakpoint)
        bps = []
        used = set()
        for chrom, pos in sorted(all_edges):
            if (chrom, pos) in used:
                continue
            hits = {}
            for mkey, es in edge_sets.items():
                for (ch2, y) in es:
                    if ch2 == chrom and abs(pos - y) <= args.tolerance:
                        hits[mkey] = y
            if len(hits) >= 2:
                spread = max(hits.values()) - min(hits.values())
                canonical = min(hits.values())
                bps.append((chrom, canonical, len(hits), spread, hits))
                # mark every edge in the cluster as used
                for (ch2, y) in all_edges:
                    if ch2 == chrom and abs(pos - y) <= args.tolerance:
                        used.add((ch2, y))
        nsamples = len({s for s, _ in members})
        summary = (f"family{fam_i} n={len(members)} samples={nsamples} "
                   f"locus={locus_str} bps={len(bps)}")
        logger.info(summary)
        rows.append({
            "family": fam_i, "n_circles": len(members), "n_samples": nsamples,
            "locus": locus_str, "members": ";".join(f"{s}:circ{c}" for s, c in members),
            "n_shared_bps": len(bps),
            "bp_detail": " | ".join(
                f"{ch}:{pos:,}(n={n},spread={spread}bp)"
                for ch, pos, n, spread, _ in bps),
            "singletons": "" if len(members) > 1 else members[0][0],
        })

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    logger.info("Wrote %s", args.out)

    # stdout: per-family detail with member lists
    for fam_i, members in enumerate(ordered, 1):
        label = ";".join(f"{s}:circ{c}" for s, c in sorted(members))
        print(f"family {fam_i}: {len(members)} circles, "
              f"{len({s for s, _ in members})} samples -> {label}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
