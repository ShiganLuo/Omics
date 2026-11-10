#!/usr/bin/env python3
"""Repeat content at eccDNA breakpoints.

Question: do eccDNA breakpoints carry repetitive sequences, or are they
generally repeat-free?  Breakpoints are the fragment boundaries of decoil
reconstruct.ecDNA.filtered.bed (cross-sample coincident breakpoints collapse
into one locus).  For every breakpoint locus the script reports

- repeat-element hit rate in windows around the breakpoint (±0/100/... bp),
- distance to the nearest repeat element,
- repeat family composition in a flanking window,

each compared against three null models (1000 position draws per locus):

- ``local``          — random position in the locus's own circle span ± flank
                       (controls regional repeat density),
- ``circle_interior``— random position inside the circle's own fragments
                       (the DNA that was actually retained on the circle),
- ``genome``         — random position on the same chromosome.

Outputs (TSV) plus a stdout summary.

Parameters
----------
--ecdna-bed : ecDNA fragment BED as ``SAMPLE=PATH`` (repeatable; the file may
    live anywhere; typically ``reconstruct.ecDNA.filtered.bed``)
--rmsk : UCSC rmsk GTF (plain or gzipped, family_id/class_id attributes)
--chrom-sizes : chromosome sizes (two-column TSV, e.g. genome .fai)
--windows : comma-separated half-window sizes in bp
    (default 0,25,50,100,250,500,1000,2000,5000)
--comp-window : half-window for family composition (default 500)
--tolerance : cross-sample breakpoint collapse distance in bp (default 200)
--local-flank : flank added around the circle span for the local null (default 50000)
--n-draws : permutations per null model (default 1000)
--seed : RNG seed (default 0)
--cache-dir : optional pickle cache dir for the parsed rmsk index
--out-dir : output directory
"""
from __future__ import annotations

import argparse
import logging
import os
import random
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repeat_utils import (RmskIndex, RmskKey, bh_fdr, empirical_p,  # noqa: E402
                          load_fragments, sample_map, sample_path_spec)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("ecdna_breakpoint_repeat")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--ecdna-bed", action="append", required=True,
                        type=sample_path_spec, metavar="SAMPLE=PATH",
                        help="ecDNA fragment BED (repeatable; any file location)")
    parser.add_argument("--rmsk", required=True,
                        help="UCSC rmsk GTF (plain or gzipped)")
    parser.add_argument("--chrom-sizes", required=True,
                        help="Two-column chromosome sizes TSV")
    parser.add_argument("--windows", default="0,25,50,100,250,500,1000,2000,5000",
                        help="Comma-separated half-window sizes (bp)")
    parser.add_argument("--comp-window", type=int, default=500,
                        help="Half-window for family composition (default: 500)")
    parser.add_argument("--tolerance", type=int, default=200,
                        help="Breakpoint collapse distance (default: 200)")
    parser.add_argument("--local-flank", type=int, default=50_000,
                        help="Flank for the local null (default: 50000)")
    parser.add_argument("--n-draws", type=int, default=1000,
                        help="Permutations per null model (default: 1000)")
    parser.add_argument("--seed", type=int, default=0,
                        help="RNG seed (default: 0)")
    parser.add_argument("--cache-dir", default=None,
                        help="Optional pickle cache dir for the rmsk index")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    args = parser.parse_args(argv)
    try:
        args.ecdna_bed = sample_map(args.ecdna_bed)
    except ValueError as exc:
        parser.error(str(exc))
    return args


def read_chrom_sizes(path: str) -> dict[str, int]:
    sizes = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            sizes[f[0]] = int(f[1])
    return sizes


def collect_breakpoints(fragments) -> list[dict]:
    """One record per fragment boundary: sample, circ, side, pos."""
    recs = []
    for sample, rows in fragments.items():
        for fr in rows:
            recs.append(dict(sample=sample, circ=fr.circ_id,
                             fragment=fr.fragment_id, side="start",
                             chrom=fr.chrom, pos=fr.start))
            recs.append(dict(sample=sample, circ=fr.circ_id,
                             fragment=fr.fragment_id, side="end",
                             chrom=fr.chrom, pos=fr.end))
    return recs


def cluster_loci(recs: list[dict], tolerance: int) -> list[dict]:
    """Collapse cross-sample breakpoints within ``tolerance`` bp into loci."""
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for r in recs:
        by_chrom[r["chrom"]].append(r)
    loci = []
    for chrom, rs in by_chrom.items():
        rs.sort(key=lambda r: r["pos"])
        group = [rs[0]]
        for r in rs[1:]:
            if r["pos"] - group[-1]["pos"] <= tolerance:
                group.append(r)
            else:
                loci.append(_make_locus(chrom, group))
                group = [r]
        loci.append(_make_locus(chrom, group))
    loci.sort(key=lambda d: (d["chrom"], d["pos"]))
    for i, loc in enumerate(loci, 1):
        loc["locus_id"] = f"BP{i:03d}"
    return loci


def _make_locus(chrom: str, group: list[dict]) -> dict:
    return dict(chrom=chrom,
                pos=int(statistics.median(r["pos"] for r in group)),
                members=sorted(f"{r['sample']}:{r['circ']}:{r['fragment']}:{r['side']}"
                               for r in group),
                n_breakpoints=len(group))


def circle_spans(fragments) -> dict[tuple[str, str], tuple[str, int, int]]:
    """{(sample, circ): (chrom, span_start, span_end)} from fragment rows."""
    spans = {}
    for sample, rows in fragments.items():
        for fr in rows:
            key = (sample, fr.circ_id)
            if key not in spans:
                spans[key] = (fr.chrom, fr.start, fr.end)
            else:
                c, s, e = spans[key]
                spans[key] = (c, min(s, fr.start), max(e, fr.end))
    return spans


def fragment_union(fragments, sample: str, circ: str) -> list[tuple[int, int]]:
    """Merged fragment intervals of one circle (sample, circ)."""
    iv = sorted((fr.start, fr.end) for fr in fragments[sample]
                if fr.circ_id == circ)
    out: list[tuple[int, int]] = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def random_position(rng, chrom: str, sizes: dict[str, int],
                    spans, circle_iv, local_flank: int, mode: str) -> int:
    """Draw one random breakpoint position under the given null model."""
    if mode == "genome":
        return rng.randrange(sizes[chrom])
    if mode == "local":
        _, s, e = spans
        lo = max(0, s - local_flank)
        hi = min(sizes[chrom], e + local_flank)
        return rng.randrange(lo, hi)
    # circle_interior: pick a fragment interval weighted by length, then a pos
    total = sum(e - s for s, e in circle_iv)
    x = rng.randrange(total)
    for s, e in circle_iv:
        if x < e - s:
            return s + x
        x -= e - s
    return circle_iv[-1][1] - 1


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)
    windows = sorted(int(w) for w in args.windows.split(","))

    fragments = load_fragments(args.ecdna_bed)
    recs = collect_breakpoints(fragments)
    loci = cluster_loci(recs, args.tolerance)
    logger.info(f"{len(recs)} fragment boundaries -> {len(loci)} breakpoint "
                f"loci (tolerance {args.tolerance} bp)")

    rmsk = RmskIndex(args.rmsk, cache_dir=args.cache_dir)
    sizes = read_chrom_sizes(args.chrom_sizes)
    spans = circle_spans(fragments)

    keys: list[RmskKey] = ["any"] + [("cls", i) for i in range(len(rmsk.classes))]
    key_name: dict[RmskKey, str] = {
        ("cls", i): f"class:{c}" for i, c in enumerate(rmsk.classes)}
    key_name["any"] = "any"

    # per-locus circle context (first member decides the circle span / union)
    for loc in loci:
        sample, circ = loc["members"][0].split(":")[:2]
        loc["span"] = spans[(sample, circ)]
        loc["circle_iv"] = fragment_union(fragments, sample, circ)

    # ---------------- observed statistics ---------------- #
    def window_hits(pos: int, chrom: str) -> dict[tuple[int, RmskKey], bool]:
        out: dict[tuple[int, RmskKey], bool] = {}
        for w in windows:
            for k in keys:
                if w == 0:
                    out[(w, k)] = rmsk.contains(chrom, pos, k)
                else:
                    out[(w, k)] = rmsk.overlap_bp(chrom, pos - w, pos + w, k) > 0
        return out

    obs_hits: list[dict[tuple[int, RmskKey], bool]] = [
        window_hits(loc["pos"], loc["chrom"]) for loc in loci]
    obs_nearest = [rmsk.nearest_distance(loc["chrom"], loc["pos"]) for loc in loci]

    comp_fams: dict[int, set[int]] = {}   # locus index -> set of fam idx
    for i, loc in enumerate(loci):
        els = rmsk.elements_overlapping(loc["chrom"], loc["pos"] - args.comp_window,
                                        loc["pos"] + args.comp_window)
        comp_fams[i] = {fi for _, _, fi, _ in els}

    # ---------------- null draws ---------------- #
    null_hits: dict[str, dict[int, dict[RmskKey, list[float]]]] = {
        mode: {w: {k: [] for k in keys} for w in windows}
        for mode in ("local", "circle_interior", "genome")}
    null_comp: dict[str, dict[int, list[float]]] = {
        mode: defaultdict(list) for mode in null_hits}
    null_nearest: dict[str, list[float]] = {mode: [] for mode in null_hits}

    modes = ("local", "circle_interior", "genome")
    for draw in range(args.n_draws):
        draw_hits: dict[str, dict[int, dict[RmskKey, int]]] = {
            mode: {w: {k: 0 for k in keys} for w in windows} for mode in modes}
        draw_nearest: dict[str, list[int]] = {mode: [] for mode in modes}
        draw_comp: dict[str, dict[int, int]] = {
            mode: defaultdict(int) for mode in modes}
        for i, loc in enumerate(loci):
            for mode in modes:
                p = random_position(rng, loc["chrom"], sizes, loc["span"],
                                    loc["circle_iv"], args.local_flank, mode)
                h = window_hits(p, loc["chrom"])
                for (w, k), v in h.items():
                    draw_hits[mode][w][k] += int(v)
                nd = rmsk.nearest_distance(loc["chrom"], p)
                if nd is not None:
                    draw_nearest[mode].append(nd)
                els = rmsk.elements_overlapping(
                    loc["chrom"], p - args.comp_window, p + args.comp_window)
                for fi in {fi for _, _, fi, _ in els}:
                    draw_comp[mode][fi] += 1
        for mode in modes:
            for w in windows:
                for k in keys:
                    null_hits[mode][w][k].append(draw_hits[mode][w][k] / len(loci))
            null_nearest[mode].append(
                statistics.mean(draw_nearest[mode]) if draw_nearest[mode] else float("nan"))
            for fi, cnt in draw_comp[mode].items():
                null_comp[mode][fi].append(cnt / len(loci))
        if (draw + 1) % 200 == 0:
            logger.info(f"null draws: {draw + 1}/{args.n_draws}")

    # ---------------- output: per-locus table ---------------- #
    loci_path = os.path.join(args.out_dir, "breakpoint_loci.tsv")
    with open(loci_path, "w") as fh:
        cols = (["locus_id", "chrom", "pos", "n_breakpoints",
                 "nearest_repeat_bp", "members"]
                + [f"hit_{key_name[k]}_w{w}" for w in windows for k in keys])
        fh.write("\t".join(cols) + "\n")
        for i, loc in enumerate(loci):
            row = [loc["locus_id"], loc["chrom"], str(loc["pos"]),
                   str(loc["n_breakpoints"]),
                   "" if obs_nearest[i] is None else str(obs_nearest[i]),
                   ",".join(loc["members"])]
            row += [str(int(obs_hits[i][(w, k)])) for w in windows for k in keys]
            fh.write("\t".join(row) + "\n")
    logger.info(f"Wrote {loci_path}")

    # ---------------- output: window summary ---------------- #
    window_path = os.path.join(args.out_dir, "breakpoint_window_summary.tsv")
    with open(window_path, "w") as fh:
        fh.write("window_bp\tkey\tobs_frac\tn_loci_hit\t"
                 "null\tnull_mean_frac\tenrichment\tp_emp\tp_emp_low\n")
        for w in windows:
            for k in keys:
                obs_frac = sum(h[(w, k)] for h in obs_hits) / len(loci)
                n_hit = sum(h[(w, k)] for h in obs_hits)
                for mode in modes:
                    draws = null_hits[mode][w][k]
                    mean_null = statistics.mean(draws)
                    enr = obs_frac / mean_null if mean_null > 0 else float("inf")
                    p = empirical_p(obs_frac, draws)
                    p_low = empirical_p(obs_frac, draws, greater=False)
                    fh.write(f"{w}\t{key_name[k]}\t{obs_frac:.4f}\t{n_hit}\t"
                             f"{mode}\t{mean_null:.4f}\t{enr:.3f}\t{p:.4f}\t"
                             f"{p_low:.4f}\n")
    logger.info(f"Wrote {window_path}")

    # ---------------- output: family composition ---------------- #
    fam_path = os.path.join(args.out_dir, "breakpoint_family_composition.tsv")
    obs_fam_frac = {fi: sum(1 for i in range(len(loci)) if fi in comp_fams[i])
                    / len(loci) for fi in range(len(rmsk.families))}
    with open(fam_path, "w") as fh:
        fh.write("family\tclass\tobs_locus_frac\tn_loci_hit\t"
                 "null\tnull_mean_frac\tenrichment\tp_emp\tp_fdr\n")
        rows: dict[str, list[list]] = {mode: [] for mode in modes}
        for fi, fam in enumerate(rmsk.families):
            obs_f = obs_fam_frac[fi]
            if obs_f == 0:
                continue
            for mode in modes:
                draws = null_comp[mode].get(fi, [])
                mean_null = statistics.mean(draws) if draws else 0.0
                enr = obs_f / mean_null if mean_null > 0 else float("inf")
                p = empirical_p(obs_f, draws) if draws else 1.0
                rows[mode].append([fam, rmsk.classes[rmsk.family_class[fi]],
                                   obs_f, int(obs_f * len(loci) + 0.5),
                                   mode, mean_null, enr, p])
        for mode in modes:
            fdrs = bh_fdr([r[-1] for r in rows[mode]])
            for r, fdr in zip(rows[mode], fdrs):
                fh.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]}\t{r[4]}\t"
                         f"{r[5]:.4f}\t{r[6]:.3f}\t{r[7]:.4f}\t{fdr:.4f}\n")
    logger.info(f"Wrote {fam_path}")

    # ---------------- stdout summary ---------------- #
    print(f"\n=== eccDNA breakpoint repeat content ({len(loci)} loci) ===")
    print(f"breakpoints with a repeat within W bp (obs vs null, any repeat):")
    for w in windows:
        obs_frac = sum(h[(w, "any")] for h in obs_hits) / len(loci)
        null_txt = "  ".join(
            f"{m}: {statistics.mean(null_hits[m][w]['any']):.2f}"
            for m in modes)
        print(f"  W=±{w:<5d} obs {obs_frac:.2f}   [{null_txt}]")
    inside = sum(1 for d in obs_nearest if d == 0)
    finite = [d for d in obs_nearest if d is not None]
    print(f"breakpoint inside a repeat element: {inside}/{len(loci)}")
    print(f"median nearest-repeat distance: {statistics.median(finite):.0f} bp")
    fam_rows = []
    for fi in range(len(rmsk.families)):
        obs_f = obs_fam_frac[fi]
        if obs_f == 0:
            continue
        draws = null_comp["local"].get(fi, [])
        if draws and statistics.mean(draws) > 0:
            fam_rows.append((rmsk.families[fi], obs_f,
                             statistics.mean(draws),
                             empirical_p(obs_f, draws)))
    fam_rows.sort(key=lambda r: -r[1])
    print(f"top repeat families within ±{args.comp_window} bp "
          f"(locus fraction, vs local null):")
    for fam, obs_f, null_f, p in fam_rows[:12]:
        print(f"  {fam:<16s} obs {obs_f:.2f}  null {null_f:.2f}  "
              f"enr {obs_f / null_f if null_f else float('inf'):.2f}  p={p:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
