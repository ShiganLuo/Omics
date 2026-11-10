#!/usr/bin/env python3
"""Repeat-element enrichment across eccDNA bodies.

Question: which repetitive elements are enriched on eccDNA overall?  The
eccDNA region set is the merged union of all decoil
``reconstruct.ecDNA.filtered.bed`` fragments (identical loci shared across
samples count once).  For every repeat class and family the script computes
the fraction of eccDNA bases covered and compares it against a matched
background: each draw re-places every eccDNA interval at a random position
of the same chromosome and length, 1000 draws.

Outputs: per-class and per-family TSVs (observed vs background, enrichment,
empirical p, BH-FDR), the merged region table, and a stdout summary.

Parameters
----------
--ecdna-bed : ecDNA fragment BED as ``SAMPLE=PATH`` (repeatable; the file may
    live anywhere; typically ``reconstruct.ecDNA.filtered.bed``)
--rmsk : UCSC rmsk GTF (plain or gzipped, family_id/class_id attributes)
--chrom-sizes : chromosome sizes (two-column TSV, e.g. genome .fai)
--n-draws : background placements per interval (default 1000)
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
logger = logging.getLogger("ecdna_repeat_enrichment")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--ecdna-bed", action="append", required=True,
                        type=sample_path_spec, metavar="SAMPLE=PATH",
                        help="ecDNA fragment BED (repeatable; any file location)")
    parser.add_argument("--rmsk", required=True,
                        help="UCSC rmsk GTF (plain or gzipped)")
    parser.add_argument("--chrom-sizes", required=True,
                        help="Two-column chromosome sizes TSV")
    parser.add_argument("--n-draws", type=int, default=1000,
                        help="Background placements per interval (default: 1000)")
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


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1] = (out[-1][0], e)
        else:
            out.append((s, e))
    return out


def eccdna_regions(fragments, per_sample: bool = False
                   ) -> dict[str, list[tuple[str, int, int]]]:
    """Merged eccDNA regions as ``{label: [(chrom, start, end), ...]}``."""
    if per_sample:
        out = {}
        for sample, rows in fragments.items():
            by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
            for fr in rows:
                by_chrom[fr.chrom].append((fr.start, fr.end))
            out[sample] = [(c, s, e) for c, iv in by_chrom.items()
                           for s, e in merge_intervals(iv)]
        return out
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for rows in fragments.values():
        for fr in rows:
            by_chrom[fr.chrom].append((fr.start, fr.end))
    return {"all": [(c, s, e) for c, iv in by_chrom.items()
                    for s, e in merge_intervals(iv)]}


def covered_bp(rmsk: RmskIndex, regions: list[tuple[str, int, int]],
               key: RmskKey) -> int:
    return sum(rmsk.overlap_bp(c, s, e, key) for c, s, e in regions)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)

    fragments = load_fragments(args.ecdna_bed)
    regions_by_label = eccdna_regions(fragments)
    regions = regions_by_label["all"]
    total_bp = sum(e - s for _, s, e in regions)
    logger.info(f"eccDNA region set: {len(regions)} merged intervals, "
                f"{total_bp / 1e6:.2f} Mb")

    region_path = os.path.join(args.out_dir, "ecdna_regions.tsv")
    with open(region_path, "w") as fh:
        fh.write("chrom\tstart\tend\tlength\n")
        for c, s, e in sorted(regions):
            fh.write(f"{c}\t{s}\t{e}\t{e - s}\n")
    logger.info(f"Wrote {region_path}")

    rmsk = RmskIndex(args.rmsk, cache_dir=args.cache_dir)
    sizes = read_chrom_sizes(args.chrom_sizes)
    genome_total = sum(sizes.values())

    keys: list[RmskKey] = ["any"]
    key_name: dict[RmskKey, str] = {"any": "any"}
    for i, c in enumerate(rmsk.classes):
        keys.append(("cls", i))
        key_name[("cls", i)] = f"class:{c}"
    for i, f in enumerate(rmsk.families):
        keys.append(("fam", i))
        key_name[("fam", i)] = f

    # observed + genome context
    obs: dict[RmskKey, int] = {}
    genome_bp: dict[RmskKey, int] = {}
    for k in keys:
        obs[k] = covered_bp(rmsk, regions, k)
        genome_bp[k] = rmsk.family_bp(k)
    logger.info("Observed coverage computed")

    # background: matched random placements, per draw aggregate coverage
    draws_bp: dict[RmskKey, list[int]] = {k: [] for k in keys}
    for draw in range(args.n_draws):
        per_key = {k: 0 for k in keys}
        for c, s, e in regions:
            length = e - s
            lo = 0
            hi = max(0, sizes.get(c, length) - length)
            rs = rng.randint(lo, hi) if hi > 0 else 0
            for k in keys:
                per_key[k] += rmsk.overlap_bp(c, rs, rs + length, k)
        for k in keys:
            draws_bp[k].append(per_key[k])
        if (draw + 1) % 100 == 0:
            logger.info(f"background draws: {draw + 1}/{args.n_draws}")

    # per-sample 'any' coverage vs the same background logic (secondary)
    per_sample = eccdna_regions(fragments, per_sample=True)

    def write_table(path: str, want_level: str) -> None:
        rows = []
        for k in keys:
            if want_level == "class" and k != "any" and k[0] != "cls":
                continue
            if want_level == "family" and k[0] != "fam":
                continue
            obs_frac = obs[k] / total_bp
            bg_fracs = [b / total_bp for b in draws_bp[k]]
            mean_bg = statistics.mean(bg_fracs)
            enr = obs_frac / mean_bg if mean_bg > 0 else float("inf")
            p = empirical_p(obs_frac, bg_fracs)
            p_low = empirical_p(obs_frac, bg_fracs, greater=False)
            rows.append((k, obs_frac, mean_bg, enr, p, p_low))
        fdrs = bh_fdr([r[4] for r in rows])
        with open(path, "w") as fh:
            fh.write("key\tclass\tecdna_bp\tecdna_frac\tbg_frac_mean\t"
                     "enrichment\tp_emp\tp_emp_low\tp_fdr\tgenome_frac\t"
                     "enrichment_vs_genome\n")
            for (k, obs_frac, mean_bg, enr, p, p_low), fdr in zip(rows, fdrs):
                cls = ("all" if k == "any"
                       else rmsk.classes[rmsk.family_class[k[1]]]
                       if k[0] == "fam" else rmsk.classes[k[1]])
                gf = genome_bp[k] / genome_total
                fh.write(f"{key_name[k]}\t{cls}\t{obs[k]}\t{obs_frac:.5f}\t"
                         f"{mean_bg:.5f}\t{enr:.3f}\t{p:.4f}\t{p_low:.4f}\t"
                         f"{fdr:.4f}\t{gf:.5f}\t"
                         f"{obs_frac / gf if gf > 0 else float('inf'):.3f}\n")
        logger.info(f"Wrote {path}")

    class_path = os.path.join(args.out_dir, "ecdna_repeat_enrichment_class.tsv")
    fam_path = os.path.join(args.out_dir, "ecdna_repeat_enrichment_family.tsv")
    write_table(class_path, "class")
    write_table(fam_path, "family")

    # secondary: per-sample any-repeat fraction
    sample_path = os.path.join(args.out_dir, "ecdna_repeat_per_sample.tsv")
    with open(sample_path, "w") as fh:
        fh.write("sample\tecdna_bp\tfrac_any\tfrac_class_LINE\tfrac_class_SINE\t"
                 "frac_class_LTR\tfrac_class_DNA\n")
        cls_of = {c: ("cls", i) for i, c in enumerate(rmsk.classes)}
        for label, regs in per_sample.items():
            tbp = sum(e - s for _, s, e in regs)
            fracs = [covered_bp(rmsk, regs, "any") / tbp]
            for cname in ("LINE", "SINE", "LTR", "DNA"):
                k = cls_of.get(cname)
                fracs.append(covered_bp(rmsk, regs, k) / tbp if k else 0.0)
            fh.write(f"{label}\t{tbp}\t" + "\t".join(f"{f:.4f}" for f in fracs)
                     + "\n")
    logger.info(f"Wrote {sample_path}")

    # ---------------- stdout summary ---------------- #
    print(f"\n=== eccDNA repeat enrichment ({len(regions)} intervals, "
          f"{total_bp / 1e6:.2f} Mb) ===")
    print(f"{'key':<20s} {'ecdna':>8s} {'bg':>8s} {'enr':>6s} {'p':>7s} "
          f"{'fdr':>7s}")
    class_rows = []
    for k in keys:
        if k != "any" and k[0] != "cls":
            continue
        obs_frac = obs[k] / total_bp
        mean_bg = statistics.mean(b / total_bp for b in draws_bp[k])
        enr = obs_frac / mean_bg if mean_bg > 0 else float("inf")
        p = empirical_p(obs_frac, [b / total_bp for b in draws_bp[k]])
        class_rows.append((key_name[k], obs_frac, mean_bg, enr, p))
    for fdr, (name, obs_frac, mean_bg, enr, p) in zip(
            bh_fdr([r[4] for r in class_rows]), class_rows):
        print(f"{name:<20s} {obs_frac:8.4f} {mean_bg:8.4f} {enr:6.2f} "
              f"{p:7.4f} {fdr:7.4f}")
    fam_rows = []
    for k in keys:
        if k[0] != "fam":
            continue
        obs_frac = obs[k] / total_bp
        mean_bg = statistics.mean(b / total_bp for b in draws_bp[k])
        if mean_bg == 0 or obs_frac == 0:
            continue
        p = empirical_p(obs_frac, [b / total_bp for b in draws_bp[k]])
        fam_rows.append((key_name[k], obs_frac, mean_bg,
                         obs_frac / mean_bg, p))
    fam_rows.sort(key=lambda r: -r[3])
    print("top enriched families (matched-background enrichment):")
    for name, obs_frac, mean_bg, enr, p in fam_rows[:15]:
        print(f"  {name:<16s} obs {obs_frac:.4f}  bg {mean_bg:.4f}  "
              f"enr {enr:.2f}  p={p:.4f}")
    fam_rows.sort(key=lambda r: r[3])
    print("top depleted families:")
    for name, obs_frac, mean_bg, enr, p in fam_rows[:8]:
        print(f"  {name:<16s} obs {obs_frac:.4f}  bg {mean_bg:.4f}  "
              f"enr {enr:.2f}  p={p:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
