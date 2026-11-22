#!/usr/bin/env python3
"""Epigenetic drift analysis (fig9 data): per-fiber m5CpG fraction variance,
group A vs group B, over random windows (subsampled fibers per side).

Port of conclusion_figures/plot_fig9_drift.py data half to a reusable module
script: group labels come from --group instead of hardcoded youth/aged.

Per window: fibers covering >=50% of the window give a m5CpG fraction
(m5CpG marks / covered bp); groups with >= --min-fibers fibers are kept and
subsampled to --subsample fibers per side; per-group mean and variance are
reported.  Drift = variance of the per-fiber m5CpG fraction.

RNG: window draws use one global RNG from --seed; the per-window fiber
subsample uses a window-keyed RNG so results do not depend on iteration
order (the original script shuffled with a global RNG — same method,
order-independent reproducibility).

Outputs:
  <windows-out>   chrom start end mean_<A> mean_<B> var_<A> var_<B>
  <summary-out>   JSON: groups, n, signs, sign_test_p, median_var_ratio

Usage:
  python epi_drift_run.py \
      --bam youth=luanchao-21310.fiberseq.fire.sorted.bam \
      --bam aged=luanchao-11236.fiberseq.fire.sorted.bam \
      --group youth=youth --group aged=aged \
      --fai genome.fa.fai \
      --windows-out out/epigenetic_drift_windows.tsv \
      --summary-out out/epigenetic_drift.json
"""
import argparse
import json
import math
import os
import random
import re
import statistics

import pysam


def mm_positions(entry, seq):
    """MM tag entry -> read positions of modified bases.

    Self-contained copy (modules.md #0: 模块脚本自包含, 不引用其他模块代码).
    """
    m = re.match(r"([ACGTN])([+-])([a-z]+)(\.?)", entry)
    if not m:
        return []
    base = m.group(1)
    nums = entry[m.end():].lstrip(",?")
    try:
        deltas = [int(x) for x in nums.split(",") if x != ""]
    except ValueError:
        return []
    idx = [i for i, c in enumerate(seq) if c.upper() == base]
    out, cur = [], -1
    for d in deltas:
        cur = d if not out else cur + 1 + d
        if 0 <= cur < len(idx):
            out.append(idx[cur])
        else:
            break
    return out


def parse_kv(items, what):
    out = {}
    for it in items or []:
        if "=" not in it:
            raise ValueError(f"--{what} expects KEY=VALUE, got: {it!r}")
        k, v = it.split("=", 1)
        out[k] = v
    return out


def parse_groups(items):
    groups = {}
    for it in items or []:
        if "=" not in it:
            raise ValueError(f"--group expects LABEL=S1,S2, got: {it!r}")
        label, s = it.split("=", 1)
        groups[label] = [x for x in s.split(",") if x]
    if len(groups) != 2:
        raise ValueError(f"exactly 2 --group entries required (got {len(groups)})")
    return groups


def fiber_fracs(bams, mm_positions, chrom, s, e):
    """Per-fiber m5CpG fraction over [s, e), pooled over the group's BAMs."""
    out = []
    for path in bams:
        af = pysam.AlignmentFile(path, "rb")
        for r in af.fetch(chrom, s, e):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            cov = 0
            marks = 0
            ref_of_read = {}
            ref, read = r.reference_start, 0
            for op, ln in (r.cigartuples or []):
                if op in (0, 7, 8):
                    for i in range(ln):
                        ref_of_read[read + i] = ref + i
                    read += ln
                    ref += ln
                elif op in (2, 3):
                    ref += ln
                else:
                    read += ln
            for a, b in r.get_blocks():
                cov += max(0, min(b, e) - max(a, s))
            if cov < 0.5 * (e - s):
                continue
            seq = r.query_sequence or ""
            if r.has_tag("MM"):
                for part in str(r.get_tag("MM")).split(";"):
                    if part[:2] in ("C+", "G-"):
                        for p in mm_positions(part, seq):
                            x = ref_of_read.get(p)
                            if x is not None and s <= x < e:
                                marks += 1
            out.append(marks / cov if cov else 0.0)
        af.close()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", action="append", required=True,
                    help="SAMPLE=PATH aligned fiberseq BAM (repeatable)")
    ap.add_argument("--group", action="append", required=True,
                    help="LABEL=S1,S2 (repeatable; exactly 2 groups)")
    ap.add_argument("--fai", required=True, help="reference .fai index")
    ap.add_argument("--windows-out", required=True)
    ap.add_argument("--summary-out", required=True)
    ap.add_argument("--win", type=int, default=5000)
    ap.add_argument("--min-fibers", type=int, default=20)
    ap.add_argument("--subsample", type=int, default=20)
    ap.add_argument("--n-windows", type=int, default=600)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-chrom-len", type=int, default=30000000)
    ap.add_argument("--chrom-prefix", default="",
                    help="only sample windows on chroms starting with this "
                         "(e.g. NC_ for NCBI accessions)")
    args = ap.parse_args()

    bams = parse_kv(args.bam, "bam")
    groups = parse_groups(args.group)
    labels = list(groups)
    label_a, label_b = labels
    for label, slist in groups.items():
        missing = [s for s in slist if s not in bams]
        if missing:
            raise ValueError(f"group {label!r} references samples without "
                             f"--bam: {missing}")
    group_bams = {label: [bams[s] for s in groups[label]] for label in labels}

    fai = {}
    for line in open(args.fai):
        f = line.split("\t")
        fai[f[0]] = int(f[1])
    chroms = [c for c in sorted(fai)
              if c.startswith(args.chrom_prefix) and fai[c] > args.min_chrom_len]
    if not chroms:
        raise ValueError(
            f"no chromosomes pass --chrom-prefix {args.chrom_prefix!r} + "
            f"--min-chrom-len {args.min_chrom_len}; check the reference")

    rng = random.Random(args.seed)
    windows = []
    for _ in range(args.n_windows):
        c = rng.choice(chroms)
        s = rng.randint(1000000, fai[c] - args.win - 1000000)
        windows.append((c, s, s + args.win))
    print(f"{len(windows)} candidate windows on {len(chroms)} chroms", flush=True)

    rows = []
    for chrom, s, e in windows:
        fa = fiber_fracs(group_bams[label_a], mm_positions, chrom, s, e)
        fb = fiber_fracs(group_bams[label_b], mm_positions, chrom, s, e)
        if len(fa) < args.min_fibers or len(fb) < args.min_fibers:
            continue
        wrng = random.Random(f"{args.seed}:{chrom}:{s}")
        wrng.shuffle(fa)
        wrng.shuffle(fb)
        fa = fa[:args.subsample]
        fb = fb[:args.subsample]
        va, vb = statistics.variance(fa), statistics.variance(fb)
        ma, mb = statistics.mean(fa), statistics.mean(fb)
        if va == 0 and vb == 0:
            continue
        rows.append((chrom, s, e, ma, mb, va, vb))

    for out_path in (args.windows_out, args.summary_out):
        out_dir = os.path.dirname(out_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
    with open(args.windows_out, "w") as fh:
        fh.write(f"chrom\tstart\tend\tmean_{label_a}\tmean_{label_b}"
                 f"\tvar_{label_a}\tvar_{label_b}\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    ratios = [(r[6] + 1e-12) / (r[5] + 1e-12) for r in rows]
    signs = sum(1 for r in rows if r[6] > r[5])
    n = len(rows)
    if n == 0:
        raise RuntimeError("no windows passed the min-fibers filter")
    k = min(signs, n - signs)
    p = 2 * sum(math.comb(n, i) for i in range(0, k + 1)) / 2 ** n
    p = min(1.0, p)

    summary = {
        "groups": labels,
        "n": n,
        "signs": signs,
        "sign_test_p": p,
        "median_var_ratio": statistics.median(ratios),
        "params": {
            "win": args.win, "min_fibers": args.min_fibers,
            "subsample": args.subsample, "n_windows": args.n_windows,
            "seed": args.seed, "min_chrom_len": args.min_chrom_len,
            "chrom_prefix": args.chrom_prefix,
        },
    }
    with open(args.summary_out, "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"windows={n} signs={signs} p={p:.3g} "
          f"median_var_ratio={summary['median_var_ratio']:.3f}")
    print(f"wrote {args.windows_out} + {args.summary_out}")


if __name__ == "__main__":
    main()
