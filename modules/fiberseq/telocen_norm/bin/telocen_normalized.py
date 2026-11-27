#!/usr/bin/env python3
"""Depth-normalized telocen metrics over telomere/centromere regions (fig2 data).

Port of conclusion_figures/compute_telocen_normalized_mp.py to a reusable
module script: no hardcoded paths, two comparison groups instead of fixed
youth/aged sample ids.

Per (sample, region): aligned_bp from get_blocks clipped to the region;
m5CpG and m6A marks from MM tags mapped to reference positions inside the
region.  Densities per aligned Mb.  Regions classified by BED column 4
(cen_core / cen_flank / subtelomere / telomere_core).

Outputs:
  <outdir>/region_normalized.tsv     per (region, sample) densities
  <outdir>/class_diff_normalized.tsv per (class, metric) group comparison,
                                     paired Wilcoxon over regions

Stats honesty: with n=1 individual per group the Wilcoxon p is a paired
region-level direction-consistency measure, NOT population significance.

Usage:
  python telocen_normalized.py \
      --regions Mmul_T2T.telocen_regions.bed \
      --bam youth=luanchao-21310.fiberseq.fire.sorted.bam \
      --bam aged=luanchao-11236.fiberseq.fire.sorted.bam \
      --group youth=youth --group aged=aged \
      --outdir out/ --threads 32
"""
import argparse
import csv
import math
import os
import re
import statistics
import sys
from multiprocessing import Pool

from scipy.stats import wilcoxon


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
    """--group LABEL=S1,S2 -> ordered dict label -> [samples]; exactly 2."""
    groups = {}
    for it in items or []:
        if "=" not in it:
            raise ValueError(f"--group expects LABEL=S1,S2, got: {it!r}")
        label, s = it.split("=", 1)
        groups[label] = [x for x in s.split(",") if x]
    if len(groups) != 2:
        raise ValueError(
            f"exactly 2 --group entries required (got {len(groups)}); "
            "the comparison is a paired two-group design")
    return groups


BAMS = {}


def _init_worker(bams):
    global BAMS
    BAMS = bams


def work(t):
    import pysam
    sample, chrom, s, e, cls = t
    af = _af.get(sample)
    if af is None:
        af = pysam.AlignmentFile(BAMS[sample], "rb")
        _af[sample] = af
    aligned = 0
    mcpg = 0
    m6a = 0
    for r in af.fetch(chrom, s, e):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        for a, b in r.get_blocks():
            lo, hi = max(a, s), min(b, e)
            if hi > lo:
                aligned += hi - lo
        if not r.has_tag("MM"):
            continue
        seq = r.query_sequence
        if not seq:
            continue
        q2r = dict(r.get_aligned_pairs(matches_only=True))
        for entry in r.get_tag("MM").split(";"):
            if len(entry) < 3:
                continue
            base, strand, mod = entry[0], entry[1], entry[2]
            if base == "C" and strand == "+" and mod in "mh":
                kind = "cpg"
            elif base == "G" and strand == "-" and mod in "mh":
                kind = "cpg"
            elif base == "A" and mod == "a":
                kind = "m6a"
            else:
                continue
            for qp in mm_positions(entry, seq):
                rp = q2r.get(qp)
                if rp is not None and s <= rp < e:
                    if kind == "cpg":
                        mcpg += 1
                    else:
                        m6a += 1
    return (cls, chrom, s, e, sample, aligned, mcpg, m6a)


_af = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", required=True, help="BED: chrom start end region_class")
    ap.add_argument("--bam", action="append", required=True,
                    help="SAMPLE=PATH aligned fiberseq BAM (repeatable)")
    ap.add_argument("--group", action="append", required=True,
                    help="LABEL=S1,S2 (repeatable; exactly 2 groups)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=32)
    args = ap.parse_args()

    bams = parse_kv(args.bam, "bam")
    groups = parse_groups(args.group)
    for label, slist in groups.items():
        missing = [s for s in slist if s not in bams]
        if missing:
            raise ValueError(f"group {label!r} references samples without "
                             f"--bam: {missing}")
    os.makedirs(args.outdir, exist_ok=True)

    regions = []
    with open(args.regions) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.split()
            if len(p) >= 4:
                regions.append((p[0], int(p[1]), int(p[2]), p[3]))
    print(f"{len(regions)} regions", flush=True)

    tasks = []
    for sample in bams:
        for chrom, s, e, cls in sorted(regions, key=lambda r: (r[0], r[1])):
            tasks.append((sample, chrom, s, e, cls))
    with Pool(processes=args.threads, initializer=_init_worker,
              initargs=(bams,)) as pool:
        rows = []
        done = 0
        for row in pool.imap_unordered(work, tasks, chunksize=1):
            rows.append(row)
            done += 1
            if done % 20 == 0:
                print(f"  {done}/{len(tasks)}", flush=True)

    sample_order = {}
    for label in groups:
        for s in groups[label]:
            sample_order.setdefault(s, len(sample_order))
    rows.sort(key=lambda r: (sample_order[r[4]], r[1], r[2]))

    region_tsv = os.path.join(args.outdir, "region_normalized.tsv")
    with open(region_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["region_class", "chrom", "start", "end", "sample",
                    "aligned_bp", "mcpg_marks", "m6a_marks",
                    "cpg_per_mb", "m6a_per_mb"])
        for cls, chrom, s, e, sample, aligned, mcpg, m6a in rows:
            w.writerow([cls, chrom, s, e, sample, aligned, mcpg, m6a,
                        f"{mcpg / aligned * 1e6:.4f}" if aligned else "NA",
                        f"{m6a / aligned * 1e6:.4f}" if aligned else "NA"])
    print("wrote region_normalized.tsv", flush=True)

    labels = list(groups)
    (label_a, samples_a), (label_b, samples_b) = \
        (labels[0], groups[labels[0]]), (labels[1], groups[labels[1]])

    by = {}
    for cls, chrom, s, e, sample, aligned, mcpg, m6a in rows:
        if aligned <= 0:
            continue
        key = (chrom, s, e)
        label = label_a if sample in samples_a else label_b
        for metric, val in (("cpg_per_mb", mcpg / aligned * 1e6),
                            ("m6a_per_mb", m6a / aligned * 1e6)):
            by.setdefault((cls, metric), {}).setdefault(label, {}) \
              .setdefault(key, []).append(val)

    out_rows = []
    for (cls, metric), per in sorted(by.items()):
        pa = {k: statistics.mean(v) for k, v in per.get(label_a, {}).items()}
        pb = {k: statistics.mean(v) for k, v in per.get(label_b, {}).items()}
        keys = sorted(set(pa) & set(pb))
        if not keys:
            continue
        av = [pa[k] for k in keys]
        bv = [pb[k] for k in keys]
        try:
            p = float(wilcoxon(av, bv).pvalue)
        except ValueError:
            p = 1.0
        ma, mb = statistics.mean(av), statistics.mean(bv)
        sa = statistics.stdev(av) if len(av) > 1 else 0.0
        sb = statistics.stdev(bv) if len(bv) > 1 else 0.0
        lfc = math.log2(mb / ma) if ma > 0 and mb > 0 else 0.0
        out_rows.append([cls, metric, len(keys), f"{ma:.6g}", f"{sa:.6g}",
                         f"{mb:.6g}", f"{sb:.6g}", f"{mb - ma:.6g}",
                         f"{lfc:.6g}", f"{p:.6g}"])

    class_tsv = os.path.join(args.outdir, "class_diff_normalized.tsv")
    with open(class_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["region_class", "metric", "n_regions",
                    f"mean_{label_a}", f"sd_{label_a}",
                    f"mean_{label_b}", f"sd_{label_b}",
                    "diff", "log2fc", "wilcoxon_p_paired_regions"])
        w.writerows(out_rows)
    print("wrote class_diff_normalized.tsv")
    for r in out_rows:
        print("  ", r)


if __name__ == "__main__":
    main()
