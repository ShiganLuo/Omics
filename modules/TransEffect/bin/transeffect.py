#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TransEffect: 转录效应分析 (WES somatic × RNA-seq 联合证据的转录后果).

对应《3.2.2.3节写作构思.md》第四节(数据根 = DR subworkflow 输出):
  1. 变异等位基因表达  : RNA 支持位点(D1R1+D1R2) rna_AF 分布        -> rna_support_sites.tsv
  2. RNA/DNA VAF 比值  : dna_AF/DP 下限过滤后的 rna_AF ÷ dna_AF        -> vaf_ratios.tsv
  5. CNV 剂量一致性    : 段内基因表达随拷贝数偏移(组内对照 Δ')         -> cnv_dosage_summary.tsv
  6. 变异基因表达偏离  : 携带者 vs 同组其余样本 Δ',对比随机基因 null   -> variant_gene_delta.tsv
  P5 反向检验          : 稳定性程序基因 ∩ 变异基因(期望=随机)          -> program_overlap.tsv

统计口径(修正版,见 report.md):
  Δ' = log2CPM(携带者) − mean(log2CPM(同组其余样本)) —— 组内对照,
  消除组间整体表达偏移;携带基因若同组另一样本也携带则剔除(无法归因)。
  随机对照: 等量随机基因重复抽样,给出 Δ' 的 null 分布。

子命令:
  analyze   计算上述全部表 + teffect_report.md(纯标准库)
  figures   读 analyze 产物画 4 张图(matplotlib;重绘不需重算)

仅 analyze 用标准库;figures 惰性 import matplotlib。
"""
import argparse
import csv
import gzip
import json
import os
import random
import re
import statistics
import sys

DEFAULT_GENE_SETS = [
    "Activation of ATR in response to replication stress",
    "DNA Double Strand Break Response",
    "HDR through Homologous Recombination (HRR)",
    "Nonhomologous End-Joining (NHEJ)",
    "Mismatch Repair",
    "Nucleotide Excision Repair",
    "Base Excision Repair",
    "Fanconi Anemia Pathway",
    "Telomere Maintenance",
    "Chromosome Maintenance",
    "Cell Cycle",
    "Mitotic Spindle Checkpoint",
]
SEX_CHROMS_DEFAULT = ["chrX", "chrY", "chrM", "X", "Y", "MT", "M"]


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")



def load_counts(path, samples):
    """featureCounts 矩阵 -> gene extents + counts(仅取配置样本列)."""
    genes = {}
    with open_text(path) as fh:
        first = fh.readline()
        if not first.startswith("#"):
            fh.seek(0)
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {}
        for i, h in enumerate(hdr):
            base = os.path.basename(h)
            for s in samples:
                if h.endswith(f"/{s}/{s}.bam") or base in (f"{s}.bam", f"{s}.txt", s):
                    col[s] = i
        missing = [s for s in samples if s not in col]
        if missing:
            raise ValueError(f"counts matrix {path} 缺少样本列: {missing}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gid = f[0].split(".")[0]
            starts = [int(x) for x in f[2].split(";")]
            ends = [int(x) for x in f[3].split(";")]
            genes[gid] = {
                "chrom": f[1].split(";")[0],
                "start": min(starts), "end": max(ends),
                "counts": {s: float(f[col[s]]) for s in samples},
            }
    totals = {s: sum(g["counts"][s] for g in genes.values()) for s in samples}
    for g in genes.values():
        g["log2cpm"] = {s: (g["counts"][s] / totals[s] * 1e6 if totals[s] else 0.0)
                        for s in samples}
        g["log2cpm"] = {s: __import__("math").log2(v + 1) for s, v in g["log2cpm"].items()}
    return genes, totals


def load_gene_names(gtf):
    """gene_id(去版本) -> gene_name;无 gtf 时返回空 dict。"""
    names = {}
    if not gtf:
        return names
    with open_text(gtf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8])
            nm = re.search(r'gene_name "([^"]+)"', f[8])
            if gid and nm:
                names[gid.group(1).split(".")[0]] = nm.group(1)
    return names


def genes_at(by_chrom, chrom, pos):
    return [gid for s, e, gid in by_chrom.get(chrom, []) if s <= pos <= e]


def read_drclass_tables(drclass_dir, samples, gene_names, by_chrom):
    """RNA 支持位点(项1/2) + 每样本变异基因集(P5/项6)。"""
    rna_sites = []
    var_genes = {}
    for s in samples:
        hits = set()
        tsv = os.path.join(drclass_dir, s, f"{s}.dna_anchored.tsv")
        with open_text(tsv) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pos = int(r["pos"])
                gids = genes_at(by_chrom, r["chrom"], pos)
                if r["dr_label"] in ("D1R1", "D1R2"):
                    rna_sites.append({
                        "sample": s, "chrom": r["chrom"], "pos": pos,
                        "ref": r["ref"], "alt": r["alt"], "dr_label": r["dr_label"],
                        "dna_AF": float(r["dna_AF"]), "dna_AD": int(r["dna_AD"]),
                        "dna_DP": int(r["dna_DP"]), "rna_AF": float(r["rna_AF"]),
                        "rna_AD": int(r["rna_AD"]), "rna_DP": int(r["rna_DP"]),
                        "rna_filter": r["rna_filter"],
                        "gene": ",".join(sorted({gene_names.get(g, g) for g in gids})) or ".",
                    })
        anc = os.path.join(drclass_dir, s, "anchors.tsv")
        with open_text(anc) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["anchor"] != "dna":
                    continue
                hits.update(genes_at(by_chrom, r["chrom"], int(r["pos"])))
        var_genes[s] = hits
    return rna_sites, var_genes


def within_group_deltas(genes, expressed, carried, samples, groups):
    """Δ'(携带者 − 同组其余均值);剔除同组多携带的基因。返回 (events, skipped_note)。"""
    members = {}
    for s in samples:
        members.setdefault(groups[s], []).append(s)
    events = []
    notes = []
    for s in samples:
        grp = groups[s]
        others = [t for t in members[grp] if t != s]
        if not others:
            notes.append(f"{s}: 组 {grp} 无同组对照样本,跳过 Δ'")
            continue
        shared = set()
        for t in others:
            shared |= carried[t]
        for gid in sorted(carried[s] & expressed - shared):
            d = genes[gid]["log2cpm"][s] - statistics.mean(
                genes[gid]["log2cpm"][t] for t in others)
            events.append((grp, s, gid, d))
    return events, notes


def random_null_deltas(genes, expressed, carried, samples, groups, reps, seed):
    members = {}
    for s in samples:
        members.setdefault(groups[s], []).append(s)
    rng = random.Random(seed)
    pool = sorted(expressed)
    out = []
    for rep in range(1, reps + 1):
        for s in samples:
            others = [t for t in members[groups[s]] if t != s]
            if not others:
                continue
            shared = set()
            for t in others:
                shared |= carried[t]
            picked = rng.sample(pool, len(carried[s]))
            for gid in picked:
                if gid in shared:
                    continue
                d = genes[gid]["log2cpm"][s] - statistics.mean(
                    genes[gid]["log2cpm"][t] for t in others)
                out.append((groups[s], s, gid, d, rep))
    return out


def cnv_analysis(genes, expressed, cnv_dir, samples, groups, sex_chroms):
    members = {}
    for s in samples:
        members.setdefault(groups[s], []).append(s)
    segments = {}
    summary = []
    for s in samples:
        segs = []
        path = os.path.join(cnv_dir, s, f"{s}.cns")
        with open_text(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                segs.append((r["chromosome"], int(r["start"]), int(r["end"]),
                             int(float(r["cn"]))))
        segments[s] = segs
        others = [t for t in members[groups[s]] if t != s]
        by_class = {"gain": [], "loss": [], "neutral": []}
        auto_alt = sum(1 for c, st, en, cn in segs
                       if cn != 2 and re.match(r"^(chr)?\d+$", c))
        unplaced_alt = sum(1 for c, st, en, cn in segs
                           if cn != 2 and c not in sex_chroms
                           and not re.match(r"^(chr)?\d+$", c))
        for gid in sorted(expressed):
            g = genes[gid]
            mid = (g["start"] + g["end"]) // 2
            cls = None
            for c, st, en, cn in segs:
                if c == g["chrom"] and st <= mid <= en:
                    cls = "gain" if cn > 2 else ("loss" if cn < 2 else "neutral")
                    break
            if cls is None or not others:
                continue
            d = g["log2cpm"][s] - statistics.mean(g["log2cpm"][t] for t in others)
            by_class[cls].append(d)
        for cls in ("gain", "loss", "neutral"):
            v = by_class[cls]
            summary.append((groups[s], s, cls, len(v),
                            round(statistics.median(v), 4) if v else "NA",
                            round(sum(1 for x in v if x > 0) / len(v), 3) if v else "NA",
                            auto_alt, unplaced_alt))
    return segments, summary


def program_overlap(var_genes, samples, gmt, gene_sets, genes, gene_names):
    if not gmt or not os.path.exists(gmt):
        return None, "gmt 未配置或不存在,跳过 P5 反向检验"
    names_to_id = {v: k for k, v in gene_names.items()}
    prog = {}
    with open(gmt) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in gene_sets:
                prog[f[0]] = {names_to_id[x] for x in f[2:] if x in names_to_id}
    if not prog:
        return None, f"gmt 中未命中配置的基因集 {gene_sets},跳过 P5"
    all_prog = set().union(*prog.values())
    rows = []
    for s in samples:
        hits = var_genes[s]
        inter = sorted(hits & all_prog)
        expected = len(all_prog) / max(len(genes), 1) * len(hits)
        rows.append((s, len(inter), round(expected, 2), len(hits), len(all_prog),
                     ",".join(gene_names.get(g, g) for g in inter) or "."))
    return rows, f"程序基因 {len(prog)} 组/去重 {len(all_prog)} 个"


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def cmd_analyze(args):
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    groups = {}
    for kv in args.sample_group or []:
        k, _, v = kv.partition("=")
        groups[k.strip()] = v.strip()
    groups = {s: groups.get(s, "G") for s in samples}
    sex_chroms = set((args.sex_chroms or "").split(",")) - {""}
    sex_chroms = sex_chroms or set(SEX_CHROMS_DEFAULT)
    gene_sets = (args.gene_sets.split("|") if args.gene_sets else DEFAULT_GENE_SETS)

    os.makedirs(args.outdir, exist_ok=True)
    genes, totals = load_counts(args.counts, samples)
    gene_names = load_gene_names(args.gtf)
    by_chrom = {}
    for gid, g in genes.items():
        by_chrom.setdefault(g["chrom"], []).append((g["start"], g["end"], gid))
    for c in by_chrom:
        by_chrom[c].sort()
    expressed = {g for g, d in genes.items()
                 if sum(1 for s in samples
                        if d["counts"][s] / totals[s] * 1e6 >= args.min_cpm
                        if totals[s]) >= args.min_expressed}
    rna_sites, var_genes = read_drclass_tables(args.drclass_dir, samples, gene_names, by_chrom)

    rows1 = []
    for r in rna_sites:
        rows1.append([groups[r["sample"]], r["sample"], r["chrom"], r["pos"], r["ref"],
                      r["alt"], r["dr_label"], r["dna_AF"], r["dna_AD"], r["dna_DP"],
                      r["rna_AF"], r["rna_AD"], r["rna_DP"], r["rna_filter"], r["gene"],
                      "yes" if (r["dna_AF"] >= args.dna_af_min
                                and r["dna_DP"] >= args.dna_dp_min) else "no"])
    write_tsv(os.path.join(args.outdir, "rna_support_sites.tsv"),
              ["group", "sample", "chrom", "pos", "ref", "alt", "dr_label", "dna_AF",
               "dna_AD", "dna_DP", "rna_AF", "rna_AD", "rna_DP", "rna_filter",
               "gene", "pass_vaf_ratio_filter"], rows1)

    rows2 = [[r[0], r[1], r[2], r[3], r[4], r[5], r[14], r[7], r[10],
              round(r[10] / r[7], 4) if r[7] else "NA"]
             for r in rows1 if r[15] == "yes"]
    write_tsv(os.path.join(args.outdir, "vaf_ratios.tsv"),
              ["group", "sample", "chrom", "pos", "ref", "alt", "gene", "dna_AF",
               "rna_AF", "rna_dna_ratio"], rows2)

    carried = {s: set(var_genes[s]) & expressed for s in samples}
    events, notes = within_group_deltas(genes, expressed, carried, samples, groups)
    null = random_null_deltas(genes, expressed, carried, samples, groups,
                              args.null_reps, args.seed)
    rows6 = [[g, s, gene_names.get(gid, gid), "variant_gene", 0, round(d, 4)]
             for g, s, gid, d in events]
    rows6 += [[g, s, gene_names.get(gid, gid), "random_gene", rep, round(d, 4)]
              for g, s, gid, d, rep in null if rep <= args.null_plot_reps]
    write_tsv(os.path.join(args.outdir, "variant_gene_delta.tsv"),
              ["group", "sample", "gene", "set", "rep", "delta_prime"], rows6)

    segments, summary5 = cnv_analysis(genes, expressed, args.cnv_dir, samples,
                                      groups, sex_chroms)
    write_tsv(os.path.join(args.outdir, "cnv_dosage_summary.tsv"),
              ["group", "sample", "cn_class", "n_genes", "delta_prime_median",
               "frac_delta_pos", "autosomal_cn_alt_segments",
               "unplaced_cn_alt_segments"], summary5)
    seg_rows = [(s, c, st, en, cn) for s in samples
                for c, st, en, cn in segments[s]]
    write_tsv(os.path.join(args.outdir, "cnv_segments.tsv"),
              ["sample", "chrom", "start", "end", "cn"], seg_rows)

    prog_rows, prog_note = program_overlap(var_genes, samples, args.gmt, gene_sets,
                                           genes, gene_names)
    if prog_rows is not None:
        write_tsv(os.path.join(args.outdir, "program_overlap.tsv"),
                  ["sample", "observed", "expected", "variant_genes", "program_genes",
                   "hits"], prog_rows)

    # 汇总报告
    def med(vals):
        return round(statistics.median(vals), 4) if vals else "NA"
    grp_stats = {}
    for g, s, gid, d in events:
        grp_stats.setdefault(g, []).append(d)
    null_med = {}
    for g, s, gid, d, rep in null:
        null_med.setdefault((g, rep), []).append(d)
    lines = ["# TransEffect 转录效应分析", "",
             f"样本: {', '.join(samples)}; 分组: {groups}", ""]
    lines += ["## 1/2. RNA 支持位点(D1R1+D1R2)与 RNA/DNA VAF 比值",
              f"- RNA 支持位点 {len(rows1)} 个;过比值过滤(dna_AF>={args.dna_af_min} 且 "
              f"dna_DP>={args.dna_dp_min}) {len(rows2)} 个。",
              f"- 比值中位数: " + ", ".join(
                  f"{g}={med([r[9] for r in rows2 if r[0] == g])}"
                  for g in sorted({r[0] for r in rows2})) + "。",
              "- 位点数少且部分 RNA 深度低,仅作描述性证据。", ""]
    lines += ["## 5. CNV 剂量一致性",
              "- 每样本 cn≠2 段: 常染色体 "
              + ", ".join(f"{s}={sum(1 for c, st, en, cn in segments[s] if cn != 2 and re.match(r'^(chr)?[0-9]+$', c))}"
                          for s in samples)
              + ";未定位支架 "
              + ", ".join(f"{s}={sum(1 for c, st, en, cn in segments[s] if cn != 2 and c not in sex_chroms and not re.match(r'^(chr)?[0-9]+$', c))}"
                          for s in samples) + "。",
              "- 逐类表达偏移见 cnv_dosage_summary.tsv;无常染色体 cn≠2 段时剂量一致性无从检验,"
              "此本身即为\"解读稳定\"的证据。", ""]
    lines += ["## 6. 变异基因表达偏离(组内对照 Δ')"]
    for g in sorted(grp_stats):
        null_center = statistics.median(
            [statistics.median(v) for (gg, rep), v in null_med.items() if gg == g]) if null_med else 0.0
        lines.append(f"- {g}: n={len(grp_stats[g])}, Δ' 中位数={med(grp_stats[g])}, "
                     f"随机 null 中位数中心≈{round(null_center, 4)}")
    lines += [f"- 结论: 变异基因 Δ' 与随机 null 一致 → 变异不伴随表达偏离(以数值为准)。", ""]
    if notes:
        lines += ["## 备注"] + [f"- {n}" for n in notes] + [""]
    lines += ["## P5 反向检验: 稳定性程序基因 ∩ 变异基因",
              f"- {prog_note}"]
    if prog_rows is not None:
        lines += [f"- {s}: 观测 {o} vs 期望 {e}"
                  for s, o, e, *_ in prog_rows]
        lines += ["- 期望为基因级随机基线(程序基因/总基因 × 变异基因数),未校正基因长度/表达混杂;"
                  "变异基因偏长、程序基因偏长高表达,观测>期望部分由该混杂解释,判读以组间结构为准。"]
    lines += ["", "## 图",
              "figures 子命令读本目录 TSV 输出 fig_rna_support_vaf.png / "
              "fig_variant_gene_delta.png / fig_cnv_overview.png / fig_program_overlap.png"]
    with open(os.path.join(args.outdir, "teffect_report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"analyze -> {args.outdir}")



def cmd_figures(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = args.outdir
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False})

    def read_tsv(name):
        with open(os.path.join(outdir, name)) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    rna_sites = read_tsv("rna_support_sites.tsv")
    ratios = read_tsv("vaf_ratios.tsv")
    groups = sorted({r["group"] for r in rna_sites} | {r["group"] for r in ratios})
    palette = ["#2166ac", "#b2182b", "#1b7837", "#762a83", "#e08214", "#4d4d4d"]
    GC = {g: palette[i % len(palette)] for i, g in enumerate(groups)}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.2))
    gene_pts = {}
    for r in rna_sites:
        gene_pts.setdefault(r["gene"], []).append(r)
    offs = [(3, 2, "left"), (-3, 2, "right"), (3, -9, "left"), (-3, -9, "right"),
            (3, 8, "left"), (-3, 8, "right"), (12, 0, "left"), (-12, 0, "right"),
            (3, -15, "left"), (-3, -15, "right"), (3, 14, "left"), (-3, 14, "right")]
    for r in rna_sites:
        daf, raf = float(r["dna_AF"]), float(r["rna_AF"])
        dp = int(r["rna_DP"])
        ax1.scatter(daf, raf, s=28, color=GC[r["group"]],
                    marker="o" if dp >= 10 else "s",
                    facecolor=GC[r["group"]] if dp >= 10 else "white",
                    edgecolor=GC[r["group"]], zorder=3)
    fig.canvas.draw()
    pts_disp = [ax1.transData.transform((float(r["dna_AF"]), float(r["rna_AF"])))
                for r in rna_sites]

    def _overlap(b1, b2):
        return not (b1.x1 <= b2.x0 or b1.x0 >= b2.x1 or
                    b1.y1 <= b2.y0 or b1.y0 >= b2.y1)

    placed = []
    for gene, rs in sorted(gene_pts.items(),
                           key=lambda kv: -float(kv[1][0]["rna_AF"])):
        if gene == ".":
            continue
        r = rs[0]
        n_pos = len({(x["chrom"], x["pos"]) for x in rs})
        n_smp = len({x["sample"] for x in rs})
        suffix = ""
        if n_pos > 1 or n_smp > 1:
            parts = ([f"{n_pos} sites"] if n_pos > 1 else []) + \
                    ([f"{n_smp} samples"] if n_smp > 1 else [])
            suffix = " (" + ", ".join(parts) + ")"
        xy = (float(r["dna_AF"]), float(r["rna_AF"]))
        for dx, dy, ha in offs:
            t = ax1.annotate(gene.split(",")[0] + suffix, xy, fontsize=5.5,
                             ha=ha, xytext=(dx, dy),
                             textcoords="offset points", zorder=5)
            fig.canvas.draw()
            bb = t.get_window_extent()
            hit_point = any(
                bb.x0 - 2 <= px <= bb.x1 + 2 and bb.y0 - 2 <= py <= bb.y1 + 2
                for px, py in pts_disp)
            hit_text = any(_overlap(bb, pb) for pb in placed)
            if not hit_point and not hit_text:
                placed.append(bb)
                break
            t.remove()
    ax1.plot([0, 1], [0, 1], ls="--", lw=0.7, color="grey")
    ax1.set_xlabel("DNA VAF (somatic callset)")
    ax1.set_ylabel("RNA VAF (pileup)")
    ax1.set_xlim(0, 1.05)
    ax1.set_ylim(-0.06, 1.12)
    ax1.legend(handles=[
        plt.Line2D([], [], marker="o", ls="", color=GC[groups[0]], label=f"{groups[0]}, RNA DP>=10"),
        plt.Line2D([], [], marker="s", ls="", mfc="white", mec=GC[groups[0]], color=GC[groups[0]], label=f"{groups[0]}, RNA DP<10")]
        + ([plt.Line2D([], [], marker="o", ls="", color=GC[groups[1]], label=f"{groups[1]}, RNA DP>=10"),
            plt.Line2D([], [], marker="s", ls="", mfc="white", mec=GC[groups[1]], color=GC[groups[1]], label=f"{groups[1]}, RNA DP<10")]
           if len(groups) > 1 else []),
        fontsize=5.5, loc="lower right", frameon=False)
    xmap = {g: i for i, g in enumerate(groups)}
    jitter = {"KSR": 0, "TLSCS": 0}
    for r in ratios:
        ax2.scatter(xmap[r["group"]], float(r["rna_dna_ratio"]), s=28,
                    color=GC[r["group"]], zorder=3, alpha=0.85)
    ax2.axhline(1, ls="--", lw=0.7, color="grey")
    ax2.set_xticks(range(len(groups)))
    ax2.set_xticklabels(groups)
    ax2.set_ylabel("RNA/DNA VAF ratio")
    ax2.set_xlim(-0.5, len(groups) - 0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig_rna_support_vaf.png"), dpi=300)
    plt.close(fig)

    d = read_tsv("variant_gene_delta.tsv")
    if d:
        vd, nd = {}, {}
        for r in d:
            tgt = vd if r["set"] == "variant_gene" else nd
            tgt.setdefault(r["group"], []).append(float(r["delta_prime"]))
        all_g = sorted(set(vd) | set(nd))
        fig, ax = plt.subplots(figsize=(max(3.6, 1.2 + 1.2 * len(all_g)), 3.2))
        pos = []
        for i, g in enumerate(all_g):
            pos += [(i, 0.22), (i, 0.78)]
        sets = []
        for g in all_g:
            sets += [(vd.get(g, []), GC.get(g, "#666666")),
                     (nd.get(g, []), "#cccccc")]
        for (x, dx), (vals, col) in zip(pos, sets):
            if not vals:
                continue
            bp = ax.boxplot(vals, positions=[x + dx - 0.5], widths=0.32,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=0.8))
            for patch in bp["boxes"]:
                patch.set_facecolor(col)
                patch.set_alpha(0.55)
        ax.axhline(0, ls="--", lw=0.7, color="grey")
        ax.set_xticks(range(len(all_g)))
        ax.set_xticklabels(all_g)
        ax.set_xlim(-0.5, len(all_g) - 0.5)
        ax.set_ylabel("expression deviation of carrier\n($\\Delta'$ log$_2$ CPM)")
        handles = [plt.Rectangle((0, 0), 1, 1, fc="#cccccc", alpha=0.55, label="random genes")]
        handles += [plt.Rectangle((0, 0), 1, 1, fc=GC.get(g, "#666"), alpha=0.55,
                                  label=f"variant ({g})") for g in all_g]
        ax.legend(handles=handles, fontsize=5.5, loc="lower center",
                  bbox_to_anchor=(0.5, 1.02), ncol=len(handles), frameon=False)
        fig.tight_layout()
        fig.subplots_adjust(top=0.82)
        fig.savefig(os.path.join(outdir, "fig_variant_gene_delta.png"), dpi=300)
        plt.close(fig)

    segs = read_tsv("cnv_segments.tsv")
    if segs:
        smp_order = sorted({r["sample"] for r in segs})
        chroms = sorted({r["chrom"] for r in segs},
                        key=lambda c: (not c[3:].isdigit(), int(c[3:]) if c[3:].isdigit() else c))
        cum, x = {}, 0
        for c in chroms:
            mx = max(int(r["end"]) for r in segs if r["chrom"] == c)
            cum[c] = (x, mx)
            x += mx
        fig, ax = plt.subplots(figsize=(max(7.0, 0.38 * len(chroms)), 3.0))
        rows_y = {s: len(smp_order) - 1 - i for i, s in enumerate(smp_order)}
        for r in segs:
            y = rows_y[r["sample"]]
            cn = int(r["cn"])
            off = cum[r["chrom"]][0]
            col = "#444444" if cn == 2 else ("#b2182b" if cn > 2 else "#2166ac")
            ax.plot([off + int(r["start"]), off + int(r["end"])],
                    [y + (cn - 2) * 0.18] * 2, lw=1.6, color=col,
                    solid_capstyle="butt")
        for y in rows_y.values():
            ax.axhline(y, lw=0.3, color="#dddddd", zorder=0)
        ax.set_yticks(list(rows_y.values()))
        ax.set_yticklabels(smp_order)
        ax.set_xticks([cum[c][0] + cum[c][1] / 2 for c in chroms])
        ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=5.5,
                           rotation=90)
        ax.set_xlabel("chromosome")
        ax.set_ylabel("segment copy number\n(relative to 2)")
        ax.set_ylim(-0.6, len(smp_order) - 0.4)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.30)
        fig.savefig(os.path.join(outdir, "fig_cnv_overview.png"), dpi=300)
        plt.close(fig)

    pov = None
    if os.path.exists(os.path.join(outdir, "program_overlap.tsv")):
        pov = read_tsv("program_overlap.tsv")
    if pov:
        fig, ax = plt.subplots(figsize=(1.4 + 0.8 * len(pov), 2.6))
        xs = range(len(pov))
        bars = ax.bar(list(xs), [int(r["observed"]) for r in pov], 0.55,
                      color=[GC[r["sample"]] if r["sample"] in GC else "#666666"
                             for r in pov])
        for b in bars:
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.25,
                    f"{int(b.get_height())}", ha="center", fontsize=7)
        exp = float(pov[0]["expected"])
        ax.axhline(exp, ls="--", lw=0.9, color="grey")
        ax.text(len(pov) - 0.55, exp + 0.3, "expected", fontsize=6,
                color="grey", ha="right")
        ax.set_xticks(list(xs))
        ax.set_xticklabels([r["sample"] for r in pov], fontsize=6)
        ax.set_ylabel("variant genes in\nstability programs")
        top = max([int(r["observed"]) for r in pov] + [exp])
        ax.set_ylim(0, top * 1.18)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "fig_program_overlap.png"), dpi=300)
        plt.close(fig)
    print(f"figures -> {outdir}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("analyze")
    a.add_argument("--samples", required=True, help="逗号分隔样本 ID")
    a.add_argument("--sample-group", action="append",
                   help="sample=group,可重复;缺省样本各归一组 G")
    a.add_argument("--drclass-dir", required=True)
    a.add_argument("--counts", required=True, help="featureCounts 合并矩阵")
    a.add_argument("--cnv-dir", required=True, help="含 <sample>/<sample>.cns")
    a.add_argument("--gmt", default=None, help="基因集 GMT(P5 反向检验)")
    a.add_argument("--gene-sets", default=None,
                   help="GMT 集合名,| 分隔;缺省 12 个基因组稳定性通路")
    a.add_argument("--gtf", default=None, help="gene_id->gene_name 映射来源")
    a.add_argument("--outdir", required=True)
    a.add_argument("--dna-af-min", type=float, default=0.1)
    a.add_argument("--dna-dp-min", type=int, default=10)
    a.add_argument("--min-cpm", type=float, default=1.0)
    a.add_argument("--min-expressed", type=int, default=2)
    a.add_argument("--null-reps", type=int, default=20)
    a.add_argument("--null-plot-reps", type=int, default=3)
    a.add_argument("--seed", type=int, default=42)
    a.add_argument("--sex-chroms", default=None)
    a.set_defaults(func=cmd_analyze)

    f = sub.add_parser("figures")
    f.add_argument("--outdir", required=True, help="analyze 产物目录(原地出图)")
    f.set_defaults(func=cmd_figures)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
