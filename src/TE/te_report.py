#!/usr/bin/env python
"""Build a cross-species TE statistics markdown report from TE GTFs.

Counts TE instances per hierarchy level (class -> family -> subfamily, from
GTF attributes class_id / family_id / gene_id) for each input GTF, then
renders one markdown document with per-species tables and a side-by-side
class-level comparison.

Usage:
    python te_report.py \
        -g /data/16T/Database/Reference/human/GENCODE/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf=Human_GRCh38 \
        -g /data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf=Mouse_GRCm39 \
        -g /data/16T/Database/Reference/mulatta/ENSEMBL/Mmul_10/rheMac10_rmsk_TE.gtf=Rhesus_rheMac10 \
        -o te_stats_report.md
"""

import argparse
import gzip
import re
from collections import defaultdict

ATTR = re.compile(r'(\S+) "([^"]*)"')
KEY_ORDER = ("class", "family", "subfamily")


def count_gtf(path):
    """Return {level: {element: [instances, bp]}}, class->families, family->subs."""
    levels = {k: defaultdict(lambda: [0, 0]) for k in KEY_ORDER}
    class_fams, fam_subs = defaultdict(set), defaultdict(set)
    sub2fam, fam2cls = {}, {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            a = dict(ATTR.findall(f[8]))
            cls, fam, sub = a.get("class_id", "NA"), a.get("family_id", "NA"), a.get("gene_id", "NA")
            bp = int(f[4]) - int(f[3]) + 1
            for key, name in zip(KEY_ORDER, (cls, fam, sub)):
                levels[key][name][0] += 1
                levels[key][name][1] += bp
            class_fams[cls].add(fam)
            fam_subs[fam].add(sub)
            fam2cls[fam], sub2fam[sub] = cls, fam
    return levels, class_fams, fam_subs, fam2cls, sub2fam


def fmt_int(n):
    return f"{n:,}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-g", "--gtf", action="append", required=True,
                    help="gtf_path=Label (repeatable)")
    ap.add_argument("-o", "--out", required=True, help="markdown output path")
    ap.add_argument("--top-families", type=int, default=15,
                    help="families shown per species table (default 15)")
    args = ap.parse_args()

    specs = []
    for item in args.gtf:
        path, _, label = item.rpartition("=")
        specs.append((label, path))

    data = {}
    for label, path in specs:
        levels, class_fams, fam_subs, fam2cls, sub2fam = count_gtf(path)
        data[label] = dict(levels=levels, class_fams=class_fams,
                            fam_subs=fam_subs, fam2cls=fam2cls, sub2fam=sub2fam)
        print(f"counted {label}: {len(levels['class'])} classes, "
              f"{len(levels['family'])} families, {len(levels['subfamily'])} subfamilies")

    lines = []
    lines.append("# 转座元件 (TE) 统计报告")
    lines.append("")
    lines.append(f"输入 GTF：{len(specs)} 个；层级链 `class_id -> family_id -> gene_id (subfamily)`。")
    lines.append("")

    # ---- overview ----
    lines.append("## 1. 总览")
    lines.append("")
    lines.append("| 物种/组装 | 记录总数 | Class 数 | Family 数 | Subfamily 数 | TE 总长 (Mb) |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    all_classes = []
    for label, _ in specs:
        lv = data[label]["levels"]
        total_n = lv["class"][max(lv["class"], key=lambda k: lv["class"][k][0])][0]
        # sum over classes is safer: recompute from subfamily? classes are disjoint; sum class rows
        total_n = sum(v[0] for v in lv["class"].values())
        total_bp = sum(v[1] for v in lv["class"].values())
        for c in lv["class"]:
            if c not in all_classes:
                all_classes.append(c)
        lines.append(f"| {label} | {fmt_int(total_n)} | {len(lv['class'])} | "
                     f"{len(lv['family'])} | {len(lv['subfamily'])} | {total_bp/1e6:.1f} |")
    lines.append("")

    # ---- per-species class tables ----
    lines.append("## 2. 各物种 Class 水平构成")
    lines.append("")
    for label, _ in specs:
        lv = data[label]["levels"]
        cf = data[label]["class_fams"]
        lines.append(f"### {label}")
        lines.append("")
        lines.append("| Class | 拷贝数 | 总长 (Mb) | Family 数 |")
        lines.append("|---|---:|---:|---:|")
        for name, (n, bp) in sorted(lv["class"].items(), key=lambda x: -x[1][0]):
            lines.append(f"| {name} | {fmt_int(n)} | {bp/1e6:.2f} | {len(cf[name])} |")
        lines.append("")

    # ---- cross-species class comparison ----
    lines.append("## 3. 跨物种 Class 拷贝数对比")
    lines.append("")
    header = ["Class"] + [lbl for lbl, _ in specs] + ["合计"]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + "---|" * len(header))
    for cls in sorted(all_classes,
                      key=lambda c: -sum(data[l]["levels"]["class"].get(c, [0, 0])[0]
                                         for l, _ in specs)):
        row = [cls if cls else "NA"]
        tot = 0
        for label, _ in specs:
            n = data[label]["levels"]["class"].get(cls, [0, 0])[0]
            tot += n
            row.append(fmt_int(n) if n else "0")
        row.append(fmt_int(tot))
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")

    # ---- per-species top families ----
    lines.append(f"## 4. 各物种主要 Family (Top {args.top_families})")
    lines.append("")
    for label, _ in specs:
        lv = data[label]["levels"]
        f2c = data[label]["fam2cls"]
        fs = data[label]["fam_subs"]
        lines.append(f"### {label}")
        lines.append("")
        lines.append("| Family | Class | 拷贝数 | 总长 (Mb) | Subfamily 数 |")
        lines.append("|---|---|---:|---:|---:|")
        for name, (n, bp) in sorted(lv["family"].items(), key=lambda x: -x[1][0])[:args.top_families]:
            lines.append(f"| {name} | {f2c[name]} | {fmt_int(n)} | {bp/1e6:.2f} | {len(fs[name])} |")
        lines.append("")

    # ---- per-species top subfamilies ----
    lines.append("## 5. 各物种主要 Subfamily (Top 15)")
    lines.append("")
    for label, _ in specs:
        lv = data[label]["levels"]
        s2f = data[label]["sub2fam"]
        f2c = data[label]["fam2cls"]
        lines.append(f"### {label}")
        lines.append("")
        lines.append("| Subfamily | Family | Class | 拷贝数 | 总长 (Mb) |")
        lines.append("|---|---|---|---:|---:|")
        for name, (n, bp) in sorted(lv["subfamily"].items(), key=lambda x: -x[1][0])[:15]:
            fam = s2f[name]
            lines.append(f"| {name} | {fam} | {f2c[fam]} | {fmt_int(n)} | {bp/1e6:.2f} |")
        lines.append("")

    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
