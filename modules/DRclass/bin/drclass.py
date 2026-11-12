#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DRclass: DNA(WES somatic) x RNA-seq 联合证据分类.

设计文档: /data/16T/luosg/output/Totipotency/Mutation/3.2.2.3节写作构思.md 第三节
状态机: REJECTED > HIT-pass > HIT-strong > HIT-weak > ABSENT > NOINFO
标签空间: 表1(DNA锚) D1R1/D1R2/D1R3/D1R4/D1R6/D1R5
          表2(RNA锚) D2R1/D3R1/D6R1/D5R1/D0R1 (+R1-ANCHOR-REJ 锚级伪影剔除)

子命令:
  anchors   从 somatic VCF(PASS) + RNA VCF 生成锚表与锚 BED
  classify  逐位点判定 DR 标签, 输出 O1(dna_anchored) O2(rna_anchored)
  summary   汇总 O1/O2 -> DRclass_summary.tsv (计数/Tier/穿透率/质控)

仅用 python3 标准库; VCF 读取兼容 bgzip (gzip 模块).
"""
import argparse
import csv
import gzip

# ---------- 通用 ----------

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def norm_key(chrom, pos, ref, alt):
    return (str(chrom), int(pos), str(ref).upper(), str(alt).upper().split(",")[0])


def parse_vcf_record(line):
    """返回 dict 或 None(表头). 样本列全部保留."""
    if line.startswith("#"):
        return None
    f = line.rstrip("\n").split("\t")
    fmt_keys = f[8].split(":") if len(f) > 9 else []
    samples = {}
    for i, s in enumerate(f[9:]):
        vals = s.split(":")
        samples["s%d" % i] = dict(zip(fmt_keys, vals)) if fmt_keys else {}
    alts = f[4].split(",") if f[4] not in (".", "") else []
    return {
        "chrom": f[0], "pos": int(f[1]), "ref": f[3], "alts": alts,
        "filter": f[6], "samples": samples,
    }


def iter_vcf(path):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            rec = parse_vcf_record(line)
            if rec:
                yield rec


def sample_ad_dp_af(sm):
    """从 FORMAT dict 取 AD/DP, 计算 AF (无 AF 字段时 AD 派生). 全部返回数值(缺省0)."""
    ad = sm.get("AD", ".")
    dp = sm.get("DP", ".")
    ref_ad, alt_ad = 0, 0
    if ad not in (".", "", None):
        parts = [int(x) if x.isdigit() else 0 for x in ad.split(",")]
        ref_ad = parts[0]
        alt_ad = sum(parts[1:])
    dp_i = int(dp) if (dp or "").isdigit() else 0
    af = sm.get("AF", ".")
    if af not in (".", "", None):
        try:
            af_f = float(af.split(",")[0])
        except ValueError:
            af_f = 0.0
    elif dp_i > 0:
        af_f = alt_ad / dp_i
    else:
        af_f = 0.0
    return ref_ad, alt_ad, dp_i, af_f


# ---------- anchors ----------

def cmd_anchors(args):
    dna_keys = set()
    for rec in iter_vcf(args.somatic_vcf):
        if rec["filter"] != "PASS":
            continue
        for alt in rec["alts"]:
            dna_keys.add(norm_key(rec["chrom"], rec["pos"], rec["ref"], alt))
    rna_seen, rna_list = set(), []
    for rec in iter_vcf(args.rna_vcf):
        for alt in rec["alts"]:
            k = norm_key(rec["chrom"], rec["pos"], rec["ref"], alt)
            if k not in dna_keys and k not in rna_seen:
                rna_seen.add(k)
                rna_list.append(k)
    rna_keys = sorted(rna_list, key=lambda x: (x[0], x[1]))
    with open(args.out_tsv, "w", newline="") as tsv, open(args.out_bed, "w") as bed:
        tsv.write("chrom\tpos\tref\talt\tanchor\n")
        for (c, p, r, a) in sorted(dna_keys, key=lambda x: (x[0], x[1])):
            tsv.write(f"{c}\t{p}\t{r}\t{a}\tdna\n")
            bed.write(f"{c}\t{p - 1}\t{p - 1 + len(r)}\n")
        for (c, p, r, a) in rna_keys:
            tsv.write(f"{c}\t{p}\t{r}\t{a}\trna\n")
            bed.write(f"{c}\t{p - 1}\t{p - 1 + len(r)}\n")
    print(f"anchors: dna={len(dna_keys)} rna={len(rna_keys)} -> {args.out_tsv}")


# ---------- classify 辅助装载 ----------

def load_vcf_by_pos(path):
    bypos = {}
    for rec in iter_vcf(path):
        bypos.setdefault((rec["chrom"], rec["pos"]), []).append(rec)
    return bypos


def load_keys_set(path):
    keys = set()
    for rec in iter_vcf(path):
        for alt in rec["alts"]:
            keys.add(norm_key(rec["chrom"], rec["pos"], rec["ref"], alt))
    return keys


def load_germline_hits(path, af_min, dp_min):
    """返回 {key: 'germline'} — 自身 WES germline 中 AF>=af_min 且 DP>=dp_min 的位点."""
    hits = {}
    if not path or path in (".", "None"):
        return hits
    for rec in iter_vcf(path):
        for alt in rec["alts"]:
            sm = rec["samples"].get("s0", {})
            _, _, dp_i, af_f = sample_ad_dp_af(sm)
            if af_f >= af_min and dp_i >= dp_min:
                hits[norm_key(rec["chrom"], rec["pos"], rec["ref"], alt)] = "germline"
    return hits


def load_pileup(path):
    """bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD\\t%DP]\\n' 的输出."""
    pile = {}
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            chrom, pos, ref, alt, ad, dp = f[0], int(f[1]), f[2], f[3], f[4], f[5]
            pile.setdefault((chrom, pos), []).append(
                {"ref": ref, "alts": alt.split(","), "ad": ad, "dp": dp})
    return pile


def pileup_counts(pile, chrom, pos, ref, alt):
    """返回 (found, ref_ad, alt_ad, dp). 未覆盖 -> (False, 0, 0, 0)."""
    recs = pile.get((chrom, pos), [])
    for r in recs:
        if r["ref"] != ref:
            continue
        dp_i = int(r["dp"]) if r["dp"].isdigit() else 0
        parts = [int(x) if x.isdigit() else 0 for x in r["ad"].split(",")]
        ref_ad = parts[0] if parts else 0
        alt_ad = 0
        for i, a in enumerate(r["alts"]):
            if a == alt and i + 1 < len(parts):
                alt_ad = parts[i + 1]
        return True, ref_ad, alt_ad, dp_i
    return False, 0, 0, 0


def load_intervals(paths):
    """BED(0-based 区间) 或 VCF(单点) 混合收集为区间列表."""
    ivals = []
    for p in paths or []:
        with open_text(p) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if not f[0]:
                    continue
                try:
                    if f[0].startswith("chr") or not f[0][0].isdigit():
                        ivals.append((f[0], int(f[1]), int(f[2])))      # BED
                    else:
                        pos = int(f[1])
                        ivals.append((f[0], pos - 1, pos))               # VCF
                except (ValueError, IndexError):
                    continue
    return ivals


def in_window(ivals, chrom, pos, w):
    return any(c == chrom and s - w <= pos <= e + w for (c, s, e) in ivals)


# ---------- classify ----------

OUT_HEADER = [
    "sample", "chrom", "pos", "ref", "alt", "anchor",
    "dna_call_status", "dna_filter", "dna_AF", "dna_AD", "dna_DP",
    "rna_call_status", "rna_filter", "rna_AD", "rna_DP", "rna_AF",
    "artifact_flags", "germline_hit_source", "dr_label", "tier",
]

TIER = {
    "D1R1": "A", "D1R2": "A", "D1R3": "sensitivity", "D1R4": "rejected",
    "D1R6": "B", "D1R5": "U",
    "D2R1": "C1", "D3R1": "C1", "D6R1": "C2", "D5R1": "C3",
    "D0R1": "germline", "R1-ANCHOR-REJ": "rejected",
}


def cmd_classify(args):
    min_af, min_ad, min_dp = args.min_af, args.min_ad, args.min_dp

    somatic = load_vcf_by_pos(args.somatic_vcf)
    rna = load_vcf_by_pos(args.rna_vcf)
    pile = load_pileup(args.pileup)
    germline = load_germline_hits(args.germline_vcf, args.germline_af, args.germline_dp)
    strain = {}
    for sv in args.strain_vcf or []:
        strain.update({k: "strain" for k in load_keys_set(sv)})
    splice = load_intervals(args.splice_bed)
    editing = load_intervals(args.editing_bed)

    anchors = []
    with open_text(args.anchors) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            anchors.append((row["chrom"], int(row["pos"]), row["ref"], row["alt"], row["anchor"]))

    out1 = open(args.out_dna, "w", newline="")
    out2 = open(args.out_rna, "w", newline="")
    w1 = csv.writer(out1, delimiter="\t")
    w2 = csv.writer(out2, delimiter="\t")
    w1.writerow(OUT_HEADER)
    w2.writerow(OUT_HEADER)

    def find_rec(bypos, chrom, pos, ref, alt):
        for rec in bypos.get((chrom, pos), []):
            if rec["ref"] == ref and alt in rec["alts"]:
                return rec
        return None

    def rna_artifacts(key, rrec):
        flags = []
        if key in germline:
            flags.append("germline")
        if key in strain:
            flags.append("strain")
        if in_window(splice, key[0], key[1], args.splice_window):
            flags.append("splice")
        if in_window(editing, key[0], key[1], 0):
            flags.append("editing")
        if rrec is not None and len(rrec["alts"]) > 1:
            flags.append("multiallelic")
        return flags

    def dna_vals(srec):
        """somatic VCF case 列(s1; 无则 s0) 的 AD/DP/AF."""
        if srec is None:
            return 0.0, 0, 0
        sm = srec["samples"].get("s1", srec["samples"].get("s0", {}))
        _, alt_ad, dp_i, af_f = sample_ad_dp_af(sm)
        return af_f, alt_ad, dp_i

    def rna_vals_from_call(rrec):
        """RNA VCF 的 AD/DP/AF(仅锚位定性用; VAF 测量以 pileup 为准)."""
        sm = rrec["samples"].get("s0", {})
        _, alt_ad, dp_i, af_f = sample_ad_dp_af(sm)
        return alt_ad, dp_i, af_f

    for (chrom, pos, ref, alt, anchor) in anchors:
        key = norm_key(chrom, pos, ref, alt)
        found, p_ref, p_alt, p_dp = pileup_counts(pile, chrom, pos, ref, alt)
        srec = find_rec(somatic, chrom, pos, ref, alt)
        rrec = find_rec(rna, chrom, pos, ref, alt)
        gsrc = germline.get(key, "") or strain.get(key, "")

        if anchor == "dna":
            # ---- 表1: D 轴恒 HIT-pass, 判 R 轴 ----
            d_af, d_alt, d_dp = dna_vals(srec)
            rna_status, r_filter = "not_detected", "."
            r_ad, r_dp, r_af = p_alt, p_dp, (p_alt / p_dp if p_dp else 0.0)
            art = []
            if rrec is not None:
                rna_status, r_filter = "detected", rrec["filter"]
                c_ad, c_dp, c_af = rna_vals_from_call(rrec)
                if found:                     # pileup 数值优先(M3 是 VAF 测量来源)
                    r_ad, r_dp, r_af = p_alt, p_dp, (p_alt / p_dp if p_dp else 0.0)
                else:
                    r_ad, r_dp, r_af = c_ad, c_dp, c_af
                art = rna_artifacts(key, rrec)
                label = "D1R4" if art else "D1R1"
            else:
                if p_dp > min_dp and r_af >= min_af and p_alt > min_ad:
                    label = "D1R2"
                elif p_dp > min_dp and r_af >= min_af:
                    label = "D1R3"
                elif p_dp > min_dp:
                    label = "D1R6"
                else:
                    label = "D1R5"
                r_ad, r_dp = p_alt, p_dp
            w1.writerow([args.sample, chrom, pos, ref, alt, "dna",
                         srec["filter"] if srec else "PASS",
                         srec["filter"] if srec else "PASS",
                         round(d_af, 4), d_alt, d_dp,
                         rna_status, r_filter, r_ad, r_dp, round(r_af, 4),
                         ";".join(art) if art else ".", gsrc,
                         label, TIER[label]])
        else:
            # ---- 表2: R 轴恒 HIT-pass, 判 D 轴; 锚级伪影先剔除 ----
            art = rna_artifacts(key, rrec)
            c_ad, c_dp, c_af = rna_vals_from_call(rrec or {"samples": {"s0": {}}})
            if found:
                r_ad, r_dp, r_af = p_alt, p_dp, (p_alt / p_dp if p_dp else 0.0)
            else:
                r_ad, r_dp, r_af = c_ad, c_dp, c_af
            d_af, d_alt, d_dp = dna_vals(srec)
            if "germline" in art or "strain" in art:
                label = "D0R1"
            elif any(f in art for f in ("splice", "editing", "multiallelic")):
                label = "R1-ANCHOR-REJ"
            elif srec is not None:
                if d_af >= min_af and d_dp > min_dp and d_alt > min_ad:
                    label = "D2R1"
                elif d_af >= min_af and d_dp > min_dp:
                    label = "D3R1"
                elif d_dp > min_dp:
                    label = "D6R1"
                else:
                    label = "D5R1"
            else:
                label = "D5R1"
            w2.writerow([args.sample, chrom, pos, ref, alt, "rna",
                         srec["filter"] if srec is not None else "not_in_callset",
                         srec["filter"] if srec is not None else ".",
                         round(d_af, 4), d_alt, d_dp,
                         "detected", rrec["filter"] if rrec else ".",
                         r_ad, r_dp, round(r_af, 4),
                         ";".join(art) if art else ".", gsrc,
                         label, TIER[label]])
    out1.close()
    out2.close()
    print(f"classify {args.sample}: tables -> {args.out_dna}, {args.out_rna}")


# ---------- summary ----------

def cmd_summary(args):
    labels = ["D1R1", "D1R2", "D1R3", "D1R4", "D1R6", "D1R5",
              "D2R1", "D3R1", "D6R1", "D5R1", "D0R1", "R1-ANCHOR-REJ"]
    rows = []
    for t in args.tables:
        counts = {l: 0 for l in labels}
        sample = ""
        with open_text(t) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sample = r["sample"]
                if r["dr_label"] in counts:
                    counts[r["dr_label"]] += 1
        pen_num = counts["D1R1"] + counts["D1R2"]
        pen_den = pen_num + counts["D1R6"]
        penetrance = (pen_num / pen_den) if pen_den else None
        rna_det = counts["D1R1"] + counts["D1R4"]
        r4_frac = (counts["D1R4"] / rna_det) if rna_det else None
        qc = "WARN:r4_fraction>%.2f" % args.max_r4 if (r4_frac is not None and r4_frac > args.max_r4) else "OK"
        rows.append([sample] + [counts[l] for l in labels] + [
            round(penetrance, 4) if penetrance is not None else "NA",
            round(r4_frac, 4) if r4_frac is not None else "NA", qc])
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample"] + labels + [
            "penetrance_(R1+R2)/(R1+R2+R6)", "r4_fraction_RNA-detected", "QC"])
        w.writerows(rows)
    print(f"summary -> {args.out}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("anchors")
    a.add_argument("--somatic-vcf", required=True)
    a.add_argument("--rna-vcf", required=True)
    a.add_argument("--out-tsv", required=True)
    a.add_argument("--out-bed", required=True)
    a.set_defaults(func=cmd_anchors)

    c = sub.add_parser("classify")
    c.add_argument("--sample", required=True)
    c.add_argument("--anchors", required=True)
    c.add_argument("--pileup", required=True)
    c.add_argument("--somatic-vcf", required=True)
    c.add_argument("--rna-vcf", required=True)
    c.add_argument("--germline-vcf", default=None)
    c.add_argument("--strain-vcf", nargs="*", default=None)
    c.add_argument("--splice-bed", nargs="*", default=None)
    c.add_argument("--editing-bed", nargs="*", default=None)
    c.add_argument("--min-af", type=float, default=0.02)
    c.add_argument("--min-ad", type=int, default=4)
    c.add_argument("--min-dp", type=int, default=20)
    c.add_argument("--germline-af", type=float, default=0.20)
    c.add_argument("--germline-dp", type=int, default=10)
    c.add_argument("--splice-window", type=int, default=5)
    c.add_argument("--out-dna", required=True)
    c.add_argument("--out-rna", required=True)
    c.set_defaults(func=cmd_classify)

    s = sub.add_parser("summary")
    s.add_argument("--tables", nargs="+", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--max-r4", type=float, default=0.30)
    s.set_defaults(func=cmd_summary)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
