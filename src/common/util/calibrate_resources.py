#!/usr/bin/env python
"""方案 D 校准脚本:benchmark tsv → resource_history.tsv → resources_suggest.tsv →(--apply) config.

用法(在 workflow 输出目录或其父目录下):
    python src/common/util/calibrate_resources.py --root OUTDIR [--root OUTDIR2 ...] \
        [--history PATH] [--config config/X.json] [--apply]
    多项目时可重复 --root 各项目输出目录,观测合并进同一历史;历史表位置由
    调用方以 --history 指定(不预设共享路径,多 --root 时必填)。

行为(对应 docs/resource-balancing-proposal.md §3.4):
  D.1 收集:扫描 --root 下所有 benchmark 输出——`**/benchmarks/*.tsv`
     或 Omics 实际布局 `**/*.benchmark.tsv` / `**/*.benchmark.txt`
     (snakemake `benchmark:` 输出,列 s / h:m:s / max_rss / max_vms,
     max_rss 单位为 MB),每 job 每 attempt 追加多行,**追加**进 --history
     (append-only;默认 <root>/resource_history.tsv,跨运行累积时显式
     传同一共享路径);
  规则名:benchmark 文件名即规则名;模块文件名与 config 资源键不一致的按
     RULE_ALIASES/SUFFIX_ALIASES 固定映射;PE/SE 规则共用同一 benchmark
     文件名的(bowtie2_align/trimming)按 raw.json 样本布局逐行归属,
     归属不了的行丢弃并计数(不猜);
  D.3 建议:按 rule 聚合历史,取 max_rss 的 p90 × 安全系数(默认 1.3)为建议值;
     ≥3 个数据点才建议;单次调整幅度限 ±30%(防震荡);OOM 行按删失处理
     ——只用其 max_rss 作下界,不压低建议;
  输出 --out(默认 <root>/resources_suggest.tsv):rule、n、observed_p90_mb、
     suggest_mb、prev_declared_mb、delta_pct、censored_n、note;
  --apply:把建议值写回 --config 的 resources.<rule>.mem_mb(默认只出表不写)。
     护栏:已有声明值的调整须 n>=MIN_POINTS(样本不足不动)且单次 ±30%
     限幅;config 里不存在的规则键 n>=1 新建锚定(首个真实数据集起步,
     家族聚合名除外,首次写入无 prev 不限幅);值取整、无变化不动。每次
     写回是一次可回滚的 config 变更,依据都在 tsv 里可审计。
"""
import argparse
import csv
import datetime
import glob
import json
import os
import re
import sys

DEFAULT_SAFETY = 1.3
DAMPING = 0.30          # 单次调整限幅 ±30%(D.4)
MIN_POINTS = 3          # ≥3 个数据点才动(D.4:样本不足不动)

# benchmark 文件名(去 .benchmark.*)→ config resources 键的固定映射
# (模块输出文件名与规则名不一致处,来源:各模块 .smk 的 benchmark: 声明)。
RULE_ALIASES = {
    "hisat2_build": "hisat2_index",
    "RNAseq_report": "generate_report",
    "gsva": "function_gsva",
    "cutadapt_paired": "trimming_Paired",
    "cutadapt_single": "trimming_Single",
    # benchmark 文件名与 rule_res 键的历史拼写差异,来源:各模块 benchmark: 声明
    "bwa_mem2-alignment": "bwaMem2_alignment",          # modules/bwa-mem2/bwa-mem2.smk
    "flagstat": "bam_flagstat",                          # modules/samtools/samtools.smk
    "haplotypeCaller": "HaplotypeCaller",                # modules/gatk/gatk_germline/gatk_germline.smk
    "haplotypeCaller-filtered": "filterHaplotypeCallerVcf",
    "DR_report": "generate_report",                      # modules/DR_report/DR_report.smk
    "telogator2": "telogator2_run",                      # length.smk 现名;存量观测为旧名
    # fibercluster fiber_cluster 的 {setting} 通配值档位,非规则名
    "fire_auto": "fiber_cluster",
    "msp_auto": "fiber_cluster",
}
# 文件名带 contrast 前缀(<contrast>.<suffix>)的按后缀映射。
SUFFIX_ALIASES = {
    ".go_kegg": "function_go_kegg",
    ".gsea": "function_gsea",
    ".manhattan": "DESeq2_TElocal_manhattan",
    ".locus_box": "DESeq2_TElocal_locus_boxplot",
    ".subfamily_bar": "DESeq2_TElocal_subfamily_bar",
    # 2026-11-26 拆开 polygenomes DESeq2 两条规则 benchmark 路径后的新命名
    ".TEcount": "DESeq2_TEcount",
    ".TElocal": "DESeq2_TElocal",
}
# PE/SE 两条规则共用同一 benchmark 文件名:按 raw.json 样本布局逐行归属。
LAYOUT_ALIASES = {
    "trimming": ("trimming_Paired", "trimming_Single"),
    "bowtie2_align": ("bowtie2_align_paired", "bowtie2_align_single"),
}


def parse_hms(v):
    """snakemake benchmark 'h:m:s' → seconds."""
    try:
        parts = [float(x) for x in str(v).split(":")]
    except ValueError:
        return None
    s = 0.0
    for p in parts:
        s = s * 60 + p
    return s


def load_sample_layout(root):
    """raw.json → {sample_id: "paired"/"single"};读不到返回 {}(归属不了的行丢弃)。"""
    for cand in ("raw.json", os.path.join("common", "raw.json")):
        p = os.path.join(root, cand)
        if not os.path.exists(p):
            continue
        try:
            with open(p) as fh:
                wf = (json.load(fh).get("Params") or {}).get("workflow") or {}
        except (OSError, ValueError):
            return {}
        layout = {}
        for sids in (wf.get("genome_paired_samples") or {}).values():
            for sid in sids or []:
                layout[str(sid)] = "paired"
        for sids in (wf.get("genome_single_samples") or {}).values():
            for sid in sids or []:
                layout[str(sid)] = "single"
        return layout
    return {}


def sample_from_path(path, layout):
    """benchmark 源路径里出现的样本 ID(路径分量与 raw.json 样本名精确匹配)。"""
    for part in path.split(os.sep):
        if part in layout:
            return part
    return None


def resolve_rule(base, path, layout):
    """benchmark 文件名 → config 资源键;归属不了返回 (None, 原因)。"""
    if base in RULE_ALIASES:
        return RULE_ALIASES[base], ""
    if base in LAYOUT_ALIASES:
        paired, single = LAYOUT_ALIASES[base]
        sid = sample_from_path(path, layout)
        if sid is None:
            return None, f"unattributable {base} row (no sample in path / raw.json)"
        return (paired if layout[sid] == "paired" else single), ""
    for suf, key in SUFFIX_ALIASES.items():
        if base.endswith(suf):
            return key, ""
    # 对比名文件名(<contrast>,如 2iL_vs_SL):主干是 comparison_id 不是规则名。
    # 规则取 log/group/<seg>/ 的模块段(如 DESeq2)——裸对比文件由 DESeq2_TEcount/
    # DESeq2_TElocal 共写同一 benchmark 路径,无法从文件区分,家族名是诚实标注。
    # (与 bioplatform OmicsRuleNames.resolve 同构,2026-11-26)
    parts = path.split(os.sep)
    seg = None
    for i in range(len(parts) - 1):
        if parts[i] == "group" and i + 1 < len(parts) - 1:
            seg = parts[i + 1]
            break
    if "_vs_" in base and seg:
        return seg, "family"
    # 已按家族名落库的历史行(如 "DESeq2"):与 group/<seg>/ 段同名 → 家族聚合,
    # 不是可写回 config 的规则键。
    if seg and base == seg:
        return base, "family"
    return base, ""


def read_benchmarks(root):
    """Collect (rule, attempt, seconds, max_rss_mb, source) from benchmark outputs.

    Accepts both layouts: benchmarks/<rule>.tsv and the Omics log layout
    <...>/<rule>.benchmark.tsv(.txt) (snakemake ``benchmark:`` output).
    """
    rows = []
    dropped = 0
    layout = load_sample_layout(root)
    seen_files = set()
    tsvs = []
    for pat in ("**/benchmarks/*.tsv", "**/*.benchmark.tsv", "**/*.benchmark.txt"):
        for tsv in sorted(glob.glob(os.path.join(root, pat), recursive=True)):
            if tsv not in seen_files:
                seen_files.add(tsv)
                tsvs.append(tsv)
    for tsv in tsvs:
        base = os.path.basename(tsv)
        if base.endswith(".benchmark.tsv"):
            base = base[:-len(".benchmark.tsv")]
        elif base.endswith(".benchmark.txt"):
            base = base[:-len(".benchmark.txt")]
        else:
            base = base[:-4]
        # OOM 删失标记:若某次运行检测到 OOM,wrapper 会写 benchmarks/<rule>.OOM;
        # 该规则本次全部观测按删失(真实需求 ≥ 观测)处理。
        oom_marker = os.path.join(os.path.dirname(tsv), base + ".OOM")
        censored_all = os.path.exists(oom_marker)
        try:
            with open(tsv, newline="") as fh:
                for rec in csv.DictReader(fh, delimiter="\t"):
                    # benchmark tsv 每 job(含每次重试)追加一行;attempt 无法从行号
                    # 唯一还原(多 wildcard job 共享文件),删失判定只用 .OOM 标记。
                    attempt = 1
                    rule, _why = resolve_rule(base, tsv, layout)
                    if rule is None:
                        dropped += 1
                        continue
                    secs = parse_hms(rec.get("s") or rec.get("h:m:s"))
                    try:
                        # snakemake benchmark 的 max_rss 单位是 MB
                        rss_mb = float(rec.get("max_rss", 0) or 0)
                    except ValueError:
                        rss_mb = 0.0
                    rows.append((rule, 1 if not censored_all else 2, secs, rss_mb, tsv))
        except OSError:
            continue
    if dropped:
        print(f"dropped {dropped} unattributable benchmark rows(不猜归属)", file=sys.stderr)
    return rows


def percentile(vals, q):
    if not vals:
        return None
    xs = sorted(vals)
    k = max(0, min(len(xs) - 1, int(round(q * (len(xs) - 1)))))
    return xs[k]


def detect_indent(text):
    """config json 的缩进宽度(写回时保持原格式)。"""
    m = re.search(r"\n( +)\"", text)
    return len(m.group(1)) if m else 2


def apply_to_config(config_path, applies):
    """把 {rule: mem_mb} 写回 config 的 resources(保持缩进与键序)。"""
    with open(config_path) as fh:
        text = fh.read()
    cfg = json.loads(text)
    res = cfg.setdefault("resources", {})
    for rule, new_mb in sorted(applies.items()):
        res.setdefault(rule, {})["mem_mb"] = new_mb
    indent = detect_indent(text)
    out = json.dumps(cfg, indent=indent, ensure_ascii=False)
    if text.endswith("\n"):
        out += "\n"
    with open(config_path, "w") as fh:
        fh.write(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", action="append", default=None, dest="roots", metavar="ROOT",
                    help="workflow 输出根目录(递归找 benchmark);可重复传多个项目 "
                         "(--root A --root B),各项目观测合并进同一历史")
    ap.add_argument("--history", default=None,
                    help="累计历史表(append-only;单 --root 时默认 <root>/resource_history.tsv;"
                         "多项目累积时须显式指定,位置由调用方决定,不预设共享路径)")
    ap.add_argument("--out", default=None,
                    help="建议表输出(默认 <root>/resources_suggest.tsv)")
    ap.add_argument("--safety", type=float, default=DEFAULT_SAFETY, help="安全系数(默认 1.3)")
    ap.add_argument("--config", default=None, help="workflow config json,读当前声明值列")
    ap.add_argument("--apply", action="store_true",
                    help="把建议值写回 --config 的 resources.<rule>.mem_mb"
                         f"(护栏:n<{MIN_POINTS} 不动、±30%% 限幅、取整、无变化不动)")
    ap.add_argument("--no-append", action="store_true", help="只出建议,不追加历史")
    args = ap.parse_args(argv)

    roots = args.roots or ["."]
    if args.history is None:
        if len(roots) > 1:
            ap.error("多 --root 须显式 --history(历史表位置由调用方指定,不预设共享路径)")
        args.history = os.path.join(roots[0], "resource_history.tsv")
    if args.out is None:
        args.out = os.path.join(roots[0], "resources_suggest.tsv")
    if args.apply and not args.config:
        ap.error("--apply 需要 --config")

    rows = []
    for r in roots:
        got = read_benchmarks(r)
        if not got:
            print(f"warning: no benchmark tsv under {r}", file=sys.stderr)
        rows += got
    if not rows:
        print("no benchmark tsv under any --root(先跑一次带 benchmark 的运行)", file=sys.stderr)
        return 1

    # ---- D.1 append to history ----
    new_history = 0
    if not args.no_append:
        exists = os.path.exists(args.history)
        seen = set()
        if exists:
            with open(args.history, newline="") as fh:
                for rec in csv.DictReader(fh, delimiter="\t"):
                    seen.add((rec["run_id"], rec["rule"], rec["attempt"], rec["max_rss_mb"]))
        run_id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        with open(args.history, "a", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            if not exists:
                w.writerow(["run_id", "ts", "rule", "attempt", "seconds", "max_rss_mb", "censored", "source"])
            ts = datetime.datetime.now().isoformat(timespec="seconds")
            for rule, attempt, secs, rss_mb, src in rows:
                key = (run_id, rule, attempt, rss_mb)
                # per-run dedupe: same (rule,attempt,rss) repeated across wildcard jobs
                if sum(1 for r in rows if r[0] == rule and r[1] == attempt and r[3] == rss_mb) > 1:
                    key = (run_id, rule, attempt, rss_mb, src)
                if key in seen:
                    continue
                w.writerow([run_id, ts, rule, attempt,
                            f"{secs:.1f}" if secs else "", f"{rss_mb:.0f}",
                            "1" if attempt > 1 else "0", src])
                new_history += 1
        print(f"history: appended {new_history} rows -> {args.history}")

    # ---- load history for suggestions ----
    # 历史行按当前别名表重新归并(append-only 不改写旧行):老的未映射名
    # (bwa_mem2-alignment 等)在此折叠进正确规则键,别名幂等。
    hist = {}
    family_rules = set()
    root_layouts = [(r, load_sample_layout(r)) for r in roots]

    def layout_for(src):
        """History rows attribute samples via the layout of the root they came from."""
        for r, lay in root_layouts:
            if src.startswith(r) or src.startswith(os.path.abspath(r)):
                return lay
        return {}

    if os.path.exists(args.history):
        with open(args.history, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                try:
                    rss = float(rec["max_rss_mb"])
                except (ValueError, KeyError):
                    continue
                censored = (rec.get("censored") == "1")
                rule, why = resolve_rule(rec["rule"], rec.get("source") or "",
                                         layout_for(rec.get("source") or ""))
                if rule is None:
                    continue
                if why == "family":
                    family_rules.add(rule)
                hist.setdefault(rule, {"all": [], "censored": []})
                (hist[rule]["censored"] if censored else hist[rule]["all"]).append(rss)

    declared = {}
    if args.config and os.path.exists(args.config):
        with open(args.config) as fh:
            declared = (json.load(fh).get("resources") or {})

    # ---- D.3 suggest ----
    out_rows = []
    applies = {}
    for rule in sorted(hist):
        obs = hist[rule]["all"]
        cen = hist[rule]["censored"]
        n_ok = len(obs)
        n_cen = len(cen)
        if n_ok == 0 and n_cen == 0:
            continue
        note = []
        low_n = n_ok + n_cen < MIN_POINTS
        if n_ok:
            p90 = percentile(obs, 0.90)
            suggest = p90 * args.safety
        else:
            # 只有删失观测:以下界给建议(真实需求 ≥ 观测)
            p90 = max(cen)
            suggest = p90 * args.safety
            note.append(f"censored-only({n_cen} oom rows), lower-bound suggestion")
        if n_cen:
            # 删失数据:建议值不得低于任何删失观测(下界约束)
            floor = max(cen) * args.safety
            if floor > suggest:
                suggest = floor
                note.append(f"censored floor applied({n_cen} oom rows)")
        prev = declared.get(rule, {}).get("mem_mb")
        if isinstance(prev, (int, float)):
            if low_n:
                note.append(f"insufficient points({n_ok + n_cen}<{MIN_POINTS}), suggestion informational")
            # 限幅:单次调整 ≤±30%(D.4)
            bounded = min(max(suggest, prev * (1 - DAMPING)), prev * (1 + DAMPING))
            if abs(bounded - suggest) > 1:
                note.append("damped ±30%")
                suggest = bounded
            delta = 100.0 * (suggest - prev) / prev
        else:
            delta = None
            if declared:
                if rule in family_rules:
                    note.append("family aggregate(家族聚合名,非规则键,不合入)")
                else:
                    note.append("not in config resources(--apply 新建条目)")
                    if low_n:
                        note.append(f"insufficient points({n_ok + n_cen}<{MIN_POINTS}), 单点锚定")
        out_rows.append({
            "rule": rule, "n_ok": n_ok, "n_censored": n_cen,
            "observed_p90_mb": f"{p90:.0f}",
            "suggest_mb": f"{suggest:.0f}",
            "prev_declared_mb": prev if prev else "",
            "delta_pct": f"{delta:+.0f}" if delta is not None else "",
            "note": "; ".join(note),
        })
        # ---- --apply 合入候选(护栏在写回前统一判定) ----
        if args.apply:
            new_mb = int(round(suggest))
            if isinstance(prev, (int, float)):
                # 已有声明值:调整须 n>=MIN_POINTS(样本不足不动),单次 ±30% 限幅
                if n_ok + n_cen >= MIN_POINTS and new_mb != int(prev):
                    applies[rule] = new_mb
            elif rule not in family_rules:
                # 新建锚定:首个真实数据集起步(n>=1,note 记录样本强度)。
                # 首次写入无 prev,不做 ±30% 限幅。
                applies[rule] = new_mb

    if out_rows:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(out_rows)
        print(f"suggest: {len(out_rows)} rules -> {args.out}")
    else:
        print("no rules with usable history")

    # ---- --apply 写回 config ----
    if args.apply:
        apply_to_config(args.config, applies)
        for rule in sorted(applies):
            old = declared.get(rule, {}).get("mem_mb")
            print(f"apply: {rule} {old if old else '(new)'} -> {applies[rule]}")
        print(f"apply: {len(applies)} rules written -> {args.config}"
              f"(n<{MIN_POINTS} 未动 {sum(1 for r in out_rows if 'insufficient' in r['note'])} 条)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
