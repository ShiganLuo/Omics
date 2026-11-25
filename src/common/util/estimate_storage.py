#!/usr/bin/env python
"""存储估算脚本:meta + config → 运行后磁盘占用(最终保留 / 运行期暂用 / 峰值)。

用法:
    # 跑之前估算(核心功能)
    python src/common/util/estimate_storage.py estimate \
        --meta data/meta_xxx.tsv --config config/PacVar.json --workflow PacVar

    # 跑完(或跑了一半)后实测:任意工作流,不需要 manifest
    python src/common/util/estimate_storage.py measure --root OUTDIR

    # 用真实运行校准比例表(append-only,和 calibrate_resources.py 同风格)
    python src/common/util/estimate_storage.py calibrate \
        --root OUTDIR --meta data/meta_xxx.tsv --workflow PacVar

行为:
  estimate
    1. meta 里每个输入文件真实 stat(不存在直接报错,不猜大小);
    2. 按 MANIFESTS[workflow] 的产物类清单展开 DAG 产物(含中间产物),
       每类大小 = 基准量(输入 bam / fastq / 基因组 / 固定值)× 比例系数;
    3. kept/temp 语义与模块一致:maybe_temp(...) 的产物在
       Params.<mod>.keep_intermediates == False 时为 temp,跑完自动删除;
    4. 峰值 = 按 stage 时间线求同时在盘的最大值。stage s 时在盘 =
       born <= s <= die 的类(temp 在其最后一个消费者所在 stage 仍占用,
       之后被 Snakemake 删除;kept 类 die=None 永久保留);
    5. 比例系数优先用 storage_ratio_history.tsv 里同 workflow 同类的
       实测校准(n>=MIN_POINTS 取 p50,p10/p90 为区间),否则用 manifest
       的经验值(low/mid/high 三档,输出区间而不是单点)。
    输出:分类明细表 + 最终保留/暂用/峰值三个数(带区间),--json 落盘。

  measure
    递归扫描 --root,按 manifest patterns 优先、文件后缀兜底归类,
    报告每类字节数/文件数/占比。软链接(输入数据)不计入,单独计数。

  calibrate
    measure + meta 输入基准 → 每类 实测/基准 比例,追加进
    storage_ratio_history.tsv(默认脚本同目录,--history 可改)。
    estimate 会自动读取该表,越跑越准。

原则(与 calibrate_resources.py 一致):
  - 不猜:输入 stat 不到就报错;没有 manifest 的工作流只做 measure/
    calibrate 后缀兜底,estimate 不凭空编造产物;
  - 区间而非单点:经验比例给 low/mid/high,校准后用实测分位数;
  - append-only 历史,estimate 的每次修正都有 tsv 依据可审计。

新增工作流:往 MANIFESTS 加一个产物类清单(字段见 ClassSpec 注释),
或用 --manifest extra.json 合并(格式: {"WorkflowName": [类 dict, ...]})。

注意(显式未建模,不在上表):
  - GATK/SnakeMake 的 java tmp-dir、.snakemake 元数据抖动:tmp-dir 在
    out/tmp,MarkDuplicates 期间可能额外写临时文件,结束后清理;
  - 文件系统块开销(小文件多时 du 略大于字节和)。
"""
import argparse
import csv
import datetime
import fnmatch
import json
import os
import re
import sys

MIN_POINTS = 3  # >=3 个校准点才覆盖经验比例(同 calibrate_resources.py)
HISTORY_NAME = "storage_ratio_history.tsv"

# ---------------------------------------------------------------------------
# ClassSpec(产物类清单,manifest 条目字段):
#   id      唯一 id,也用于校准表归并
#   desc    说明
#   scope   "sample" = 每样本一份(× n_samples); "run" = 全局一份
#   basis   "input_bam" | "input_fastq" | "input_total" | "genome_fasta" | "fixed"
#   ratio   (low, mid, high):basis 非 fixed 时为「输出字节/基准字节」比例;
#           basis == fixed 时为「每 scope 单位的绝对字节数」三档
#   gate    [(dotted_cfg_path, expected, default_if_missing), ...] 全部满足才
#           产生产物;default 与 node.py 语义一致(_enabled 缺省 True,
#           sawfish 缺省 False)
#   require ["genome_repeat_bed"] 等特殊前置:参考库里必须有对应文件
#   need_comparisons  True = 需要 meta design 里有 ctr_/exp_ 对(sv 分析类)
#   temp_when (dotted_cfg_path, expected, default):满足时该类是 temp
#           (跑完自动删);不满足按 kept
#   born    产生该类的 stage 序号; die  最后一个消费者所在 stage(temp
#           在该 stage 仍占用,之后释放);kept 类 die=None
#   patterns  measure/calibrate 归类用的相对路径 glob(先匹配先赢)
# ---------------------------------------------------------------------------

PACVAR_CLASSES = [
    dict(id="pbmm2_sorted_bam", desc="pbmm2 比对 BAM(裸 output,不标 temp,永久保留)",
         scope="sample", basis="input_bam", ratio=(0.95, 1.00, 1.10),
         gate=[], temp_when=None, born=2, die=None,
         patterns=["common/2_sorted_bam/**"]),
    dict(id="addReadsGroup_bam", desc="GATK 加 read group BAM(maybe_temp,keep_intermediates=false 时删除)",
         scope="sample", basis="input_bam", ratio=(0.95, 1.00, 1.10),
         gate=[], temp_when=("Params.gatk_prepare.keep_intermediates", False, True),
         born=3, die=4,
         patterns=["common/3_markdup_bam/**/*.addReadsGroup*"]),
    dict(id="markdup_bam", desc="MarkDuplicates 最终 BAM + metrics(交付主产物)",
         scope="sample", basis="input_bam", ratio=(0.95, 1.00, 1.10),
         gate=[], temp_when=None, born=4, die=None,
         patterns=["common/3_markdup_bam/**/*.sorted_markdup*"]),
    dict(id="qc_sexcheck", desc="性别检查 idxstats + 表",
         scope="run", basis="fixed", ratio=(1e6, 5e7, 5e8),
         gate=[], temp_when=None, born=4, die=None,
         patterns=["QC/sexcheck/**"]),
    dict(id="qc_summary", desc="QC 汇总表 + HTML",
         scope="run", basis="fixed", ratio=(1e6, 2e7, 2e8),
         gate=[], temp_when=None, born=4, die=None,
         patterns=["QC/qc_summary/**"]),
    dict(id="deepvariant_vcf", desc="DeepVariant SNV/indel VCF",
         scope="sample", basis="input_bam", ratio=(0.0002, 0.001, 0.003),
         gate=[("Params.snv.enabled", True, True),
               ("Params.workflow.snv_caller", "deepvariant", "deepvariant")],
         temp_when=None, born=5, die=None,
         patterns=["variation/germline_snv_indel/**/*.vcf.gz*"]),
    dict(id="deepvariant_examples", desc="DeepVariant make_examples 中间件(keep_intermediates=false 时删除)",
         scope="sample", basis="input_bam", ratio=(0.005, 0.02, 0.06),
         gate=[("Params.snv.enabled", True, True),
               ("Params.workflow.snv_caller", "deepvariant", "deepvariant")],
         temp_when=("Params.deepvariant.keep_intermediates", False, True),
         born=5, die=5,
         patterns=["variation/germline_snv_indel/**/make_examples*",
                   "variation/germline_snv_indel/**/*.tfrecord*"]),
    dict(id="pbsv_svsig", desc="pbsv discover 信号文件(registered deliverable,保留)",
         scope="sample", basis="fixed", ratio=(5e6, 3e7, 2e8),
         gate=[("Params.sv.enabled", True, True)],
         temp_when=None, born=5, die=None,
         patterns=["variation/germline_sv/**/*.svsig.gz*"]),
    dict(id="pbsv_vcf", desc="pbsv call SV VCF",
         scope="sample", basis="fixed", ratio=(2e6, 2e7, 1e8),
         gate=[("Params.sv.enabled", True, True)],
         temp_when=None, born=5, die=None,
         patterns=["variation/germline_sv/**/*.vcf.gz*"]),
    dict(id="hiphase_phased_vcf", desc="HiPhase 定相 VCF(SNV+SV)",
         scope="sample", basis="fixed", ratio=(1e7, 1e8, 5e8),
         gate=[("Params.phase.enabled", True, True),
               ("Params.snv.enabled", True, True),
               ("Params.sv.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["variation/**/*phased.vcf.gz*"]),
    dict(id="sawfish_joint", desc="sawfish 联合 CNV call(队列级)",
         scope="run", basis="fixed", ratio=(1e8, 5e8, 2e9),
         gate=[("Params.sawfish.enabled", True, False)],
         temp_when=None, born=6, die=None,
         patterns=["variation/germline_cnv/cnv/**"]),
    dict(id="sawfish_workdir", desc="sawfish per-sample 中间目录(keep_intermediates=false 时删除)",
         scope="sample", basis="fixed", ratio=(1e8, 1e9, 5e9),
         gate=[("Params.sawfish.enabled", True, False)],
         temp_when=("Params.sawfish.keep_intermediates", False, True),
         born=5, die=6,
         patterns=["variation/germline_cnv/**"]),
    dict(id="bamCoverage_bw", desc="覆盖度 bigWig(decoil/regionplot 共用)",
         scope="sample", basis="input_bam", ratio=(0.0005, 0.002, 0.006),
         gate=[],  # gate_any 注入:decoil 或 regionplot 任一开启
         temp_when=None, born=5, die=None,
         patterns=["eccdna/coverage/**"]),
    dict(id="sv_go_kegg", desc="SV 基因 GO/KEGG 分析",
         scope="run", basis="fixed", ratio=(1e7, 5e7, 2e8),
         gate=[("Params.function.sv_go_kegg.enabled", True, True),
               ("Params.sv.enabled", True, True)],
         need_comparisons=True, temp_when=None, born=6, die=None,
         patterns=["**/sv_go_kegg/**"]),
    dict(id="decoil_out", desc="Decoil eccDNA 重建 + 可视化(含 IGV HTML)",
         scope="sample", basis="fixed", ratio=(2e7, 1e8, 1e9),
         gate=[("Params.decoil.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["eccdna/**"]),
    dict(id="sv_analysis", desc="SV 比较分析(merged/annotated/MEI/RepeatMasker/circos)",
         scope="run", basis="fixed", ratio=(2e8, 1e9, 5e9),
         gate=[("Params.sv.enabled", True, True)],
         need_comparisons=True, temp_when=None, born=6, die=None,
         patterns=["variation/germline_sv_analysis/**"]),
    dict(id="trgt_out", desc="TRGT 重复序列分型(需 repeat_bed)",
         scope="sample", basis="fixed", ratio=(1e7, 1e8, 1e9),
         gate=[("Params.repeat.enabled", True, True)],
         require=["genome_repeat_bed"], temp_when=None, born=5, die=None,
         patterns=["repeat/trgt/**"]),
    dict(id="hifiasm_assembly", desc="hifiasm 组装目录(contigs + 未声明的 .bin/.gfa 中间件,一并算)",
         scope="sample", basis="input_bam", ratio=(0.05, 0.15, 0.30),
         gate=[("Params.telomere.enabled", True, True)],  # node.py:着丝粒/组装挂在 telomere 块内
         temp_when=None, born=2, die=None,
         patterns=["repeat/centromere/**/assembly/**"]),
    dict(id="centromere_out", desc="着丝粒分析(RepeatMasker + 统计)",
         scope="sample", basis="fixed", ratio=(1e7, 1e8, 1e9),
         gate=[("Params.telomere.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["repeat/centromere/**"]),
    dict(id="telomere_out", desc="端粒长度/组成/区域全套",
         scope="run", basis="fixed", ratio=(5e7, 5e8, 2e9),
         gate=[("Params.telomere.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["telomere/**"]),
    dict(id="regionplot_out", desc="区域图",
         scope="run", basis="fixed", ratio=(5e6, 5e7, 5e8),
         gate=[("Params.regionplot.enabled", True, True)],
         temp_when=None, born=6, die=None,
         patterns=["regionplot/**"]),
    dict(id="run_overhead", desc="raw.json/log/tmp/.snakemake/输入软链接目录等固定开销",
         scope="run", basis="fixed", ratio=(5e8, 2e9, 1e10),
         gate=[], temp_when=None, born=1, die=None,
         patterns=["log/**", "tmp/**", ".snakemake/**", "raw.json",
                   "common/raw.json", "common/1_raw_data/**"]),
]

# bamCoverage bigWig 的 gate 是 decoil 或 regionplot 任一开启(node.py:两者共用)
PACVAR_GATE_ANY = {
    "bamCoverage_bw": [("Params.decoil.enabled", True, True),
                       ("Params.regionplot.enabled", True, True)],
}

MANIFESTS = {"PacVar": PACVAR_CLASSES}

# measure/calibrate 的后缀兜底类(无 manifest 也可归类;先匹配 manifest 类)
FALLBACK_PATTERNS = [
    ("fb_bam", ["**/*.bam", "**/*.bam.bai", "**/*.bai", "**/*.cram", "**/*.crai", "**/*.pbi"]),
    ("fb_vcf", ["**/*.vcf", "**/*.vcf.gz", "**/*.vcf.gz.tbi", "**/*.vcf.gz.csi",
                "**/*.bcf", "**/*.bcf.csi"]),
    ("fb_bigwig", ["**/*.bw", "**/*.bigWig", "**/*.bigwig", "**/*.wig", "**/*.wig.gz"]),
    ("fb_fasta", ["**/*.fa", "**/*.fasta", "**/*.fa.gz", "**/*.fasta.gz",
                  "**/*.fa.fai", "**/*.fai", "**/*.gfa", "**/*.gfa.gz"]),
    ("fb_fastq", ["**/*.fastq", "**/*.fastq.gz", "**/*.fq", "**/*.fq.gz"]),
    ("fb_bed", ["**/*.bed", "**/*.bed.gz", "**/*.gtf", "**/*.gtf.gz",
                "**/*.gff", "**/*.gff3", "**/*.gff.gz"]),
    ("fb_table", ["**/*.tsv", "**/*.csv", "**/*.txt"]),
    ("fb_figure", ["**/*.png", "**/*.pdf", "**/*.svg", "**/*.html", "**/*.htm"]),
    ("fb_gz", ["**/*.gz", "**/*.gzip", "**/*.zip", "**/*.tar", "**/*.tgz"]),
]


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------

def glob_to_re(glob_pat):
    """glob → regex(支持 ** 跨目录;fnmatch 不支持 **,自己转)。"""
    out = []
    i = 0
    while i < len(glob_pat):
        c = glob_pat[i]
        if glob_pat.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif glob_pat.startswith("**", i):
            out.append(".*")
            i += 2
        elif c == "*":
            out.append("[^/]*")
            i += 1
        elif c == "?":
            out.append("[^/]")
            i += 1
        else:
            out.append(re.escape(c))
            i += 1
    return re.compile("^" + "".join(out) + "$")


def human(nbytes):
    """字节数 → 人类可读(十进制单位)。"""
    x = float(nbytes)
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(x) < 1000.0 or unit == "PB":
            return f"{x:.1f} {unit}" if unit != "B" else f"{x:.0f} B"
        x /= 1000.0


def cfg_get(cfg, dotted, default=None):
    obj = cfg
    for part in dotted.split("."):
        if isinstance(obj, dict) and part in obj:
            obj = obj[part]
        else:
            return default
    return obj


def parse_meta(path):
    """meta tsv/csv → 行 dict 列表(分隔符按扩展名,兜底嗅探)。"""
    delim = "," if path.lower().endswith(".csv") else "\t"
    with open(path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        if delim == "\t" and "\t" not in sample.splitlines()[0] and "," in sample.splitlines()[0]:
            delim = ","
        rows = list(csv.DictReader(fh, delimiter=delim))
    if not rows or "sample_id" not in rows[0]:
        raise SystemExit(f"meta 解析失败(缺 sample_id 列或空文件): {path}")
    return rows


INPUT_FILE_COLS = ("fastq_1", "fastq_2", "bam", "pbi", "bai", "ms_file")


def sample_inputs(rows, allow_missing=False):
    """meta → {sample_id: {basis: bytes}},输入文件真实 stat。

    basis: input_bam(bam 列) / input_fastq(fastq_1+2) / input_total(全部)。
    同一 sample 多行(data_id 一对多)累加。缺失文件默认报错(不猜)。
    """
    per = {}
    missing = []
    for row in rows:
        sid = (row.get("sample_id") or "").strip()
        if not sid:
            continue
        d = per.setdefault(sid, {"input_bam": 0, "input_fastq": 0, "input_total": 0,
                                 "n_files": 0})
        for col in INPUT_FILE_COLS:
            p = (row.get(col) or "").strip()
            if not p:
                continue
            if not os.path.exists(p):
                missing.append(p)
                continue
            size = os.path.getsize(p)  # 软链接按目标大小(真实数据量)
            d["input_total"] += size
            if col == "bam":
                d["input_bam"] += size
            elif col.startswith("fastq"):
                d["input_fastq"] += size
            d["n_files"] += 1
    if missing and not allow_missing:
        raise SystemExit("meta 输入文件不存在(不猜大小),共 %d 个,例如: %s"
                         % (len(missing), missing[0]))
    return per


def genome_size(cfg):
    """config genome.<default> 的 fasta 大小(0 = 未配置/不存在)。"""
    g = cfg.get("genome") or {}
    name = g.get("default")
    ref = (g.get("references") or {}).get(name) or {}
    fasta = ref.get("fasta")
    if fasta and os.path.exists(fasta):
        return os.path.getsize(fasta)
    return 0


def comparisons_from_meta(rows):
    """meta design 列 → ctr_/exp_ 比较对(同 MetaUtil.build_design_pairs 语义)。

    design 格式 ctr_TAG / ctrl_TAG / exp_TAG,TAG 是下划线分词集合;
    ctr 与 exp 共享任一 token 即成对(ctr_PlaBP6_DMSO ∩ exp_PlaBP6_PlaB = {PlaBP6})。
    返回 {对比名: {"ctr": [...], "exp": [...]}}。
    """
    pat = re.compile(r"^(ctr|ctrl|exp)_(.+)$")
    ctr_tags, exp_tags = {}, {}
    for row in rows:
        m = pat.match((row.get("design") or "").strip())
        if not m:
            continue
        role, tag = m.group(1), m.group(2)
        sid = (row.get("sample_id") or "").strip()
        d = ctr_tags if role in ("ctr", "ctrl") else exp_tags
        entry = d.setdefault(tag, {"tokens": set(tag.split("_")), "samples": []})
        entry["samples"].append(sid)
    out = {}
    for exp_tag, exp in exp_tags.items():
        for ctr_tag, ctr in ctr_tags.items():
            if exp["tokens"] & ctr["tokens"]:  # token 交集非空 → 成对
                out[f"{ctr_tag}_vs_{exp_tag}"] = {
                    "ctr": ctr["samples"], "exp": exp["samples"]}
    return out


# ---------------------------------------------------------------------------
# gate / 尺寸 / 校准
# ---------------------------------------------------------------------------

def gate_ok(cls, cfg, n_comparisons, gate_any=None):
    for cond in cls.get("gate", []):
        path, expected, default = cond
        if cfg_get(cfg, path, default) != expected:
            return False
    if gate_any:
        # any-of:一条满足即可,全不满足 → 不产生产物
        if not any(cfg_get(cfg, p, d) == e for p, e, d in gate_any):
            return False
    if cls.get("need_comparisons") and n_comparisons <= 0:
        return False
    for req in cls.get("require", []):
        if req == "genome_repeat_bed":
            g = cfg.get("genome") or {}
            ref = (g.get("references") or {}).get(g.get("default")) or {}
            if not ref.get("repeat_bed"):
                return False
        elif req == "genome_subtel_ref":
            g = cfg.get("genome") or {}
            ref = (g.get("references") or {}).get(g.get("default")) or {}
            if not ref.get("subtel_ref"):
                return False
    return True


def is_temp(cls, cfg):
    tw = cls.get("temp_when")
    if not tw:
        return False
    path, expected, default = tw
    return cfg_get(cfg, path, default) == expected


def load_history(path):
    """校准历史 → {(workflow, class_id): [ratio, ...]}。"""
    hist = {}
    if path and os.path.exists(path):
        with open(path, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                try:
                    hist.setdefault((rec["workflow"], rec["class_id"]), []).append(
                        float(rec["ratio"]))
                except (KeyError, ValueError):
                    continue
    return hist


def calibrated_ratio(hist, workflow, cls_id, manifest_lo, manifest_mid, manifest_hi):
    """有 >=MIN_POINTS 实测点 → p50/p10/p90,否则 manifest 经验值。"""
    obs = hist.get((workflow, cls_id)) or []
    if len(obs) < MIN_POINTS:
        return manifest_lo, manifest_mid, manifest_hi, len(obs), False
    xs = sorted(obs)

    def pct(q):
        k = max(0, min(len(xs) - 1, int(round(q * (len(xs) - 1)))))
        return xs[k]
    return pct(0.10), pct(0.50), pct(0.90), len(obs), True


# ---------------------------------------------------------------------------
# estimate
# ---------------------------------------------------------------------------

def do_estimate(args):
    with open(args.config) as fh:
        cfg = json.load(fh)
    if args.manifest:
        with open(args.manifest) as fh:
            MANIFESTS.update(json.load(fh))
    if args.workflow not in MANIFESTS:
        raise SystemExit(f"无 {args.workflow} 的产物清单(现有: {sorted(MANIFESTS)});"
                         "用 --manifest 提供,或先用 measure/calibrate 摸底")
    classes = MANIFESTS[args.workflow]
    rows = parse_meta(args.meta)
    per_sample = sample_inputs(rows, allow_missing=args.allow_missing)
    samples = sorted(per_sample)
    if not samples:
        raise SystemExit("meta 里没有可用样本")
    n_cmp = len(comparisons_from_meta(rows))
    gsize = genome_size(cfg)
    hist = load_history(args.calibration or
                        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     HISTORY_NAME))

    basis_all = {k: sum(d[k] for d in per_sample.values())
                 for k in ("input_bam", "input_fastq", "input_total")}

    table = []
    final = [0.0, 0.0, 0.0]
    transient = [0.0, 0.0, 0.0]
    timeline = []  # (born, die|None, lo, mid, hi),峰值时间线
    skipped = []
    for cls in classes:
        if not gate_ok(cls, cfg, n_cmp, gate_any=PACVAR_GATE_ANY.get(cls["id"])):
            skipped.append((cls["id"], cls["desc"]))
            continue
        # ratio 是「每单位」量:非 fixed = 输出/基准 比例;fixed = 每单位绝对
        # 字节。校准历史记的也是每单位量 → 先定每单位再乘基准,避免二次缩放。
        r_lo, r_mid, r_hi = cls["ratio"]
        r_lo, r_mid, r_hi, n_obs, used = calibrated_ratio(
            hist, args.workflow, cls["id"], r_lo, r_mid, r_hi)
        basis_key = cls["basis"] if cls["basis"] != "fixed" else None

        def sid_basis(sid):
            if basis_key == "genome_fasta":
                return gsize
            return per_sample[sid].get(basis_key, 0)

        if cls["scope"] == "sample":
            if basis_key:
                lo = sum(r_lo * sid_basis(s) for s in samples)
                mid = sum(r_mid * sid_basis(s) for s in samples)
                hi = sum(r_hi * sid_basis(s) for s in samples)
                unit_desc = f"{len(samples)} 样本 × 输入比例"
            else:
                lo, mid, hi = (r_lo * len(samples), r_mid * len(samples),
                               r_hi * len(samples))
                unit_desc = f"{len(samples)} 样本 × 固定"
        else:
            if basis_key:
                b = gsize if basis_key == "genome_fasta" else basis_all.get(basis_key, 0)
                lo, mid, hi = r_lo * b, r_mid * b, r_hi * b
                unit_desc = "全局 × 输入比例"
            else:
                lo, mid, hi = r_lo, r_mid, r_hi
                unit_desc = "全局固定"
        temp = is_temp(cls, cfg)
        for i in range(3):
            (transient if temp else final)[i] += (lo, mid, hi)[i]
        born = cls.get("born", 1)
        die = cls.get("die") if temp else None
        # kept 类 die=None:born 之后一直占用;temp 类占用到 die(含)
        timeline.append((born, die, lo, mid, hi))
        table.append(dict(
            id=cls["id"], desc=cls["desc"], scope=cls["scope"], unit=unit_desc,
            low=lo, mid=mid, hi=hi, temp=temp,
            calibrated=used, n_obs=n_obs))

    # 峰值:stage s 在盘 = born <= s 且 (die is None 或 s <= die) 的类之和
    checkpoints = set()
    for b, d, _lo, _mid, _hi in timeline:
        checkpoints.add(b)
        if d is not None:
            checkpoints.add(d)
            checkpoints.add(d + 1)
    peak_stage, peak = None, [0.0, 0.0, 0.0]
    for s in sorted(checkpoints):
        acc = [0.0, 0.0, 0.0]
        for b, d, lo_i, mid_i, hi_i in timeline:
            if b <= s and (d is None or s <= d):
                acc[0] += lo_i
                acc[1] += mid_i
                acc[2] += hi_i
        if acc[1] > peak[1]:
            peak_stage, peak = s, acc

    # ---- 报告 ----
    print(f"workflow: {args.workflow}   样本数: {len(samples)}   "
          f"SV 比较对: {n_cmp}   基因组 fasta: {human(gsize)}")
    print(f"输入合计: bam {human(basis_all['input_bam'])}, "
          f"fastq {human(basis_all['input_fastq'])}, "
          f"全部 {human(basis_all['input_total'])}")
    print()
    hdr = f"{'class':<22}{'scope':<8}{'kept/temp':<10}{'low':>10}{'mid':>10}{'high':>10}  note"
    print(hdr)
    print("-" * len(hdr))
    for r in table:
        note = f"校准 n={r['n_obs']}" if r["calibrated"] else f"经验比例 n={r['n_obs']}"
        print(f"{r['id']:<22}{r['scope']:<8}{'temp' if r['temp'] else 'kept':<10}"
              f"{human(r['low']):>10}{human(r['mid']):>10}{human(r['hi']):>10}  {note}")
    if skipped:
        print()
        print("本次未产生(被 config 关闭或缺前置):")
        for cid, desc in skipped:
            print(f"  - {cid}: {desc}")
    print()
    print("=" * 60)
    print(f"最终保留(跑完后长期占用): {human(final[1])}  "
          f"[{human(final[0])} ~ {human(final[2])}]")
    print(f"运行期暂用(temp,之后自动删): {human(transient[1])}  "
          f"[{human(transient[0])} ~ {human(transient[2])}]")
    print(f"峰值(stage {peak_stage},须按此预留): {human(peak[1])}  "
          f"[{human(peak[0])} ~ {human(peak[2])}]")
    print(f"                    (十进制 GB/GiB 换算: {final[1]/1e9:.1f} GB = "
          f"{final[1]/2**30:.1f} GiB)")
    print("=" * 60)
    print("未建模:GATK java tmp-dir 运行期临时写入(out/tmp,结束后清理)、"
          "小文件块开销。")

    if args.json:
        report = dict(
            workflow=args.workflow, generated=datetime.datetime.now().isoformat(
                timespec="seconds"),
            meta=args.meta, config=args.config,
            n_samples=len(samples), n_comparisons=n_cmp,
            input_bytes=basis_all, genome_fasta_bytes=gsize,
            classes=table, skipped=[c for c, _ in skipped],
            final_bytes=dict(zip(("low", "mid", "high"), final)),
            transient_bytes=dict(zip(("low", "mid", "high"), transient)),
            peak_bytes=dict(zip(("low", "mid", "high"), peak)),
            peak_stage=peak_stage,
            unmodeled=["gatk_java_tmpdir", "fs_block_overhead"])
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=2, ensure_ascii=False)
        print(f"json: {args.json}")
    return 0


# ---------------------------------------------------------------------------
# measure / calibrate
# ---------------------------------------------------------------------------

def classify_matcher(workflow=None):
    """(relpath) → class_id;manifest patterns 先赢,后缀兜底。"""
    entries = []
    if workflow and workflow in MANIFESTS:
        for cls in MANIFESTS[workflow]:
            for pat in cls["patterns"]:
                entries.append((glob_to_re(pat), cls["id"]))
    for cid, pats in FALLBACK_PATTERNS:
        for pat in pats:
            entries.append((glob_to_re(pat), cid))

    def match(relpath):
        for rx, cid in entries:
            if rx.match(relpath):
                return cid
        return "fb_other"
    return match


def walk_root(root, matcher):
    """扫描 root → {class: {"bytes": n, "files": n}} + 软链接计数。"""
    agg = {}
    n_links = 0
    for dirpath, dirnames, filenames in os.walk(root):
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            if os.path.islink(p):
                n_links += 1  # 输入软链接:不计目标大小,避免和真实数据重复计
                continue
            try:
                size = os.path.getsize(p)
            except OSError:
                continue
            rel = os.path.relpath(p, root)
            cid = matcher(rel)
            d = agg.setdefault(cid, {"bytes": 0, "files": 0})
            d["bytes"] += size
            d["files"] += 1
    return agg, n_links


def do_measure(args):
    matcher = classify_matcher(args.workflow)
    agg, n_links = walk_root(args.root, matcher)
    total = sum(d["bytes"] for d in agg.values())
    print(f"root: {args.root}")
    hdr = f"{'class':<22}{'files':>8}{'bytes':>14}{'size':>10}{'pct':>7}"
    print(hdr)
    print("-" * len(hdr))
    for cid, d in sorted(agg.items(), key=lambda kv: -kv[1]["bytes"]):
        pct = 100.0 * d["bytes"] / total if total else 0
        print(f"{cid:<22}{d['files']:>8}{d['bytes']:>14}{human(d['bytes']):>10}{pct:>6.1f}%")
    print("-" * len(hdr))
    print(f"{'TOTAL':<22}{sum(d['files'] for d in agg.values()):>8}{total:>14}"
          f"{human(total):>10}{100.0:>6.1f}%")
    if n_links:
        print(f"(软链接 {n_links} 个,不计入 —— 输入数据软链接不重复计)")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(dict(root=args.root, generated=datetime.datetime.now().isoformat(
                timespec="seconds"), total_bytes=total, symlinks_skipped=n_links,
                classes=agg), fh, indent=2, ensure_ascii=False)
        print(f"json: {args.json}")
    return 0


def do_calibrate(args):
    if not args.manifest and args.workflow not in MANIFESTS:
        raise SystemExit(f"无 {args.workflow} 的产物清单,无法把文件归到校准类;"
                         "用 --manifest 提供")
    if args.manifest:
        with open(args.manifest) as fh:
            MANIFESTS.update(json.load(fh))
    matcher = classify_matcher(args.workflow)
    agg, n_links = walk_root(args.root, matcher)
    rows = parse_meta(args.meta)
    per_sample = sample_inputs(rows, allow_missing=args.allow_missing)
    samples = sorted(per_sample)
    basis_all = {k: sum(d[k] for d in per_sample.values())
                 for k in ("input_bam", "input_fastq", "input_total")}
    gsize = genome_size(json.load(open(args.config)) if args.config else {})

    # 样本归属:路径分量与 meta sample_id 精确匹配(同 calibrate_resources)
    cls_by_id = {c["id"]: c for c in MANIFESTS[args.workflow]}
    history = args.history or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), HISTORY_NAME)

    # measure 的 bytes 是全 root 混在一起的;按样本重新扫一遍才能逐样本归因
    per_class_sample = {}  # (class_id, sample) -> bytes
    per_class_run = {}     # class_id -> bytes(非样本类)
    n_links = 0
    for dirpath, dirnames, filenames in os.walk(args.root):
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            if os.path.islink(p):
                n_links += 1
                continue
            try:
                size = os.path.getsize(p)
            except OSError:
                continue
            rel = os.path.relpath(p, args.root)
            cid = matcher(rel)
            sid = None
            for part in rel.split(os.sep):
                if part in per_sample:
                    sid = part
                    break
            cls = cls_by_id.get(cid)
            if cls is not None and cls["scope"] == "sample" and sid:
                per_class_sample[(cid, sid)] = per_class_sample.get((cid, sid), 0) + size
            else:
                per_class_run[cid] = per_class_run.get(cid, 0) + size

    run_id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    ts = datetime.datetime.now().isoformat(timespec="seconds")
    new_rows = 0
    exists = os.path.exists(history)
    with open(history, "a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        if not exists:
            w.writerow(["run_id", "ts", "workflow", "class_id", "scope",
                        "sample_id", "basis_bytes", "observed_bytes", "ratio", "source"])
        def emit(cid, scope, sid, basis_b, obs_b, cls):
            nonlocal new_rows
            if obs_b <= 0:
                return
            if cls is not None and cls["basis"] == "fixed":
                # fixed 类:校准「每单位绝对字节」,ratio 列记字节数
                ratio = obs_b
                basis_b = ""
            else:
                if basis_b <= 0:
                    return
                ratio = obs_b / basis_b
            w.writerow([run_id, ts, args.workflow, cid, scope, sid or "",
                        f"{basis_b:.0f}" if basis_b != "" else "",
                        f"{obs_b:.0f}", f"{ratio:.6g}", args.root])
            new_rows += 1

        for (cid, sid), obs in sorted(per_class_sample.items()):
            cls = cls_by_id.get(cid)
            basis_key = cls["basis"] if cls and cls["basis"] != "fixed" else None
            basis_b = per_sample[sid].get(basis_key, 0) if basis_key else 0
            emit(cid, "sample", sid, basis_b, obs, cls)
        for cid, obs in sorted(per_class_run.items()):
            cls = cls_by_id.get(cid)
            if cls is not None and cls["scope"] == "sample":
                continue  # 样本类没归因到样本的丢弃(不猜归属)
            basis_key = cls["basis"] if cls and cls["basis"] != "fixed" else None
            basis_b = basis_all.get(basis_key, 0) if basis_key else 0
            emit(cid, "run", None, basis_b, obs, cls)
    print(f"calibrate: appended {new_rows} rows -> {history}")
    print(f"(estimate 读该表需同类 >= {MIN_POINTS} 个点才覆盖经验比例)")
    return 0


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("estimate", help="跑之前估算占用")
    e.add_argument("--meta", required=True, help="meta tsv(csv),输入文件真实 stat")
    e.add_argument("--config", required=True, help="workflow config json(读 enabled/keep_intermediates)")
    e.add_argument("--workflow", required=True, help=f"工作流名(现有: {sorted(MANIFESTS)})")
    e.add_argument("--manifest", default=None, help="额外产物清单 json(合并进 MANIFESTS)")
    e.add_argument("--calibration", default=None,
                   help=f"校准历史 tsv(默认脚本同目录 {HISTORY_NAME})")
    e.add_argument("--json", default=None, help="报告另存 json")
    e.add_argument("--allow-missing", action="store_true",
                   help="输入文件缺失时不报错(该文件按 0 计,会低估)")
    e.set_defaults(func=do_estimate)

    m = sub.add_parser("measure", help="实测已有运行目录(任意工作流)")
    m.add_argument("--root", default=".", help="运行输出根目录")
    m.add_argument("--workflow", default=None, help="用该工作流的 manifest 精细归类(可选)")
    m.add_argument("--json", default=None, help="报告另存 json")
    m.set_defaults(func=do_measure)

    c = sub.add_parser("calibrate", help="实测 + 输入基准 → 追加校准历史")
    c.add_argument("--root", required=True, help="运行输出根目录")
    c.add_argument("--meta", required=True, help="该次运行的 meta tsv")
    c.add_argument("--config", default=None, help="config json(读基因组 fasta 作基准)")
    c.add_argument("--workflow", required=True, help="工作流名")
    c.add_argument("--manifest", default=None, help="额外产物清单 json")
    c.add_argument("--history", default=None,
                   help=f"校准历史 tsv(默认脚本同目录 {HISTORY_NAME},append-only)")
    c.add_argument("--allow-missing", action="store_true", help="meta 输入缺失不报错")
    c.set_defaults(func=do_calibrate)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
