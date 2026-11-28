#!/usr/bin/env python
"""benchmark 校准闭环脚本(对应 docs/resource-balancing-proposal.md §6):benchmark tsv → resource_history.tsv → resources_suggest.tsv →(--apply) config。

用法(在 workflow 输出目录或其父目录下):
    python src/common/util/calibrate_resources.py --root OUTDIR [--root OUTDIR2 ...] \
        [--history PATH] [--config config/X.json] [--apply]
    多项目时可重复 --root 各项目输出目录,观测合并进同一历史;历史表位置由
    调用方以 --history 指定(不预设共享路径,多 --root 时必填)。

行为(对应 docs/resource-balancing-proposal.md §6):
  6.1 收集:扫描 --root 下所有 benchmark 输出——`**/benchmarks/*.tsv`
     或 Omics 实际布局 `**/*.benchmark.tsv` / `**/*.benchmark.txt`
     (snakemake `benchmark:` 输出,列 s / h:m:s / max_rss / max_vms,
     max_rss 单位为 MB),每 job 每 attempt 追加多行,**追加**进 --history
     (append-only;默认 <root>/resource_history.tsv,跨运行累积时显式
     传同一共享路径);同一 benchmark 行按内容键去重,重复校准不重复追加;
  规则名:benchmark 文件名即规则名;模块文件名与 config 资源键不一致的按
     RULE_ALIASES/SUFFIX_ALIASES 固定映射;PE/SE 规则共用同一 benchmark
     文件名的(bowtie2_align/trimming)按 raw.json 样本布局逐行归属,
     归属不了的行丢弃并计数(不猜);
  规模(§6.2 拟合所需,与运行期 rule_res input_size 同口径):优先取
     `.snakemake/metadata` 里该 job 的真实输入文件(benchmark 路径与
     metadata 的 log 路径同干换后缀关联,次选 (log 主名, 样本) 匹配),
     按 bam / fastq 分型求和;取不到回退 `common/1_raw_data/<sample>/`
     原始链接规模;仍取不到则留空——该行不参与 scale 拟合,不猜。
     sample_count 口径记录当次样本数。历史表新增 input_size_mb /
     size_metric 两列(旧表自动迁移,旧行留空)。
  6.3 建议:静态规则 = 历史 max_rss p90 × 安全系数(默认 1.3);
     配 scale 的规则校准对象是整条曲线 mem = A + B×size:
     稳健仿射拟合 a+b×size(OOM 删失只抬升不压低),安全系数乘在拟合
     曲线上,基准点等值重算换回三参数——ref_input_size = 观测规模中位数,
     base_mem_mb = A + B×ref,fixed_frac = A / base;
     动 fixed_frac/ref 的门槛:n≥MIN_POINTS、≥2 个不同规模、跨度
     ≥ --fit-min-span(默认 3×,窄带斜率不可辨识);门槛不足时按旧
     scale 形状归一化只校 base。
  输出 --out(默认 <root>/resources_suggest.tsv):rule、n、observed_p90_mb、
     suggest_mb、prev_declared_mb、delta_pct、note,scale 拟合块
     (fit_size_metric、fit_n、size_min/max_mb、fit_a_mb、fit_b_mb_per_mb、
     suggest_ref_input_size_mb、suggest_fixed_frac、prev_ref_input_size_mb、
     prev_fixed_frac、curve_delta_max_pct),外加 io 列 io_n、
     observed_io_p90_mb_s、suggest_io_bandwidth_mb_s、io_delta_pct
     (设计 §10.1:io 权重 = p90((io_in+io_out)/s),±30% 限幅,只出建议
     人工合入)与 io 规模归一拟合列 io_fit_a_mb_s / io_fit_b_mb_s_per_mb
     (§8.3 输入跨数量级规则 io 配 scale 的备数据,只出建议)。
  --io-probe(默认关):实测 sustainable_disk_bandwidth_mb_s(设计 §10.2),
     优先现成工具 fio/hdparm,缺失时内置 256MB 测速;结果写 <probe-dir>/io_probe.json。
  --apply:把建议值写回 --config 的 resources.<rule>.mem_mb,拟合达标时同步
     写 scale.ref_input_size / scale.fixed_frac(默认只出表不写)。
     护栏:已有声明值的 mem 调整须 n>=MIN_POINTS(样本不足不动)且单次
     ±30% 限幅;config 里不存在的规则键 n>=1 新建锚定(首个真实数据集
     起步,家族聚合名除外,首次写入无 prev 不限幅);scale 更新/新建须拟合
     达标,且新旧曲线在各观测规模(含 ref)处 ±30% 限幅(不达标按二分向
     旧曲线回缩),OOM 删失点处新曲线 ≥ 删失 rss × 安全系数;新建 scale
     的 min/max 暂定 0.5×/4× base(出表注明,请人工审定);值取整、
     无变化不动。每次写回是一次可回滚的 config 变更,依据都在 tsv 里
     可审计。
"""
import argparse
import csv
import datetime
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

DEFAULT_SAFETY = 1.3
DAMPING = 0.30
MIN_POINTS = 3
SCALE_FIT_MIN_SPAN = 3.0
SCALE_NEW_MIN_MEM_FRAC = 0.5
SCALE_NEW_MAX_MEM_FRAC = 4.0
RAW_INDEX_SUFFIXES = (".pbi", ".idx", ".tbi", ".fai", ".gzi", ".crai",
                      ".csi", ".bai", ".dict")

RULE_ALIASES = {
    "hisat2_build": "hisat2_index",
    "RNAseq_report": "generate_report",
    "gsva": "function_gsva",
    "cutadapt_paired": "trimming_Paired",
    "cutadapt_single": "trimming_Single",
    "bwa_mem2-alignment": "bwaMem2_alignment",
    "flagstat": "bam_flagstat",
    "haplotypeCaller": "HaplotypeCaller",
    "haplotypeCaller-filtered": "filterHaplotypeCallerVcf",
    "DR_report": "generate_report",
    "telogator2": "telogator2_run",
    "fire_auto": "fiber_cluster",
    "msp_auto": "fiber_cluster",
}
SUFFIX_ALIASES = {
    ".go_kegg": "function_go_kegg",
    ".gsea": "function_gsea",
    ".manhattan": "DESeq2_TElocal_manhattan",
    ".locus_box": "DESeq2_TElocal_locus_boxplot",
    ".subfamily_bar": "DESeq2_TElocal_subfamily_bar",
    ".TEcount": "DESeq2_TEcount",
    ".TElocal": "DESeq2_TElocal",
}
LAYOUT_ALIASES = {
    "trimming": ("trimming_Paired", "trimming_Single"),
    "bowtie2_align": ("bowtie2_align_paired", "bowtie2_align_single"),
}

HIST_HEADER = ["run_id", "ts", "rule", "attempt", "seconds", "max_rss_mb",
               "io_mb_s", "censored", "source", "input_size_mb", "size_metric",
               "size_src"]


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


def load_raw_json(root):
    """raw.json 解析(候选 <root>/raw.json 与 <root>/common/raw.json);读不到返回 {}。"""
    for cand in ("raw.json", os.path.join("common", "raw.json")):
        p = os.path.join(root, cand)
        if not os.path.exists(p):
            continue
        try:
            with open(p) as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return {}
    return {}


def load_sample_layout(root):
    """raw.json → {sample_id: "paired"/"single"};读不到返回 {}(归属不了的行丢弃)。"""
    wf = (load_raw_json(root).get("Params") or {}).get("workflow") or {}
    layout = {}
    for sids in (wf.get("genome_paired_samples") or {}).values():
        for sid in sids or []:
            layout[str(sid)] = "paired"
    for sids in (wf.get("genome_single_samples") or {}).values():
        for sid in sids or []:
            layout[str(sid)] = "single"
    return layout


def load_sample_meta(root):
    """(sample_ids, n_samples, raw_sizes);raw_sizes = {sample: {"bam_mb", "fastq_mb"}}。

    样本清单来自 sample_params.samples / workflow.samples / genome_*_samples;
    原始输入规模来自 common/1_raw_data/<sample>/ 原始链接(次选 raw.json
    raw_files 文件名/路径精确匹配),取不到的样本不出现,不猜。
    """
    raw = load_raw_json(root)
    sp = raw.get("sample_params") or {}
    wf = (raw.get("Params") or {}).get("workflow") or {}
    ids = []
    for sids in (sp.get("samples"), wf.get("samples")):
        for sid in sids or []:
            if str(sid) not in ids:
                ids.append(str(sid))
    for key in ("genome_paired_samples", "genome_single_samples"):
        for sids in (wf.get(key) or {}).values():
            for sid in sids or []:
                if str(sid) not in ids:
                    ids.append(str(sid))
    raw_sizes = {}
    raw_dir = os.path.join(root, "common", "1_raw_data")
    if os.path.isdir(raw_dir):
        for sid in ids:
            sd = os.path.join(raw_dir, sid)
            if not os.path.isdir(sd):
                continue
            sums = _sum_files(os.path.join(sd, f) for f in os.listdir(sd)
                             if not f.lower().endswith(RAW_INDEX_SUFFIXES))
            if sums:
                raw_sizes[sid] = sums
    if len(raw_sizes) < len(ids):
        files = [f for f in (sp.get("raw_files") or wf.get("raw_files") or []) if f]
        for sid in ids:
            if sid in raw_sizes:
                continue
            mine = [f for f in files
                    if sid in os.path.basename(f) or sid in f.split(os.sep)]
            sums = _sum_files(mine)
            if sums:
                raw_sizes[sid] = sums
    return ids, len(ids), raw_sizes


def _classify_input(path):
    name = os.path.basename(path).lower()
    if name.endswith(".bam"):
        return "bam"
    for ext in (".fastq", ".fq"):
        for comp in ("", ".gz", ".bz2"):
            if name.endswith(ext + comp):
                return "fastq"
    return "other"


def _sum_files(paths):
    """{"bam_mb", "fastq_mb"}(MB,0 记 None);一个文件都取不到返回 None。"""
    bam = fastq = 0.0
    seen = 0
    for p in paths:
        try:
            sz = os.path.getsize(p)
        except OSError:
            continue
        seen += 1
        cls = _classify_input(p)
        if cls == "bam":
            bam += sz / 1e6
        elif cls == "fastq":
            fastq += sz / 1e6
    if not seen:
        return None
    return {"bam_mb": bam or None, "fastq_mb": fastq or None}


def load_job_inputs(root):
    """(by_log, by_stem):benchmark 路径 → 该 job 真实输入规模(.snakemake/metadata)。

    by_log 键 = log 路径换 .benchmark.tsv/.benchmark.txt 后缀的绝对路径;
    by_stem 键 = (log 文件主名, 路径中的样本 ID)。值 = {"bam_mb", "fastq_mb"}。
    同键指向不同输入的 job 标记 None(歧义,不猜)。
    """
    by_log = {}
    by_stem = {}
    ids, _, _ = load_sample_meta(root)
    meta_dir = os.path.join(root, ".snakemake", "metadata")
    for mf in glob.glob(os.path.join(meta_dir, "*")):
        try:
            with open(mf) as fh:
                meta = json.load(fh)
        except (OSError, ValueError):
            continue
        inputs = meta.get("input") or []
        if isinstance(inputs, dict):
            inputs = [v for vs in inputs.values()
                      for v in (vs if isinstance(vs, list) else [vs])]
        sums = _sum_files(str(p) for p in inputs if p)
        if not sums:
            continue
        for lp in meta.get("log") or []:
            lp = str(lp)
            if not lp.endswith(".log"):
                continue
            stem = os.path.basename(lp)[:-4]
            key = (stem, sample_from_path(lp, ids))
            if key in by_stem and by_stem[key] != sums:
                by_stem[key] = None
            else:
                by_stem.setdefault(key, sums)
            for suf in (".benchmark.tsv", ".benchmark.txt"):
                k = os.path.abspath(lp[:-4] + suf)
                if k in by_log and by_log[k] != sums:
                    by_log[k] = None
                else:
                    by_log.setdefault(k, sums)
    return by_log, by_stem


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
    parts = path.split(os.sep)
    seg = None
    for i in range(len(parts) - 1):
        if parts[i] == "group" and i + 1 < len(parts) - 1:
            seg = parts[i + 1]
            break
    if "_vs_" in base and seg:
        return seg, "family"
    if seg and base == seg:
        return base, "family"
    return base, ""


def _row_size(metric, sums, tsv, ids, raw_sizes, n_samples, size_src):
    """单行 (input_size_mb, size_metric, size_src);取不到 → (None, "", "")。"""
    if metric == "sample_count":
        return (n_samples, "sample_count", "job_input") if n_samples else (None, "", "")
    src = size_src
    if sums is None:
        sums = raw_sizes.get(sample_from_path(tsv, ids)) or {}
        src = "raw" if sums else ""
    bam = (sums or {}).get("bam_mb")
    fastq = (sums or {}).get("fastq_mb")
    if metric == "bam_mb":
        return (bam, "bam_mb", src) if bam else (None, "", "")
    if metric == "fastq_mb":
        return (fastq, "fastq_mb", src) if fastq else (None, "", "")
    if bam and fastq:
        return (bam, "bam_mb", src) if bam >= fastq else (fastq, "fastq_mb", src)
    if bam:
        return bam, "bam_mb", src
    if fastq:
        return fastq, "fastq_mb", src
    return None, "", ""


def read_benchmarks(root, declared=None):
    """Collect (rule, attempt, seconds, max_rss_mb, io_mb_s, size, size_metric,
    size_src, source) from benchmark outputs.

    Accepts both layouts: benchmarks/<rule>.tsv and the Omics log layout
    <...>/<rule>.benchmark.tsv(.txt) (snakemake ``benchmark:`` output).
    io_mb_s = (io_in + io_out) / wall seconds (设计 §8.3); None when the tsv
    has no io columns (old benchmarks) so zeros never fake "no I/O".
    size/size_metric = per-row input size in the rule's scale metric
    (设计 §6.2,与运行期 input_size 同口径;size_src = job_input/raw/"" 为
    来源标记);metric 取该规则声明的 scale.size_metric,无声明时取 job
    输入里 bam/fastq 的大者;取不到规模的行 size=None,不猜。
    """
    rows = []
    dropped = 0
    declared = declared or {}
    layout = load_sample_layout(root)
    ids, n_samples, raw_sizes = load_sample_meta(root)
    by_log, by_stem = load_job_inputs(root)
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
        oom_marker = os.path.join(os.path.dirname(tsv), base + ".OOM")
        censored_all = os.path.exists(oom_marker)
        sums = by_log.get(os.path.abspath(tsv))
        size_src = "job_input" if sums else ""
        if sums is None:
            sums = by_stem.get((base, sample_from_path(tsv, ids)))
            size_src = "job_input" if sums else ""
        try:
            with open(tsv, newline="") as fh:
                for rec in csv.DictReader(fh, delimiter="\t"):
                    attempt = 1
                    rule, _why = resolve_rule(base, tsv, layout)
                    if rule is None:
                        dropped += 1
                        continue
                    secs = parse_hms(rec.get("s") or rec.get("h:m:s"))
                    try:
                        rss_mb = float(rec.get("max_rss", 0) or 0)
                    except ValueError:
                        rss_mb = 0.0
                    if "io_in" in rec or "io_out" in rec:
                        try:
                            io_mb = float(rec.get("io_in") or 0) + float(rec.get("io_out") or 0)
                        except ValueError:
                            io_mb = 0.0
                        io_rate = (io_mb / secs) if (secs and secs > 0) else None
                    else:
                        io_rate = None
                    metric = ((declared.get(rule) or {}).get("scale") or {}).get("size_metric")
                    size, size_label, src = _row_size(
                        metric, sums, tsv, ids, raw_sizes, n_samples, size_src)
                    rows.append((rule, 1 if not censored_all else 2, secs, rss_mb,
                                 io_rate, size, size_label, src, tsv))
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


def fit_affine(xs, ys):
    """least-squares a + b×x;返回 (a, b),不可拟合返回 None。"""
    n = len(xs)
    if n < 2:
        return None
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if denom == 0:
        return None
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    return a, b


def curve_ok(a, b, xs, old_a, old_b):
    """新曲线在各 x 处 ∈ [0.7, 1.3]×旧曲线。"""
    for x in xs:
        y = a + b * x
        old = old_a + old_b * x
        if old > 0 and not (0.7 * old - 1 <= y <= 1.3 * old + 1):
            return False
    return True


def blend_curve(a_new, b_new, old_a, old_b, xs):
    """二分最大 keep∈[0,1],使 keep×新+(1-keep)×旧满足曲线护栏;返回 (a, b, keep)。"""
    if curve_ok(a_new, b_new, xs, old_a, old_b):
        return a_new, b_new, 1.0
    lo, hi = 0.0, 1.0
    for _ in range(20):
        mid = (lo + hi) / 2
        a = old_a + mid * (a_new - old_a)
        b = old_b + mid * (b_new - old_b)
        if curve_ok(a, b, xs, old_a, old_b):
            lo = mid
        else:
            hi = mid
    keep = lo
    return old_a + keep * (a_new - old_a), old_b + keep * (b_new - old_b), keep


def detect_indent(text):
    """config json 的缩进宽度(写回时保持原格式)。"""
    m = re.search(r"\n( +)\"", text)
    return len(m.group(1)) if m else 2


def apply_to_config(config_path, applies):
    """把 {rule: {"mem_mb", "scale", "create_scale"}} 写回 config 的 resources(保持缩进与键序)。"""
    with open(config_path) as fh:
        text = fh.read()
    cfg = json.loads(text)
    res = cfg.setdefault("resources", {})
    for rule in sorted(applies):
        upd = applies[rule]
        entry = res.setdefault(rule, {})
        if upd.get("mem_mb") is not None:
            entry["mem_mb"] = upd["mem_mb"]
        sc = upd.get("scale")
        if sc:
            if upd.get("create_scale"):
                entry["scale"] = dict(sc)
            else:
                entry.setdefault("scale", {}).update(sc)
    indent = detect_indent(text)
    out = json.dumps(cfg, indent=indent, ensure_ascii=False)
    if text.endswith("\n"):
        out += "\n"
    with open(config_path, "w") as fh:
        fh.write(out)


def _median(vals):
    xs = sorted(vals)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def _majority_metric(pts):
    counts = {}
    for _x, _r, _c, label, _io, _src in pts:
        if label:
            counts[label] = counts.get(label, 0) + 1
    return max(counts, key=counts.get) if counts else ""


def probe_disk_bandwidth(probe_dir, size_mb=256, rounds=3):
    """sustainable_disk_bandwidth_mb_s 实测(设计 §10.2)。

    优先现成 Linux 测速工具(2026-11-27 已决):fio 顺序写+读 → hdparm(只读
    上界估计)→ 内置 256 MB 写+读回退;每轮取写/读较慢者(混合读写口径),
    3 轮取中位。返回 (mb_s, method)。测速本身写盘,只在空闲时人工开启。
    """
    fn = os.path.join(probe_dir, ".io_probe.bin")

    def _cleanup():
        try:
            os.remove(fn)
        except OSError:
            pass

    if shutil.which("fio"):
        rates = []
        try:
            for rw in ("write", "read"):
                out = subprocess.run(
                    ["fio", "--name=io_probe", f"--filename={fn}",
                     f"--size={size_mb}M", f"--rw={rw}", "--bs=1M",
                     "--ioengine=psync", "--direct=1", "--numjobs=1",
                     "--output-format=json"],
                    capture_output=True, text=True, check=True)
                data = json.loads(out.stdout)
                rates.append(data["jobs"][0][rw]["bw_bytes"] / 1e6)
            _cleanup()
            return min(rates), "fio(write+read)"
        except (subprocess.CalledProcessError, KeyError, ValueError, OSError,
                IndexError) as e:
            print(f"io-probe: fio failed ({e}), falling back", file=sys.stderr)
            _cleanup()

    if shutil.which("hdparm"):
        try:
            st_dev = os.stat(probe_dir).st_dev
            dev_dir = os.path.realpath(
                f"/sys/dev/block/{os.major(st_dev)}:{os.minor(st_dev)}")
            devname = None
            with open(os.path.join(dev_dir, "uevent")) as fh:
                for line in fh:
                    if line.startswith("DEVNAME="):
                        devname = line.strip().split("=", 1)[1]
                        break
            if devname:
                out = subprocess.run(["hdparm", "-t", f"/dev/{devname}"],
                                     capture_output=True, text=True, check=True)
                m = re.search(r"=\s*([\d.]+)\s*MB/sec", out.stdout)
                if m:
                    return float(m.group(1)), "hdparm(read-only 上界估计)"
        except (subprocess.CalledProcessError, OSError, ValueError) as e:
            print(f"io-probe: hdparm failed ({e}), falling back", file=sys.stderr)

    # 内置回退:256 MB 顺序写+读,3 轮取中位
    rates = []
    block = os.urandom(1 << 20)  # 1 MB,避免全零被压缩/稀疏化
    try:
        for _ in range(rounds):
            t0 = time.monotonic()
            with open(fn, "wb") as fh:
                for _ in range(size_mb):
                    fh.write(block)
                fh.flush()
                os.fsync(fh.fileno())
            w_rate = size_mb / max(time.monotonic() - t0, 1e-9)
            t0 = time.monotonic()
            with open(fn, "rb") as fh:
                while fh.read(1 << 20):
                    pass
            r_rate = size_mb / max(time.monotonic() - t0, 1e-9)
            rates.append(min(w_rate, r_rate))
    finally:
        _cleanup()
    return _median(rates), "builtin(write+read, 256MB×3 取中位)"


def _content_key(rule, attempt, rss, secs, src):
    return (rule, str(attempt), f"{rss:.0f}", f"{secs:.1f}" if secs else "", src)


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
    ap.add_argument("--fit-min-span", type=float, default=SCALE_FIT_MIN_SPAN,
                    help="动 fixed_frac/ref_input_size 的最小规模跨度(默认 3×;窄带斜率不可辨识,"
                         "门槛内只校 base)")
    ap.add_argument("--config", default=None, help="workflow config json,读当前声明值列")
    ap.add_argument("--apply", action="store_true",
                    help="把建议值写回 --config 的 resources.<rule>.mem_mb,拟合达标时同步写 "
                         "scale.ref_input_size / scale.fixed_frac"
                         f"(护栏:n<{MIN_POINTS} 不动、±30%% 限幅、取整、无变化不动)")
    ap.add_argument("--no-append", action="store_true", help="只出建议,不追加历史")
    ap.add_argument("--io-probe", action="store_true",
                    help="实测 sustainable_disk_bandwidth_mb_s(设计 §10.2,默认关;优先 "
                         "fio/hdparm,缺失时内置 256MB 测速;测速写盘,只在空闲时人工开启)")
    ap.add_argument("--probe-dir", default=None,
                    help="--io-probe 的测速目录(默认第一个 --root)")
    args = ap.parse_args(argv)

    roots = args.roots or ["."]

    if args.io_probe:
        probe_dir = os.path.abspath(args.probe_dir or roots[0])
        mb_s, method = probe_disk_bandwidth(probe_dir)
        print(f"io-probe: sustainable_disk_bandwidth_mb_s = {mb_s:.0f} "
              f"({method}; dir={probe_dir})")
        try:
            with open(os.path.join(probe_dir, "io_probe.json"), "w") as fh:
                json.dump({"sustainable_disk_bandwidth_mb_s": mb_s,
                           "method": method, "dir": probe_dir,
                           "ts": datetime.datetime.now().isoformat(timespec="seconds")},
                          fh, indent=2)
        except OSError:
            pass
    if args.history is None:
        if len(roots) > 1:
            ap.error("多 --root 须显式 --history(历史表位置由调用方指定,不预设共享路径)")
        args.history = os.path.join(roots[0], "resource_history.tsv")
    if args.out is None:
        args.out = os.path.join(roots[0], "resources_suggest.tsv")
    if args.apply and not args.config:
        ap.error("--apply 需要 --config")

    declared = {}
    if args.config and os.path.exists(args.config):
        with open(args.config) as fh:
            declared = (json.load(fh).get("resources") or {})

    rows = []
    for r in roots:
        got = read_benchmarks(r, declared)
        if not got:
            print(f"warning: no benchmark tsv under {r}", file=sys.stderr)
        rows += got
    if not rows:
        print("no benchmark tsv under any --root(先跑一次带 benchmark 的运行)", file=sys.stderr)
        return 1

    new_history = 0
    if not args.no_append:
        exists = os.path.exists(args.history)
        seen = set()
        if exists:
            with open(args.history, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                old_rows = list(reader)
                old_header = reader.fieldnames or []
            if old_header != HIST_HEADER:
                with open(args.history, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=HIST_HEADER, delimiter="\t")
                    w.writeheader()
                    for rec in old_rows:
                        w.writerow({k: rec.get(k, "") for k in HIST_HEADER})
            for rec in old_rows:
                try:
                    seen.add(_content_key(rec["rule"], rec.get("attempt", "1"),
                                          float(rec["max_rss_mb"]),
                                          float(rec["seconds"]) if rec.get("seconds") else None,
                                          rec.get("source") or ""))
                except (ValueError, KeyError):
                    continue
        run_id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        with open(args.history, "a", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            if not exists:
                w.writerow(HIST_HEADER)
            ts = datetime.datetime.now().isoformat(timespec="seconds")
            for rule, attempt, secs, rss_mb, io_rate, size, size_label, _src, src in rows:
                key = _content_key(rule, attempt, rss_mb, secs, src)
                if key in seen:
                    continue
                w.writerow([run_id, ts, rule, attempt,
                            f"{secs:.1f}" if secs else "", f"{rss_mb:.0f}",
                            f"{io_rate:.1f}" if io_rate is not None else "",
                            "1" if attempt > 1 else "0", src,
                            f"{size:.1f}" if size is not None else "", size_label,
                            _src])
                new_history += 1
        print(f"history: appended {new_history} rows -> {args.history}")

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
                hist.setdefault(rule, {"all": [], "censored": [], "io": [], "pts": []})
                (hist[rule]["censored"] if censored else hist[rule]["all"]).append(rss)
                raw_io = rec.get("io_mb_s")
                io_val = None
                if raw_io not in (None, ""):
                    try:
                        io_val = float(raw_io)
                    except ValueError:
                        io_val = None
                if io_val is not None and not censored:
                    hist[rule]["io"].append(io_val)
                size = None
                if rec.get("input_size_mb") not in (None, ""):
                    try:
                        size = float(rec["input_size_mb"])
                    except ValueError:
                        size = None
                if size is not None:
                    hist[rule]["pts"].append(
                        (size, rss, censored, rec.get("size_metric") or "", io_val,
                         rec.get("size_src") or ""))

    out_rows = []
    applies = {}
    for rule in sorted(hist):
        obs = hist[rule]["all"]
        cen = hist[rule]["censored"]
        pts = hist[rule]["pts"]
        n_ok = len(obs)
        n_cen = len(cen)
        if n_ok == 0 and n_cen == 0:
            continue
        note = []
        low_n = n_ok + n_cen < MIN_POINTS
        entry = declared.get(rule) if isinstance(declared.get(rule), dict) else {}
        sc = entry.get("scale") if isinstance(entry.get("scale"), dict) else None
        prev = entry.get("mem_mb")
        prev_num = isinstance(prev, (int, float))
        p90 = percentile(obs, 0.90) if n_ok else max(cen)
        scale_rows = {"fit_size_metric": "", "fit_n": 0, "size_min_mb": "",
                      "size_max_mb": "", "fit_a_mb": "", "fit_b_mb_per_mb": "",
                      "suggest_ref_input_size_mb": "", "suggest_fixed_frac": "",
                      "prev_ref_input_size_mb": "", "prev_fixed_frac": "",
                      "curve_delta_max_pct": ""}

        metric = (sc or {}).get("size_metric") or _majority_metric(pts)
        pts_m = [p for p in pts if p[3] == metric]
        excluded = len([p for p in pts if p[3] and p[3] != metric])
        if excluded:
            note.append(f"{excluded} 行规模口径与目标 metric({metric})不符,未参与拟合")
        src_counts = {}
        for p in pts_m:
            src_counts[p[5]] = src_counts.get(p[5], 0) + 1
        if len(src_counts) > 1:
            fit_src = "job_input" if src_counts.get("job_input") \
                else max(src_counts, key=src_counts.get)
            note.append("规模口径来源混杂("
                        + "/".join(f"{k or 'unknown'}:{v}" for k, v in src_counts.items())
                        + f"),拟合仅用 {fit_src or 'unknown'}")
        else:
            fit_src = next(iter(src_counts), "")
        pts_m = [p for p in pts_m if p[5] == fit_src] or pts_m
        fit_pts = [p for p in pts_m if not p[2]]
        xs = [p[0] for p in fit_pts]
        ys = [p[1] for p in fit_pts]
        distinct = sorted(set(xs))
        span = (max(distinct) / min(distinct)) if len(distinct) >= 2 and min(distinct) > 0 else 1.0
        fit = fit_affine(xs, ys) if len(xs) >= 2 else None
        fit_ready = bool(fit) and len(fit_pts) >= MIN_POINTS and len(distinct) >= 2 \
            and span >= args.fit_min_span
        floor_pts = [(p[0], p[1] * args.safety) for p in pts_m if p[2]]
        a = b = None
        A_raw = B_raw = 0.0
        if fit:
            a, b = max(fit[0], 0.0), max(fit[1], 0.0)
            A_raw, B_raw = a * args.safety, b * args.safety
            for x, floor in floor_pts:
                deficit = floor - (A_raw + B_raw * x)
                if deficit > 0:
                    A_raw += deficit
                    note.append(f"censored floor applied at size {x:.0f}")
            scale_rows["fit_size_metric"] = metric
            scale_rows["fit_n"] = len(fit_pts)
            scale_rows["size_min_mb"] = f"{min(distinct):.0f}"
            scale_rows["size_max_mb"] = f"{max(distinct):.0f}"
            scale_rows["fit_a_mb"] = f"{a:.0f}"
            scale_rows["fit_b_mb_per_mb"] = f"{b:.5f}"

        def curve_delta(a, b, guard_xs, old_a, old_b):
            deltas = [100.0 * ((a + b * x) / (old_a + old_b * x) - 1)
                      for x in guard_xs if (old_a + old_b * x) > 0]
            return f"{max(deltas, key=abs):+.0f}" if deltas else ""

        suggest = None
        scale_upd = None
        create_scale = False
        damped_value = True
        if sc and not prev_num:
            note.append("scale 存在但无 mem_mb 声明,跳过曲线校准(仅静态建议)")
        if sc and prev_num:
            base_old = float(prev)
            ref_old = float(sc.get("ref_input_size") or 0)
            frac_old = float(sc.get("fixed_frac") or 0)
            scale_rows["prev_ref_input_size_mb"] = f"{ref_old:.0f}" if ref_old else ""
            scale_rows["prev_fixed_frac"] = f"{frac_old:.3f}"
            if fit_ready:
                old_a = base_old * frac_old
                old_b = base_old * (1 - frac_old) / ref_old if ref_old > 0 else 0.0
                guard_xs = list(distinct) + ([ref_old] if ref_old > 0 else [])
                A, B, keep = blend_curve(A_raw, B_raw, old_a, old_b, guard_xs)
                if keep < 1.0:
                    note.append(f"scale fit damped to {keep:.0%} keep(曲线护栏 ±30%)")
                unmet = [(x, f) for x, f in floor_pts if A + B * x < f - 1]
                if unmet:
                    note.append(f"censored floor 未达({len(unmet)} 点,±30% 限幅内不可及,"
                                f"连续同向后放宽)")
                ref_new = _median(xs)
                base_new = max(A + B * ref_new, 1.0)
                frac_new = min(max(A / base_new, 0.0), 1.0)
                suggest = base_new
                scale_upd = {"ref_input_size": int(round(ref_new)),
                             "fixed_frac": round(frac_new, 3)}
                scale_rows["suggest_ref_input_size_mb"] = f"{ref_new:.0f}"
                scale_rows["suggest_fixed_frac"] = f"{frac_new:.3f}"
                scale_rows["curve_delta_max_pct"] = curve_delta(A, B, guard_xs, old_a, old_b)
                damped_value = False
            else:
                factor = (lambda x: frac_old + (1 - frac_old) * x / ref_old) \
                    if ref_old > 0 else (lambda x: 1.0)
                vals = [rss / max(factor(x), 1e-9) for x, rss, _c, _m, _io, _s in fit_pts]
                if vals:
                    suggest = percentile(vals, 0.90) * args.safety
                if floor_pts:
                    floor = max(f / max(factor(x), 1e-9) for x, f in floor_pts)
                    if suggest is None or floor > suggest:
                        suggest = floor
                        note.append(f"censored floor applied({n_cen} oom rows)")
                note.append(f"scale shape kept(span {span:.2f}× < {args.fit_min_span:.0f}× "
                            f"或 fit 点 {len(fit_pts)}<{MIN_POINTS}:只校 base)")
        else:
            if n_ok:
                suggest = p90 * args.safety
            else:
                suggest = p90 * args.safety
                note.append(f"censored-only({n_cen} oom rows), lower-bound suggestion")
            if n_cen:
                floor = max(cen) * args.safety
                if floor > suggest:
                    suggest = floor
                    note.append(f"censored floor applied({n_cen} oom rows)")
            if fit_ready and not sc:
                old_a = float(prev) if prev_num else 0.0
                guard_xs = list(distinct)
                A, B, keep = blend_curve(A_raw, B_raw, old_a, 0.0, guard_xs)
                if prev_num and keep < 1.0:
                    note.append(f"scale fit damped to {keep:.0%} keep(曲线护栏 ±30%)")
                unmet = [(x, f) for x, f in floor_pts if A + B * x < f - 1]
                if unmet:
                    note.append(f"censored floor 未达({len(unmet)} 点,±30% 限幅内不可及,"
                                f"连续同向后放宽)")
                ref_new = _median(xs)
                base_new = max(A + B * ref_new, 1.0)
                base_i = int(round(base_new))
                frac_new = min(max(A / base_new, 0.0), 1.0)
                suggest = base_new
                metric_label = metric if metric in ("bam_mb", "fastq_mb", "sample_count") \
                    else "bam_mb"
                scale_upd = {
                    "size_metric": metric_label,
                    "ref_input_size": int(round(ref_new)),
                    "fixed_frac": round(frac_new, 3),
                    "min_mem_mb": int(round(base_i * SCALE_NEW_MIN_MEM_FRAC)),
                    "max_mem_mb": int(round(base_i * SCALE_NEW_MAX_MEM_FRAC)),
                }
                create_scale = True
                scale_rows["suggest_ref_input_size_mb"] = f"{ref_new:.0f}"
                scale_rows["suggest_fixed_frac"] = f"{frac_new:.3f}"
                if prev_num:
                    scale_rows["curve_delta_max_pct"] = curve_delta(A, B, guard_xs, old_a, 0.0)
                note.append("scale 未配,拟合达标(--apply 新建 scale;"
                            "min/max 暂定 0.5×/4× base,请人工审定)")
                damped_value = False
            elif pts_m:
                note.append("scale 未配;拟合不足(跨度/点数门槛),仅静态 mem")
        if suggest is None:
            suggest = float(prev) if prev_num else 0.0

        if prev_num:
            if low_n:
                note.append(f"insufficient points({n_ok + n_cen}<{MIN_POINTS}), suggestion informational")
            if damped_value:
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
        io_obs = hist[rule]["io"]
        io_n = len(io_obs)
        io_p90_raw = percentile(io_obs, 0.90) if io_obs else None
        io_suggest = None
        io_delta = None
        if io_p90_raw is not None:
            io_suggest = float(io_p90_raw)
            prev_io = entry.get("io_bandwidth_mb_s")
            if isinstance(prev_io, (int, float)) and not isinstance(prev_io, bool) \
                    and prev_io > 0:
                io_bounded = min(max(io_suggest, prev_io * (1 - DAMPING)),
                                 prev_io * (1 + DAMPING))
                if abs(io_bounded - io_suggest) > 1e-6:
                    note.append("io damped ±30%")
                    io_suggest = float(io_bounded)
                io_delta = 100.0 * (io_suggest - prev_io) / prev_io
            if io_n < MIN_POINTS:
                note.append(f"insufficient io points({io_n}<{MIN_POINTS}), informational")
        io_pairs = [(p[0], p[4]) for p in fit_pts if p[4] is not None]
        io_fit = fit_affine([p[0] for p in io_pairs], [p[1] for p in io_pairs])
        row = {
            "rule": rule, "n_ok": n_ok, "n_censored": n_cen,
            "observed_p90_mb": f"{p90:.0f}",
            "suggest_mb": f"{suggest:.0f}",
            "prev_declared_mb": prev if prev else "",
            "delta_pct": f"{delta:+.0f}" if delta is not None else "",
        }
        row.update(scale_rows)
        row.update({
            "io_n": io_n,
            "observed_io_p90_mb_s": f"{io_p90_raw:.0f}" if io_p90_raw is not None else "",
            "suggest_io_bandwidth_mb_s": f"{io_suggest:.0f}" if io_suggest is not None else "",
            "io_delta_pct": f"{io_delta:+.0f}" if io_delta is not None else "",
            "io_fit_a_mb_s": f"{io_fit[0]:.2f}" if io_fit else "",
            "io_fit_b_mb_s_per_mb": f"{io_fit[1]:.5f}" if io_fit else "",
            "note": "; ".join(note),
        })
        out_rows.append(row)
        if args.apply:
            new_mb = int(round(suggest))
            upd = {}
            if prev_num:
                if n_ok + n_cen >= MIN_POINTS and new_mb != int(prev):
                    upd["mem_mb"] = new_mb
            elif rule not in family_rules:
                upd["mem_mb"] = new_mb
            if scale_upd and rule not in family_rules \
                    and (rule in declared or "mem_mb" in upd):
                upd["scale"] = scale_upd
                upd["create_scale"] = create_scale
            if upd:
                applies[rule] = upd

    if out_rows:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(out_rows)
        print(f"suggest: {len(out_rows)} rules -> {args.out}")
    else:
        print("no rules with usable history")

    if args.apply:
        apply_to_config(args.config, applies)
        for rule in sorted(applies):
            upd = applies[rule]
            old = declared.get(rule, {}).get("mem_mb")
            if upd.get("mem_mb") is not None:
                print(f"apply: {rule} {old if old else '(new)'} -> {upd['mem_mb']}")
            if upd.get("scale"):
                kind = "new" if upd.get("create_scale") else "update"
                print(f"apply: {rule} scale({kind}) -> {upd['scale']}")
        print(f"apply: {len(applies)} rules written -> {args.config}"
              f"(n<{MIN_POINTS} 未动 {sum(1 for r in out_rows if 'insufficient' in r['note'])} 条)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
