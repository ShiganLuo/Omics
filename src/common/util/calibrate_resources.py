#!/usr/bin/env python
"""benchmark 校准闭环脚本(对应 docs/resource-balancing-proposal.md §6):benchmark tsv → resource_history.tsv → resources_suggest.tsv →(--apply) config.

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
     传同一共享路径);
  规则名:benchmark 文件名即规则名;模块文件名与 config 资源键不一致的按
     RULE_ALIASES/SUFFIX_ALIASES 固定映射;PE/SE 规则共用同一 benchmark
     文件名的(bowtie2_align/trimming)按 raw.json 样本布局逐行归属,
     归属不了的行丢弃并计数(不猜);
  6.3 建议:按 rule 聚合历史,取 max_rss 的 p90 × 安全系数(默认 1.3)为建议值;
     ≥3 个数据点才建议;单次调整幅度限 ±30%(防震荡);OOM 行按删失处理
     ——只用其 max_rss 作下界,不压低建议;
  输出 --out(默认 <root>/resources_suggest.tsv):rule、n、observed_p90_mb、
     suggest_mb、prev_declared_mb、delta_pct、censored_n、note,外加 io 列
     io_n、observed_io_p90_mb_s、suggest_io_bandwidth_mb_s、io_delta_pct
     (设计 §10.1:io 权重 = p90((io_in+io_out)/s),±30% 限幅,只出建议人工合入);
  --io-probe(默认关):实测 sustainable_disk_bandwidth_mb_s(设计 §10.2),
     优先现成工具 fio/hdparm,缺失时内置 256MB 测速;结果写 <probe-dir>/io_probe.json。
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
import shutil
import subprocess
import sys
import time

DEFAULT_SAFETY = 1.3
DAMPING = 0.30
MIN_POINTS = 3

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


def read_benchmarks(root):
    """Collect (rule, attempt, seconds, max_rss_mb, io_mb_s, source) from benchmark outputs.

    Accepts both layouts: benchmarks/<rule>.tsv and the Omics log layout
    <...>/<rule>.benchmark.tsv(.txt) (snakemake ``benchmark:`` output).
    io_mb_s = (io_in + io_out) / wall seconds (设计 §8.3); None when the tsv
    has no io columns (old benchmarks) so zeros never fake "no I/O".
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
        oom_marker = os.path.join(os.path.dirname(tsv), base + ".OOM")
        censored_all = os.path.exists(oom_marker)
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
                    rows.append((rule, 1 if not censored_all else 2, secs, rss_mb, io_rate, tsv))
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


def _median(vals):
    xs = sorted(vals)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


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

    rows = []
    for r in roots:
        got = read_benchmarks(r)
        if not got:
            print(f"warning: no benchmark tsv under {r}", file=sys.stderr)
        rows += got
    if not rows:
        print("no benchmark tsv under any --root(先跑一次带 benchmark 的运行)", file=sys.stderr)
        return 1

    hist_header = ["run_id", "ts", "rule", "attempt", "seconds", "max_rss_mb",
                   "io_mb_s", "censored", "source"]
    new_history = 0
    if not args.no_append:
        exists = os.path.exists(args.history)
        seen = set()
        if exists:
            with open(args.history, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                old_rows = list(reader)
                old_header = reader.fieldnames or []
            if "io_mb_s" not in old_header:
                # 历史表 schema 升级(设计 §10.1):补 io_mb_s 列,旧行留空,内容与行序保留
                with open(args.history, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=hist_header, delimiter="\t")
                    w.writeheader()
                    for rec in old_rows:
                        rec.setdefault("io_mb_s", "")
                        w.writerow({k: rec.get(k, "") for k in hist_header})
            for rec in old_rows:
                seen.add((rec["run_id"], rec["rule"], rec["attempt"], rec["max_rss_mb"]))
        run_id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        with open(args.history, "a", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            if not exists:
                w.writerow(hist_header)
            ts = datetime.datetime.now().isoformat(timespec="seconds")
            for rule, attempt, secs, rss_mb, io_rate, src in rows:
                key = (run_id, rule, attempt, rss_mb)
                if sum(1 for r in rows if r[0] == rule and r[1] == attempt and r[3] == rss_mb) > 1:
                    key = (run_id, rule, attempt, rss_mb, src)
                if key in seen:
                    continue
                w.writerow([run_id, ts, rule, attempt,
                            f"{secs:.1f}" if secs else "", f"{rss_mb:.0f}",
                            f"{io_rate:.1f}" if io_rate is not None else "",
                            "1" if attempt > 1 else "0", src])
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
                hist.setdefault(rule, {"all": [], "censored": [], "io": []})
                (hist[rule]["censored"] if censored else hist[rule]["all"]).append(rss)
                raw_io = rec.get("io_mb_s")
                if raw_io not in (None, "") and not censored:
                    try:
                        hist[rule]["io"].append(float(raw_io))
                    except ValueError:
                        pass

    declared = {}
    if args.config and os.path.exists(args.config):
        with open(args.config) as fh:
            declared = (json.load(fh).get("resources") or {})

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
            p90 = max(cen)
            suggest = p90 * args.safety
            note.append(f"censored-only({n_cen} oom rows), lower-bound suggestion")
        if n_cen:
            floor = max(cen) * args.safety
            if floor > suggest:
                suggest = floor
                note.append(f"censored floor applied({n_cen} oom rows)")
        prev = declared.get(rule, {}).get("mem_mb")
        if isinstance(prev, (int, float)):
            if low_n:
                note.append(f"insufficient points({n_ok + n_cen}<{MIN_POINTS}), suggestion informational")
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
        # io 权重建议(设计 §10.1):suggest = p90((io_in+io_out)/s),无安全系数;
        # ±30% 限幅对 prev 声明值;<3 点仅信息性;只出建议,--apply 不写 io。
        io_obs = hist[rule]["io"]
        io_n = len(io_obs)
        io_p90_raw = percentile(io_obs, 0.90) if io_obs else None
        io_suggest = None
        io_delta = None
        if io_p90_raw is not None:
            io_suggest = float(io_p90_raw)
            prev_io = declared.get(rule, {}).get("io_bandwidth_mb_s")
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
        out_rows.append({
            "rule": rule, "n_ok": n_ok, "n_censored": n_cen,
            "observed_p90_mb": f"{p90:.0f}",
            "suggest_mb": f"{suggest:.0f}",
            "prev_declared_mb": prev if prev else "",
            "delta_pct": f"{delta:+.0f}" if delta is not None else "",
            "io_n": io_n,
            "observed_io_p90_mb_s": f"{io_p90_raw:.0f}" if io_p90_raw is not None else "",
            "suggest_io_bandwidth_mb_s": f"{io_suggest:.0f}" if io_suggest is not None else "",
            "io_delta_pct": f"{io_delta:+.0f}" if io_delta is not None else "",
            "note": "; ".join(note),
        })
        if args.apply:
            new_mb = int(round(suggest))
            if isinstance(prev, (int, float)):
                if n_ok + n_cen >= MIN_POINTS and new_mb != int(prev):
                    applies[rule] = new_mb
            elif rule not in family_rules:
                applies[rule] = new_mb

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
            old = declared.get(rule, {}).get("mem_mb")
            print(f"apply: {rule} {old if old else '(new)'} -> {applies[rule]}")
        print(f"apply: {len(applies)} rules written -> {args.config}"
              f"(n<{MIN_POINTS} 未动 {sum(1 for r in out_rows if 'insufficient' in r['note'])} 条)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
