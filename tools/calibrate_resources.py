#!/usr/bin/env python
"""方案 D 校准脚本:benchmark tsv → resource_history.tsv → resources_suggest.tsv.

用法(在某个 workflow 输出目录或其父目录下,或用 --glob 指定):
    python tools/calibrate_resources.py [--root OUTDIR] [--history PATH] [--out PATH]

行为(对应 docs/resource-balancing-proposal.md §5.4):
  D.1 收集:扫描 --root 下所有 benchmarks/*.tsv(snakemake `benchmark:` 输出,
     列 s / h:m:s / max_rss / max_vms,每 job 每 attempt 追加多行),连同最近的
     run_summary.json(规模/预算上下文),**追加**进 --history(默认
     resource_history.tsv,append-only,跨运行累积);
  D.3 建议:按 rule 聚合历史,取 max_rss 的 p90 × 安全系数(默认 1.3)为建议值;
     ≥3 个数据点才建议;单次调整幅度限 ±30%(防震荡);OOM 行(attempt>1 或
     退出码 137 上下文)按删失处理——只用其 max_rss 作下界,不压低建议;
  输出 --out(默认 resources_suggest.tsv):rule、n、observed_p90_mb、suggest_mb、
     prev_declared_mb、delta_pct、censored_n、note。**人工审核后**合入 config,
     本脚本绝不改写 config。
"""
import argparse
import csv
import datetime
import glob
import json
import os
import sys

DEFAULT_SAFETY = 1.3
DAMPING = 0.30          # 单次调整限幅 ±30%(D.4)
MIN_POINTS = 3          # ≥3 个数据点才动形状(D.3)


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


def read_benchmarks(root):
    """Yield (rule, attempt, seconds, max_rss_kb, source) from benchmarks/*.tsv."""
    rows = []
    for tsv in glob.glob(os.path.join(root, "**", "benchmarks", "*.tsv"), recursive=True):
        rule = os.path.basename(tsv)[:-4]
        # OOM 删失标记:若某次运行检测到 OOM,wrapper 会写 benchmarks/<rule>.OOM;
        # 该规则本次全部观测按删失(真实需求 ≥ 观测)处理。
        oom_marker = os.path.join(os.path.dirname(tsv), rule + ".OOM")
        censored_all = os.path.exists(oom_marker)
        try:
            with open(tsv, newline="") as fh:
                for i, rec in enumerate(csv.DictReader(fh, delimiter="\t")):
                    # benchmark tsv 每 job(含每次重试)追加一行;attempt 无法从行号
                    # 唯一还原(多 wildcard job 共享文件),删失判定只用 .OOM 标记。
                    attempt = 1
                    secs = parse_hms(rec.get("s") or rec.get("h:m:s"))
                    try:
                        rss_kb = float(rec.get("max_rss", 0) or 0)
                    except ValueError:
                        rss_kb = 0.0
                    rows.append((rule, 1 if not censored_all else 2, secs, rss_kb, tsv))
        except OSError:
            continue
    return rows


def percentile(vals, q):
    if not vals:
        return None
    xs = sorted(vals)
    k = max(0, min(len(xs) - 1, int(round(q * (len(xs) - 1)))))
    return xs[k]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="workflow 输出根目录(递归找 benchmarks/*.tsv)")
    ap.add_argument("--history", default="resource_history.tsv", help="累计历史表(append-only)")
    ap.add_argument("--out", default="resources_suggest.tsv", help="建议表输出")
    ap.add_argument("--safety", type=float, default=DEFAULT_SAFETY, help="安全系数(默认 1.3)")
    ap.add_argument("--config", default=None, help="可选:workflow config json,读当前声明值列")
    ap.add_argument("--no-append", action="store_true", help="只出建议,不追加历史")
    args = ap.parse_args(argv)

    rows = read_benchmarks(args.root)
    if not rows:
        print(f"no benchmark tsv under {args.root}(先跑一次带 benchmark 的运行)", file=sys.stderr)
        return 1

    # ---- D.1 append to history ----
    new_history = 0
    if not args.no_append:
        exists = os.path.exists(args.history)
        seen = set()
        if exists:
            with open(args.history, newline="") as fh:
                for rec in csv.DictReader(fh, delimiter="\t"):
                    seen.add((rec["run_id"], rec["rule"], rec["attempt"], rec["max_rss_kb"]))
        run_id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        with open(args.history, "a", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            if not exists:
                w.writerow(["run_id", "ts", "rule", "attempt", "seconds", "max_rss_kb", "censored", "source"])
            ts = datetime.datetime.now().isoformat(timespec="seconds")
            for rule, attempt, secs, rss_kb, src in rows:
                key = (run_id, rule, attempt, rss_kb)
                # per-run dedupe: same (rule,attempt,rss) repeated across wildcard jobs
                if sum(1 for r in rows if r[0] == rule and r[1] == attempt and r[3] == rss_kb) > 1:
                    key = (run_id, rule, attempt, rss_kb, src)
                if key in seen:
                    continue
                w.writerow([run_id, ts, rule, attempt,
                            f"{secs:.1f}" if secs else "", f"{rss_kb:.0f}",
                            "1" if attempt > 1 else "0", src])
                new_history += 1
        print(f"history: appended {new_history} rows -> {args.history}")

    # ---- load history for suggestions ----
    hist = {}
    if os.path.exists(args.history):
        with open(args.history, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                try:
                    rss = float(rec["max_rss_kb"])
                except (ValueError, KeyError):
                    continue
                censored = (rec.get("censored") == "1")
                rule = rec["rule"]
                hist.setdefault(rule, {"all": [], "censored": []})
                (hist[rule]["censored"] if censored else hist[rule]["all"]).append(rss)

    declared = {}
    if args.config and os.path.exists(args.config):
        with open(args.config) as fh:
            declared = (json.load(fh).get("resources") or {})

    # ---- D.3 suggest ----
    out_rows = []
    for rule in sorted(hist):
        obs = hist[rule]["all"]
        cen = hist[rule]["censored"]
        n_ok = len(obs)
        n_cen = len(cen)
        if n_ok == 0 and n_cen == 0:
            continue
        note = []
        if n_ok:
            p90 = percentile(obs, 0.90)
            suggest = p90 * args.safety
            if n_ok + n_cen < MIN_POINTS:
                note.append(f"insufficient points({n_ok + n_cen}<{MIN_POINTS}), suggestion informational")
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
            damped = prev * (1 + DAMPING) if suggest > prev else prev * (1 - DAMPING)
            # 限幅:单次调整 ≤±30%(D.4)
            bounded = min(max(suggest, prev * (1 - DAMPING)), prev * (1 + DAMPING))
            if abs(bounded - suggest) > 1:
                note.append("damped ±30%")
                suggest = bounded
            delta = 100.0 * (suggest - prev) / prev
        else:
            delta = None
        out_rows.append({
            "rule": rule, "n_ok": n_ok, "n_censored": n_cen,
            "observed_p90_kb": f"{p90:.0f}",
            "suggest_mb": f"{suggest / 1024:.0f}",
            "prev_declared_mb": prev if prev else "",
            "delta_pct": f"{delta:+.0f}" if delta is not None else "",
            "note": "; ".join(note),
        })

    if out_rows:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(out_rows)
        print(f"suggest: {len(out_rows)} rules -> {args.out}(人工审核后合入 config)")
    else:
        print("no rules with usable history")
    return 0


if __name__ == "__main__":
    sys.exit(main())
