include: "../../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# 图片格式:转发给 wrapper 的 -f/--format,经 MIMSEQ_PLOT_FORMATS 交给
# vendored 绘图代码(ggsave/savefig)。figures 是 touch() 规则的副产物,不进
# rule output(② 展开不适用);默认 ["pdf"] 保持上游行为。
# ① parse 期读入 + 规范化（str 收成单元素 list；合法值由 schema 声明）
# 参数精确到 rule:config 全量路径读取 Params.mimseq.mimseq_align.formats
# (无 module 级共享回退)。
_plot_formats_mimseq_align = config.get("Params", {}).get("mimseq", {}).get("mimseq_align", {}).get("formats") or ["pdf"]
if isinstance(_plot_formats_mimseq_align, str):
    _plot_formats_mimseq_align = [_plot_formats_mimseq_align]

rule mimseq_align:
    """Align reads using GSNAP."""
    input:
        sample_data = outdir + "/sample_data.tsv",
        state_dir = outdir + "/state",
    output:
        align_done = touch(outdir + "/align.done"),
    log:
        logdir + "/mimseq_align.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "align", "run.py"),
        name = config.get("Params", {}).get("mimseq", {}).get("mimseq_align", {}).get("name", "tRNAseq"),
        mismatches = config.get("Params", {}).get("mimseq", {}).get("mimseq_align", {}).get("max_mismatches", 0.075),
        remap = config.get("Params", {}).get("mimseq", {}).get("mimseq_align", {}).get("remap", True),
        keep_temp = config.get("Params", {}).get("mimseq", {}).get("mimseq_align", {}).get("keep_temp", False),
        formats = _plot_formats_mimseq_align,   # ③ 转发
    threads: rule_res("mimseq_align").get("threads", 10)
    benchmark: logdir + "/mimseq_align.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_align", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../mimseq.yaml"
    container:
        sif("../mimseq.yaml", "mimseq_align")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(RECEIPT_DIR, f"mimseq_align_{current_time}.sh")
            plain_script = os.path.join(outdir, f"mimseq_align.sh")

            cmd = [
                "python", params.script,
                "--sample-data", input.sample_data,
                "--name", params.name,
                "--out", outdir + "/",
                "--threads", str(threads),
                "--mismatches", str(params.mismatches),
            ]
            if params.remap:
                cmd += ["--remap"]
            if params.keep_temp:
                cmd += ["--keep-temp"]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]

            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule mimseq_align_result:
    """Result rule for align module."""
    input:
        align_done = outdir + "/align.done",
