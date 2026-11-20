include: "../../../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# 图片格式:转发给 wrapper 的 -f/--format,经 MIMSEQ_PLOT_FORMATS 交给
# vendored 绘图代码(ggsave/savefig)。figures 是 touch() 规则的副产物,不进
# rule output(② 展开不适用);默认 ["pdf"] 保持上游行为。
# ① parse 期读入 + 规范化（str 收成单元素 list；合法值由 schema 声明）
# 参数精确到 rule:config 全量路径读取 Params.mimseq.mimseq_coverage.formats
# (无 module 级共享回退)。
_plot_formats_mimseq_coverage = config.get("Params", {}).get("mimseq", {}).get("mimseq_coverage", {}).get("formats") or ["pdf"]
if isinstance(_plot_formats_mimseq_coverage, str):
    _plot_formats_mimseq_coverage = [_plot_formats_mimseq_coverage]

rule mimseq_coverage:
    """Calculate coverage and generate plots."""
    input:
        mods_done = outdir + "/{genome}/mods.done",
    output:
        coverage_done = touch(outdir + "/{genome}/coverage.done"),
    log:
        logdir + "/{genome}/mimseq_coverage.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "coverage", "run.py"),
        control_cond = config.get("Params", {}).get("mimseq", {}).get("mimseq_coverage", {}).get("control_cond", ""),
        mito_trnas = lambda wildcards: config["genome"][wildcards.genome].get("mito_trnas", ""),
        misinc_thresh = config.get("Params", {}).get("mimseq", {}).get("mimseq_coverage", {}).get("misinc_thresh", 0.1),
        formats = _plot_formats_mimseq_coverage,   # ③ 转发
    threads: rule_res("mimseq_coverage").get("threads", 4)
    benchmark: logdir + "/{genome}/mimseq_coverage.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_coverage", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "../../mimseq.yaml"
    container:
        sif("../../mimseq.yaml", "mimseq_coverage")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(RECEIPT_DIR, f"mimseq_coverage_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, wildcards.genome, f"mimseq_coverage.sh")

            cmd = [
                "python", params.script,
                "--out", outdir + f"/{wildcards.genome}/",
                "--misinc-thresh", str(params.misinc_thresh),
            ]
            if params.control_cond:
                cmd += ["--control-cond", params.control_cond]
            if params.mito_trnas:
                cmd += ["--mito-trnas", params.mito_trnas]
            if not config.get("Params", {}).get("mimseq", {}).get("mimseq_coverage", {}).get("no_cca", False):
                cmd += ["--cca"]
            if config.get("Params", {}).get("mimseq", {}).get("mimseq_coverage", {}).get("double_cca", False):
                cmd += ["--double-cca"]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]

            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule mimseq_coverage_result:
    """Result rule for coverage module."""
    input:
        coverage_done = outdir + "/{genome}/coverage.done",
