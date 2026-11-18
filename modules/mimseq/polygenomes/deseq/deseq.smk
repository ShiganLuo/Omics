include: "../../../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# 图片格式:转发给 wrapper 的 -f/--format,经 MIMSEQ_PLOT_FORMATS 交给
# vendored 绘图代码(ggsave/savefig)。figures 是 touch() 规则的副产物,不进
# rule output(② 展开不适用);默认 ["pdf"] 保持上游行为。
# ① parse 期读入 + 规范化（str 收成单元素 list；合法值由 schema 声明）
# 参数精确到 rule:config 全量路径读取 Params.mimseq.mimseq_deseq.formats
# (无 module 级共享回退)。
_plot_formats_mimseq_deseq = config.get("Params", {}).get("mimseq", {}).get("mimseq_deseq", {}).get("formats") or ["pdf"]
if isinstance(_plot_formats_mimseq_deseq, str):
    _plot_formats_mimseq_deseq = [_plot_formats_mimseq_deseq]

rule mimseq_deseq:
    """Run DESeq2 differential expression analysis."""
    input:
        coverage_done = outdir + "/{genome}/coverage.done",
    output:
        deseq_done = touch(outdir + "/{genome}/deseq.done"),
    log:
        logdir + "/{genome}/mimseq_deseq.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "deseq", "run.py"),
        control_cond = config.get("Params", {}).get("mimseq", {}).get("mimseq_deseq", {}).get("control_cond", ""),
        p_adj = config.get("Params", {}).get("mimseq", {}).get("mimseq_deseq", {}).get("p_adj", 0.05),
        formats = _plot_formats_mimseq_deseq,   # ③ 转发
        mito_trnas = lambda wildcards: config["genome"][wildcards.genome].get("mito_trnas", ""),
    threads: rule_res("mimseq_deseq").get("threads", 4)
    benchmark: logdir + "/{genome}/mimseq_deseq.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_deseq", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "../../mimseq.yaml"
    container:
        sif("../../mimseq.yaml", "mimseq_deseq")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(RECEIPT_DIR, f"mimseq_deseq_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, wildcards.genome, f"mimseq_deseq.sh")

            cmd = [
                "python", params.script,
                "--out", outdir + f"/{wildcards.genome}/",
                "--p-adj", str(params.p_adj),
            ]
            if params.control_cond:
                cmd += ["--control-cond", params.control_cond]
            if params.mito_trnas:
                cmd += ["--mito-trnas", params.mito_trnas]
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

rule mimseq_deseq_result:
    """Result rule for deseq module."""
    input:
        deseq_done = outdir + "/{genome}/deseq.done",
