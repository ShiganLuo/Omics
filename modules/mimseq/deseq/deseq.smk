include: "../../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

_plot_formats_mimseq_deseq = config.get("Params", {}).get("mimseq", {}).get("mimseq_deseq", {}).get("formats") or ["pdf"]
if isinstance(_plot_formats_mimseq_deseq, str):
    _plot_formats_mimseq_deseq = [_plot_formats_mimseq_deseq]

rule mimseq_deseq:
    """Run DESeq2 differential expression analysis."""
    input:
        coverage_done = outdir + "/coverage.done",
    output:
        deseq_done = touch(outdir + "/deseq.done"),
    log:
        logdir + "/mimseq_deseq.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "deseq", "run.py"),
        control_cond = config.get("Params", {}).get("mimseq", {}).get("mimseq_deseq", {}).get("control_cond", ""),
        p_adj = config.get("Params", {}).get("mimseq", {}).get("mimseq_deseq", {}).get("p_adj", 0.05),
        formats = _plot_formats_mimseq_deseq,
        mito_trnas = config.get("genome", {}).get("mito_trnas", ""),
    threads: rule_res("mimseq_deseq").get("threads", 4)
    benchmark: logdir + "/mimseq_deseq.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_deseq", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "../mimseq.yaml"
    container:
        sif("../mimseq.yaml", "mimseq_deseq")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = receipt_script(f"mimseq_deseq_{current_time}.sh")
            plain_script = os.path.join(outdir, f"mimseq_deseq.sh")

            cmd = [
                "python", params.script,
                "--out", outdir + "/",
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
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule mimseq_deseq_result:
    """Result rule for deseq module."""
    input:
        deseq_done = outdir + "/deseq.done",
