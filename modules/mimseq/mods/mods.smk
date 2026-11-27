include: "../../common/common.smk"


ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

_plot_formats_mimseq_mods = config.get("Params", {}).get("mimseq", {}).get("mimseq_mods", {}).get("formats") or ["pdf"]
if isinstance(_plot_formats_mimseq_mods, str):
    _plot_formats_mimseq_mods = [_plot_formats_mimseq_mods]

rule mimseq_mods:
    """Quantify modifications and misincorporation."""
    input:
        clusters_done = outdir + "/clusters.done",
    output:
        mods_done = touch(outdir + "/mods.done"),
    log:
        logdir + "/mimseq_mods.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "mods", "run.py"),
        name = config.get("Params", {}).get("mimseq", {}).get("mimseq_mods", {}).get("name", "tRNAseq"),
        min_cov = config.get("Params", {}).get("mimseq", {}).get("mimseq_mods", {}).get("min_cov", 0.0005),
        misinc_thresh = config.get("Params", {}).get("mimseq", {}).get("mimseq_mods", {}).get("misinc_thresh", 0.1),
        formats = _plot_formats_mimseq_mods,
    threads: rule_res("mimseq_mods").get("threads", 4)
    benchmark: logdir + "/mimseq_mods.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_mods", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../mimseq.yaml"
    container:
        sif("../mimseq.yaml", "mimseq_mods")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = receipt_script(f"mimseq_mods_{current_time}.sh")
            plain_script = os.path.join(outdir, f"mimseq_mods.sh")

            cmd = [
                "python", params.script,
                "--out", outdir + "/",
                "--name", params.name,
                "--threads", str(threads),
                "--min-cov", str(params.min_cov),
                "--misinc-thresh", str(params.misinc_thresh),
            ]
            if config.get("Params", {}).get("mimseq", {}).get("mimseq_mods", {}).get("no_cca", False):
                cmd += ["--cca"]
            if config.get("Params", {}).get("mimseq", {}).get("mimseq_mods", {}).get("remap", True):
                cmd += ["--remap"]
            if config.get("Params", {}).get("mimseq", {}).get("mimseq_mods", {}).get("crosstalks", False):
                cmd += ["--crosstalks"]
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

rule mimseq_mods_result:
    """Result rule for mods module."""
    input:
        mods_done = outdir + "/mods.done",
