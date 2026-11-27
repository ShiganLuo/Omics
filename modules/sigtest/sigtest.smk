include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")

sigtest_params = config.get("Params", {}).get("sigtest", {})
sigtest_contrasts = [str(c) for c in (config.get("contrasts") or [])]

# Required: the long table to test. Raise at DAG-build time (input
# function) when the key is missing, naming the JSON key to fill.
def _sigtest_table(wildcards):
    """Resolve the required input table path or fail naming the config key."""
    return _require(
        config.get("sigtest_table"),
        "Params.sigtest.table",
        tool="sigtest",
    )


_plot_formats = (sigtest_params.get("plot", {}) or {}).get("formats")
if _plot_formats is None:
    _plot_formats = ["png"]
if isinstance(_plot_formats, str):
    _plot_formats = [_plot_formats]


rule sigtest_run:
    input:
        table = _sigtest_table,
    output:
        tsv = outdir + "/{contrast}/significance_test.tsv",
        plots = ([outdir + "/{contrast}/significance_bars.done"]
                 if _plot_formats else []),
    log:
        logdir + "/sigtest/{contrast}.log"
    threads: rule_res("sigtest_run").get("threads", 1)
    benchmark: logdir + "/sigtest/{contrast}.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sigtest_run", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 4000),
    params:
        sigtest_script = os.path.join(ROOT_DIR, "modules", "sigtest", "bin", "sigtest.py"),
        p = sigtest_params,
        formats = _plot_formats,
    conda:
        "sigtest.yaml"
    container:
        sif("sigtest.yaml", "sigtest_run")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sigtest_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sigtest for contrast {wildcards.contrast} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"sigtest_run_{wildcards.contrast}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, "sigtest_run.sh")
            cmd = [
                "python", params.sigtest_script,
                "-i", str(input.table),
                "-o", sample_outdir,
                "-f", (params.p.get("feature_col_by_contrast", {})
                       .get(str(wildcards.contrast))
                       or params.p.get("feature_col") or "feature"),
                "--group-col", params.p.get("group_col") or "group",
                "-v", params.p.get("value_col") or "value",
                "-m", params.p.get("method") or "auto",
                "-c", params.p.get("correction") or "none",
                "-a", str(0.05 if params.p.get("alpha") is None else params.p.get("alpha")),
                "-y", params.p.get("ylabel") or "value",
            ]
            groups = (params.p.get("groups_by_contrast", {}).get(str(wildcards.contrast))
                      or params.p.get("groups") or [])
            if len(groups) == 2:
                cmd += ["-g", str(groups[0]), str(groups[1])]
            if params.p.get("paired"):
                cmd += ["--paired"]
            if params.p.get("annotate_ns"):
                cmd += ["--annotate-ns"]
            if params.p.get("no_errorbars"):
                cmd += ["--no-errorbars"]
            for fmt in params.formats:
                cmd += ["--formats", str(fmt)]
            cmd += ["--threads", str(threads)]
            legend = params.p.get("legend") or []
            if len(legend) == 2:
                cmd += ["--legend", str(legend[0]), str(legend[1])]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write('echo "sigtest completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(str(x) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"sigtest completed for contrast {wildcards.contrast}")
        except Exception as e:
            rule_logger.error(f"Error occurred during sigtest for {wildcards.contrast}: {e}\n")
            raise RuntimeError(f"Error occurred during sigtest for {wildcards.contrast}: {e}, you can check the log file {log_path} for more details.\n")


rule sigtest_result:
    input:
        tsv = expand(outdir + "/{contrast}/significance_test.tsv", contrast=sigtest_contrasts),
        plots = expand(outdir + "/{contrast}/significance_bars.done", contrast=sigtest_contrasts),
