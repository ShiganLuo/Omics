include: "../common/common.smk"

# ---------------------------------------------------------------------------
# sigtest: per-feature two-group significance testing with star labels
# and bracket-annotated bar plots.
#
# Consumes one long table (feature x group x value) via config key
# "sigtest_table" (a file path). Group-level module: no {sample_id}
# wildcard, log goes to logdir directly (subworkflow passes
# {logdir}/group for group-level modules).
# ---------------------------------------------------------------------------
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")

sigtest_params = config.get("Params", {}).get("sigtest", {})
# Contrast names come from the subworkflow (derived from its group_pairs);
# the panel rule expands over them at parse time (scanpy_result pattern).
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


rule sigtest_run:
    input:
        table = _sigtest_table,
    output:
        tsv = outdir + "/{contrast}/significance_test.tsv",
        png = outdir + "/{contrast}/significance_bars.png",
    log:
        logdir + "/sigtest/{contrast}.log"
    threads: rule_res("sigtest_run").get("threads", 1)
    resources:
        mem_mb = rule_res("sigtest_run").get("mem_mb", 4000),
    params:
        sigtest_script = os.path.join(ROOT_DIR, "modules", "sigtest", "bin", "sigtest.py"),
        p = sigtest_params,
    conda:
        "sigtest.yaml"
    container:
        sif("sigtest.yaml", "sigtest_run")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sigtest_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sigtest for contrast {wildcards.contrast} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"sigtest_run_{wildcards.contrast}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "sigtest_run.sh")
            cmd = [
                "python", params.sigtest_script,
                "-i", str(input.table),
                "-o", sample_outdir,
                "-f", params.p.get("feature_col", "feature"),
                "--group-col", params.p.get("group_col", "group"),
                "-v", params.p.get("value_col", "value"),
                "-m", params.p.get("method", "auto"),
                "-c", params.p.get("correction", "none"),
                "-a", str(params.p.get("alpha", 0.05)),
                "-y", params.p.get("ylabel", "value"),
            ]
            groups = params.p.get("groups") or []
            if len(groups) == 2:
                cmd += ["-g", str(groups[0]), str(groups[1])]
            if params.p.get("paired"):
                cmd += ["--paired"]
            if params.p.get("annotate_ns"):
                cmd += ["--annotate-ns"]
            if params.p.get("no_errorbars"):
                cmd += ["--no-errorbars"]
            if params.p.get("no_plot"):
                cmd += ["--no-plot"]
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
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"sigtest completed for contrast {wildcards.contrast}")
        except Exception as e:
            rule_logger.error(f"Error occurred during sigtest for {wildcards.contrast}: {e}\n")
            raise RuntimeError(f"Error occurred during sigtest for {wildcards.contrast}: {e}, you can check the log file {log_path} for more details.\n")


rule sigtest_result:
    input:
        tsv = expand(outdir + "/{contrast}/significance_test.tsv", contrast=sigtest_contrasts),
        png = expand(outdir + "/{contrast}/significance_bars.png", contrast=sigtest_contrasts),
