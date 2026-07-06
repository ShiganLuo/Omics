from snakemake.logging import logger
import time
import os
import sys

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule mimseq_coverage:
    """Calculate coverage and generate plots."""
    input:
        mods_done = outdir + "/mods.done",
    output:
        coverage_done = touch(outdir + "/coverage.done"),
    log:
        logdir + "/mimseq_coverage.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "coverage", "run.py"),
        control_cond = config.get("Params", {}).get("mimseq", {}).get("control_cond", ""),
        mito_trnas = config.get("genome", {}).get("mito_trnas", ""),
        misinc_thresh = config.get("Params", {}).get("mimseq", {}).get("misinc_thresh", 0.1),
    threads: 4
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"mimseq_coverage_{current_time}.sh")

            cmd = [
                "/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin/python", params.script,
                "--out", outdir + "/",
                "--misinc-thresh", str(params.misinc_thresh),
            ]
            if params.control_cond:
                cmd += ["--control-cond", params.control_cond]
            if params.mito_trnas:
                cmd += ["--mito-trnas", params.mito_trnas]
            if not config.get("Params", {}).get("mimseq", {}).get("no_cca", False):
                cmd += ["--cca"]
            if config.get("Params", {}).get("mimseq", {}).get("double_cca", False):
                cmd += ["--double-cca"]

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("export PATH=/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin:$PATH\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"mimseq_coverage failed: {e}\n")
            raise

rule mimseq_coverage_result:
    """Result rule for coverage module."""
    input:
        coverage_done = outdir + "/coverage.done",
