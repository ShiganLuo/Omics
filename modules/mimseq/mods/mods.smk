from snakemake.logging import logger
import time
import os
import sys

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

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
        name = config.get("Params", {}).get("mimseq", {}).get("name", "tRNAseq"),
        min_cov = config.get("Params", {}).get("mimseq", {}).get("min_cov", 0.0005),
        misinc_thresh = config.get("Params", {}).get("mimseq", {}).get("misinc_thresh", 0.1),
    threads: 4
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"mimseq_mods_{current_time}.sh")

            cmd = [
                "/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin/python", params.script,
                "--out", outdir + "/",
                "--name", params.name,
                "--threads", str(threads),
                "--min-cov", str(params.min_cov),
                "--misinc-thresh", str(params.misinc_thresh),
            ]
            if config.get("Params", {}).get("mimseq", {}).get("no_cca", False):
                cmd += ["--cca"]
            if config.get("Params", {}).get("mimseq", {}).get("remap", True):
                cmd += ["--remap"]
            if config.get("Params", {}).get("mimseq", {}).get("crosstalks", False):
                cmd += ["--crosstalks"]

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("export PATH=/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin:$PATH\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"mimseq_mods failed: {e}\n")
            raise

rule mimseq_mods_result:
    """Result rule for mods module."""
    input:
        mods_done = outdir + "/mods.done",
