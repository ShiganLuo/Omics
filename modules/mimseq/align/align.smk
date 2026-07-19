include: "../../common/common.smk"

from snakemake.logging import logger
import time
import os
import sys

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

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
        name = config.get("Params", {}).get("mimseq", {}).get("name", "tRNAseq"),
        mismatches = config.get("Params", {}).get("mimseq", {}).get("max_mismatches", 0.075),
        remap = config.get("Params", {}).get("mimseq", {}).get("remap", True),
        keep_temp = config.get("Params", {}).get("mimseq", {}).get("keep_temp", False),
    threads: 10
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"mimseq_align_{current_time}.sh")

            cmd = [
                "/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin/python", params.script,
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

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("export PATH=/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin:$PATH\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"mimseq_align failed: {e}\n")
            raise

rule mimseq_align_result:
    """Result rule for align module."""
    input:
        align_done = outdir + "/align.done",
