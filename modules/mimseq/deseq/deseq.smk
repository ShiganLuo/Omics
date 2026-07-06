from snakemake.logging import logger
import time
import os
import sys

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

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
        control_cond = config.get("Params", {}).get("mimseq", {}).get("control_cond", ""),
        p_adj = config.get("Params", {}).get("mimseq", {}).get("p_adj", 0.05),
        mito_trnas = config.get("genome", {}).get("mito_trnas", ""),
    threads: 4
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"mimseq_deseq_{current_time}.sh")

            cmd = [
                "/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin/python", params.script,
                "--out", outdir + "/",
                "--p-adj", str(params.p_adj),
            ]
            if params.control_cond:
                cmd += ["--control-cond", params.control_cond]
            if params.mito_trnas:
                cmd += ["--mito-trnas", params.mito_trnas]

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("export PATH=/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin:$PATH\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"mimseq_deseq failed: {e}\n")
            raise

rule mimseq_deseq_result:
    """Result rule for deseq module."""
    input:
        deseq_done = outdir + "/deseq.done",
