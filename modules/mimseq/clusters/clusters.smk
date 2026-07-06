from snakemake.logging import logger
import time
import os
import sys

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule mimseq_clusters:
    """Split isodecoders and perform cluster deconvolution."""
    input:
        align_done = outdir + "/align.done",
    output:
        clusters_done = touch(outdir + "/clusters.done"),
    log:
        logdir + "/mimseq_clusters.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "clusters", "run.py"),
        name = config.get("Params", {}).get("mimseq", {}).get("name", "tRNAseq"),
        cluster_id = config.get("Params", {}).get("mimseq", {}).get("cluster_id", 0.97),
        cov_diff = config.get("Params", {}).get("mimseq", {}).get("cov_diff", 0.5),
    threads: 4
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"mimseq_clusters_{current_time}.sh")

            cmd = [
                "/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin/python", params.script,
                "--out", outdir + "/",
                "--name", params.name,
                "--threads", str(threads),
                "--cluster-id", str(params.cluster_id),
                "--cov-diff", str(params.cov_diff),
            ]
            if config.get("Params", {}).get("mimseq", {}).get("cluster", True):
                cmd += ["--cluster"]

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("export PATH=/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin:$PATH\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"mimseq_clusters failed: {e}\n")
            raise

rule mimseq_clusters_result:
    """Result rule for clusters module."""
    input:
        clusters_done = outdir + "/clusters.done",
