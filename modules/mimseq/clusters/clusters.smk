include: "../../common/common.smk"


ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# 参数精确到 rule:Params.mimseq.mimseq_clusters.<key>(rule_param 见 common.smk;
# 无 module 级共享回退)。
_mimseq_params = config.get("Params", {}).get("mimseq", {})

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
        name = rule_param(_mimseq_params, "mimseq_clusters", "name", "tRNAseq"),
        cluster_id = rule_param(_mimseq_params, "mimseq_clusters", "cluster_id", 0.97),
        cov_diff = rule_param(_mimseq_params, "mimseq_clusters", "cov_diff", 0.5),
    threads: rule_res("mimseq_clusters").get("threads", 4)
    benchmark: logdir + "/mimseq_clusters.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_clusters", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../mimseq.yaml"
    container:
        sif("../mimseq.yaml", "mimseq_clusters")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(RECEIPT_DIR, f"mimseq_clusters_{current_time}.sh")
            plain_script = os.path.join(outdir, f"mimseq_clusters.sh")

            cmd = [
                "python", params.script,
                "--out", outdir + "/",
                "--name", params.name,
                "--threads", str(threads),
                "--cluster-id", str(params.cluster_id),
                "--cov-diff", str(params.cov_diff),
            ]
            if rule_param(_mimseq_params, "mimseq_clusters", "cluster", True):
                cmd += ["--cluster"]

            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule mimseq_clusters_result:
    """Result rule for clusters module."""
    input:
        clusters_done = outdir + "/clusters.done",

    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_clusters_result", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
