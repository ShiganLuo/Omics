include: "../../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# 参数精确到 rule:config 全量路径读取 Params.mimseq.mimseq_clusters.<key>
# (无 module 级共享回退)。

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
        name = config.get("Params", {}).get("mimseq", {}).get("mimseq_clusters", {}).get("name", "tRNAseq"),
        cluster_id = config.get("Params", {}).get("mimseq", {}).get("mimseq_clusters", {}).get("cluster_id", 0.97),
        cov_diff = config.get("Params", {}).get("mimseq", {}).get("mimseq_clusters", {}).get("cov_diff", 0.5),
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
        open(log_path, "a").close()
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
            if config.get("Params", {}).get("mimseq", {}).get("mimseq_clusters", {}).get("cluster", True):
                cmd += ["--cluster"]

            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule mimseq_clusters_result:
    """Result rule for clusters module."""
    input:
        clusters_done = outdir + "/clusters.done",
