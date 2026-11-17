
include: "../common/common.smk"


ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
input_h5ad = config.get("input_h5ad", "")
visium_h5 = config.get("visium_h5", "")
spatial_dir = config.get("spatial_dir", "")
params = config.get("Params", {}).get("spatial", {})
python = config.get("Procedure", {}).get("python") or "python"
script = os.path.join(ROOT_DIR, "modules", "spatial_scanpy", "bin", "spatial_transcriptomics.py")
qc_h5ad = outdir + "/spatial/qc/filtered.h5ad"
cluster_h5ad = outdir + "/spatial/cluster/clustered.h5ad"
advanced_h5ad = outdir + "/spatial/advanced/spatial_autocorrelation.h5ad"

# 参数精确到 rule:Params.spatial.<rule_name>.<key>(rule_param 见 common.smk;
# 原 Params.spatial 平铺 key 与 advanced 子段均挂到消费 rule 名下,无共享回退:
# min_genes 等 → spatial_scanpy_qc;n_pcs/resolution 等 → spatial_scanpy_cluster;
# spatial_autocorrelation 开关 → spatial_scanpy_advanced)。

if not input_h5ad and not visium_h5:
    raise ValueError("spatial transcriptomics requires input_h5ad or visium_h5")

def _common_command(mode, output, extra):
    cmd = [python, script, "--mode", mode]
    if input_h5ad:
        cmd += ["--input-h5ad", input_h5ad]
    else:
        cmd += ["--visium-h5", visium_h5]
        if spatial_dir:
            cmd += ["--spatial-dir", spatial_dir]
    cmd += ["--output", output]
    cmd += extra
    return cmd

rule spatial_scanpy_qc:
    input:
        h5ad=input_h5ad if input_h5ad else visium_h5,
        spatial=lambda wc: spatial_dir if spatial_dir else []
    output:
        h5ad=qc_h5ad,
        metrics=outdir + "/spatial/qc/qc_metrics.tsv"
    log:
        logdir + "/spatial/spatial_scanpy_qc.log"
    threads: rule_res("spatial_scanpy_qc").get("threads", 4)
    benchmark: logdir + "/spatial/spatial_scanpy_qc.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "spatial_scanpy_qc", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "spatial_scanpy.yaml"
    container:
        sif("spatial_scanpy.yaml", "spatial_scanpy_qc")
    params:
        # 参数精确到 rule:Params.spatial.spatial_scanpy_qc.<key>
        command=lambda wc: _common_command("qc", qc_h5ad, ["--metrics", outdir + "/spatial/qc/qc_metrics.tsv", "--min-genes", str(rule_param(params, "spatial_scanpy_qc", "min_genes", 100)), "--max-pct-mt", str(rule_param(params, "spatial_scanpy_qc", "max_pct_mt", 25)), "--n-top-genes", str(rule_param(params, "spatial_scanpy_qc", "n_top_genes", 3000))])
    run:
        log_path=str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True); os.makedirs(os.path.dirname(output.h5ad), exist_ok=True)
            stamp=time.strftime("%Y%m%d_%H%M%S", time.localtime()); command_script = os.path.join(RECEIPT_DIR, f"spatial_qc_{stamp}.sh");plain_script=os.path.join(os.path.dirname(output.h5ad), f"spatial_qc.sh")
            with open(command_script,"w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(shlex.quote(str(x)) for x in params.command)+"\necho \"script completed successfully\"\n")
            shutil.copyfile(command_script, plain_script)
            shell(f"bash {shlex.quote(plain_script)} >> {shlex.quote(log_path)} 2>&1")
        except Exception as exc:
            rule_logger.error(f"spatial QC failed: {exc}")
            raise RuntimeError(f"spatial QC failed: {exc}, see {log_path} for details.")

rule spatial_scanpy_cluster:
    input: h5ad=qc_h5ad
    output: h5ad=cluster_h5ad, genes=outdir + "/spatial/cluster/spatial_markers.tsv"
    log: logdir + "/spatial/spatial_scanpy_cluster.log"
    threads: rule_res("spatial_scanpy_cluster").get("threads", 4)
    benchmark: logdir + "/spatial/spatial_scanpy_cluster.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "spatial_scanpy_cluster", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda: "spatial_scanpy.yaml"
    container:
        sif("spatial_scanpy.yaml", "spatial_scanpy_cluster")
    params:
        # 参数精确到 rule:Params.spatial.spatial_scanpy_cluster.<key>
        command=lambda wc: _common_command("cluster", cluster_h5ad, ["--genes", outdir + "/spatial/cluster/spatial_markers.tsv", "--n-pcs", str(rule_param(params, "spatial_scanpy_cluster", "n_pcs", 50)), "--n-neighbors", str(rule_param(params, "spatial_scanpy_cluster", "n_neighbors", 15)), "--resolution", str(rule_param(params, "spatial_scanpy_cluster", "resolution", 0.8))])
    run:
        log_path=str(log)
        try:
            os.makedirs(os.path.dirname(log_path),exist_ok=True);os.makedirs(os.path.dirname(output.h5ad),exist_ok=True);stamp=time.strftime("%Y%m%d_%H%M%S",time.localtime());command_script = os.path.join(RECEIPT_DIR, f"spatial_cluster_{stamp}.sh");plain_script=os.path.join(os.path.dirname(output.h5ad), f"spatial_cluster.sh")
            with open(command_script,"w") as handle: handle.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(shlex.quote(str(x)) for x in params.command)+"\necho \"script completed successfully\"\n")
            shutil.copyfile(command_script, plain_script)
            shell(f"bash {shlex.quote(plain_script)} >> {shlex.quote(log_path)} 2>&1")
        except Exception as exc:
            rule_logger.error(f"spatial clustering failed: {exc}")
            raise RuntimeError(f"spatial clustering failed: {exc}, see {log_path} for details.")

# 参数精确到 rule:开关挂在 spatial_scanpy_advanced 名下(原 advanced 子段),
# parse 期决定该 rule 是否定义。
if rule_param(params, "spatial_scanpy_advanced", "spatial_autocorrelation", True):
    rule spatial_scanpy_advanced:
        input: h5ad=cluster_h5ad
        output: h5ad=advanced_h5ad
        log: logdir + "/spatial/spatial_scanpy_advanced.log"
        threads: rule_res("spatial_scanpy_advanced").get("threads", 4)
        conda: "spatial_scanpy.yaml"
        container:
            sif("spatial_scanpy.yaml", "spatial_scanpy_advanced")
        params: command=lambda wc: _common_command("advanced", advanced_h5ad, [])
        run:
            log_path=str(log)
            try:
                os.makedirs(os.path.dirname(log_path),exist_ok=True);os.makedirs(os.path.dirname(output.h5ad),exist_ok=True);stamp=time.strftime("%Y%m%d_%H%M%S",time.localtime());command_script = os.path.join(RECEIPT_DIR, f"spatial_advanced_{stamp}.sh");plain_script=os.path.join(os.path.dirname(output.h5ad), f"spatial_advanced.sh")
                with open(command_script,"w") as handle: handle.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(shlex.quote(str(x)) for x in params.command)+"\necho \"script completed successfully\"\n")
                shutil.copyfile(command_script, plain_script)
                shell(f"bash {shlex.quote(plain_script)} >> {shlex.quote(log_path)} 2>&1")
            except Exception as exc:
                rule_logger.error(f"spatial advanced analysis failed: {exc}")
                raise RuntimeError(f"spatial advanced analysis failed: {exc}, see {log_path} for details.")
