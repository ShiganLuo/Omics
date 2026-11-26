include: "../common/common.smk"

# ---------------------------------------------------------------------------
# tedominant: TE 主导 cluster 鉴定(v2 反分母塌缩门槛)。
# 输入 = scRNAseq 工作流 combine 阶段产出的 {tissue}_scTE_{merged,auto}.h5ad
# 对(merged 带基因+TE 特征与 counts 层,auto 带 leiden/cell_type 注释与
# UMAP);一次 group 级运行覆盖全部组织,输出鉴定表 + 三类图。
# 区域/组织名在数据里,产物由 done 标记跟踪(modules.md:禁用 directory())。
# 本文件只做接线;判据与算法见 bin/te_dominant.py 与 README §2。
# ---------------------------------------------------------------------------
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "")
tissues = config.get("tissues", [])

tedom_params = config.get("Params", {}).get("tedominant", {})


def _tedominant_inputs(wildcards):
    """Declare the h5ad pair per tissue (existence checked at DAG build)."""
    files = {"merged": [], "auto": []}
    for t in tissues:
        files["merged"].append(_require(
            os.path.join(indir, t, f"{t}_scTE_merged.h5ad"),
            f"indir/{t}_scTE_merged", tool="tedominant"))
        files["auto"].append(_require(
            os.path.join(indir, t, f"{t}_scTE_auto.h5ad"),
            f"indir/{t}_scTE_auto", tool="tedominant"))
    return files


rule te_dominant:
    """Identify TE-dominant clusters per tissue with v2 gates; one run for
    all tissues (tables + burden bars + UMAP + marker dotplot)."""
    input:
        unpack(_tedominant_inputs),
    output:
        table = outdir + "/te_dominant_clusters.csv",
        burden = expand(outdir + "/cluster_burden_{tissue}.png", tissue=tissues) if tissues else [],
        umap = expand(outdir + "/umap_te_{tissue}.png", tissue=tissues) if tissues else [],
        dotplot = outdir + "/markers_dotplot.png",
        done = outdir + "/te_dominant.done",
    log:
        logdir + "/te_dominant.log"
    threads: rule_res("te_dominant").get("threads", 4)
    benchmark: logdir + "/te_dominant.benchmark.tsv"
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "te_dominant", input_size=input_size_named(input, 'merged'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "tedominant.yaml"
    container:
        sif("tedominant.yaml", "te_dominant")
    params:
        script = os.path.join(ROOT_DIR, "modules", "tedominant", "bin", "te_dominant.py"),
        p = tedom_params,
        tissues = tissues,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("te_dominant", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start te_dominant at {current_time}")
            script = receipt_script(f"te_dominant_{current_time}.sh")
            plain_script = os.path.join(outdir, "te_dominant.sh")
            cmd = [
                "python", params.script,
                "--indir", str(indir),
                "--outdir", str(outdir),
                "--tissues", ",".join(params.tissues),
            ]
            for key, flag in (("method", "--method"), ("min_cells", "--min-cells"), ("gate_frac", "--gate-frac")):
                val = params.p.get(key)
                if val is not None:
                    cmd += [flag, str(val)]
            if params.p.get("skip_dotplot"):
                cmd += ["--skip-dotplot"]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write("echo 'te_dominant was completed successfully'\n")
            rule_logger.info("Executing: " + " ".join(str(x) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                open(str(output.done), "w").close()
            rule_logger.info("te_dominant completed")
        except Exception as e:
            rule_logger.error(f"Error occurred during te_dominant: {e}")
            raise RuntimeError(f"te_dominant error: {e}, you can check the log file {log_path} for more details.")


rule tedominant_result:
    """Panel rule: aggregate deliverables as the node.py registration point."""

    input:
        table = outdir + "/te_dominant_clusters.csv",
        burden = expand(outdir + "/cluster_burden_{tissue}.png", tissue=tissues) if tissues else [],
        umap = expand(outdir + "/umap_te_{tissue}.png", tissue=tissues) if tissues else [],
        dotplot = outdir + "/markers_dotplot.png",
        done = outdir + "/te_dominant.done",
