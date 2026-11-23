include: "../common/common.smk"

# ---------------------------------------------------------------------------
# teprofile: TE 表达三层刻画(layer 1 细胞类型本源特征 / layer 2 亚群离群
# 检测 / layer 3 分组重塑)。输入 = scRNAseq 工作流 combine 阶段的
# {tissue}_scTE_{merged,auto}.h5ad 对;一次 group 级运行覆盖全部组织。
# 产物由 done 标记跟踪(modules.md:禁用 directory())。本文件只做接线;
# 三层算法与指标定义见 bin/te_profile.py 与 README §2。
# ---------------------------------------------------------------------------
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "")
tissues = config.get("tissues", [])

tep_params = config.get("Params", {}).get("teprofile", {})


def _teprofile_inputs(wildcards):
    """Declare the h5ad pair per tissue (existence checked at DAG build)."""
    files = {"merged": [], "auto": []}
    for t in tissues:
        files["merged"].append(_require(
            os.path.join(indir, t, f"{t}_scTE_merged.h5ad"),
            f"indir/{t}_scTE_merged", tool="teprofile"))
        files["auto"].append(_require(
            os.path.join(indir, t, f"{t}_scTE_auto.h5ad"),
            f"indir/{t}_scTE_auto", tool="teprofile"))
    return files


rule te_profile:
    """Three-layer TE profile: baseline features, outlier subpopulations,
    group remodeling; one run for all tissues."""
    input:
        unpack(_teprofile_inputs),
    output:
        baseline = outdir + "/te_profile_baseline.csv",
        clusters = outdir + "/te_profile_clusters.csv",
        remodeling = outdir + "/te_profile_remodeling.csv",
        base_fig = expand(outdir + "/baseline_te_fraction_{tissue}.png", tissue=tissues) if tissues else [],
        comp_fig = expand(outdir + "/baseline_family_comp_{tissue}.png", tissue=tissues) if tissues else [],
        dev_fig = expand(outdir + "/cluster_deviation_{tissue}.png", tissue=tissues) if tissues else [],
        remod_fig = expand(outdir + "/remodeling_shift_{tissue}.png", tissue=tissues) if (tissues and (tep_params.get("groups") or tep_params.get("contrasts"))) else [],
        done = outdir + "/te_profile.done",
    log:
        logdir + "/te_profile.log"
    threads: rule_res("te_profile").get("threads", 4)
    benchmark: logdir + "/te_profile.benchmark.tsv"
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "te_profile", input_size=input_size_named(input, 'merged'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "teprofile.yaml"
    container:
        sif("teprofile.yaml", "te_profile")
    params:
        script = os.path.join(ROOT_DIR, "modules", "teprofile", "bin", "te_profile.py"),
        p = tep_params,
        tissues = tissues,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("te_profile", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start te_profile at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"te_profile_{current_time}.sh")
            plain_script = os.path.join(outdir, "te_profile.sh")
            cmd = [
                "python", params.script,
                "--indir", str(indir),
                "--outdir", str(outdir),
                "--tissues", ",".join(params.tissues),
            ]
            for key, flag in (("fe_cut", "--fe-cut"), ("n_boot", "--n-boot"),
                              ("seed", "--seed"), ("groups", "--groups")):
                val = params.p.get(key)
                if val:
                    cmd += [flag, str(val)]
            for c in (params.p.get("contrasts") or []):
                cmd += ["--contrasts", str(c)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write("echo 'te_profile was completed successfully'\n")
            rule_logger.info("Executing: " + " ".join(str(x) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                open(str(output.done), "w").close()
            rule_logger.info("te_profile completed")
        except Exception as e:
            rule_logger.error(f"Error occurred during te_profile: {e}")
            raise RuntimeError(f"te_profile error: {e}, you can check the log file {log_path} for more details.")


rule teprofile_result:
    """Panel rule: aggregate deliverables as the node.py registration point."""

    input:
        baseline = outdir + "/te_profile_baseline.csv",
        clusters = outdir + "/te_profile_clusters.csv",
        remodeling = outdir + "/te_profile_remodeling.csv",
        base_fig = expand(outdir + "/baseline_te_fraction_{tissue}.png", tissue=tissues) if tissues else [],
        comp_fig = expand(outdir + "/baseline_family_comp_{tissue}.png", tissue=tissues) if tissues else [],
        dev_fig = expand(outdir + "/cluster_deviation_{tissue}.png", tissue=tissues) if tissues else [],
        remod_fig = expand(outdir + "/remodeling_shift_{tissue}.png", tissue=tissues) if (tissues and (tep_params.get("groups") or tep_params.get("contrasts"))) else [],
        done = outdir + "/te_profile.done",
