include: "../common/common.smk"

import json

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "")
tissues = config.get("tissues", [])
pseudobulk_dir = config.get("outdir_pseudobulk", "")

tep_params = config.get("Params", {}).get("teprofile", {})
design_source = config.get("design", {}) or {}
tissue_samples = config.get("tissue_samples", {}) or {}


def _tissue_has_design(t):
    return bool((design_source.get(t) or {}).get("comparisons"))


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
            script = receipt_script(f"te_profile_{current_time}.sh")
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


def _te_fraction_inputs(wildcards):
    """Pseudo-bulk report per design tissue (declared output of
    scanpy_pseudobulk; its counts tables sit in the same out_root and are
    what the engine reads)."""
    files = []
    for t in tissues:
        if _tissue_has_design(t):
            files.append(_require(
                os.path.join(pseudobulk_dir, t, "scTE", "pseudobulk_report.html"),
                f"outdir_pseudobulk/{t}/scTE/pseudobulk_report", tool="teprofile"))
    return files


rule te_fraction_cross_tissue:
    """Per-comparison TE fraction change across tissues (compare design):
    one sorted bar panel per single-comparison tissue, one heatmap (columns
    keyed 1..n + text legend) per multi-comparison tissue. Effect = log2FC
    of the group-level TE fraction (ratio of summed pseudo-bulk counts)."""
    input:
        _te_fraction_inputs,
    output:
        fig = outdir + "/te_fraction_cross_tissue.png",
        design = outdir + "/te_fraction_compare_design.json",
    log:
        logdir + "/te_fraction_cross_tissue.log"
    threads: rule_res("te_fraction_cross_tissue").get("threads", 4)
    benchmark: logdir + "/te_fraction_cross_tissue.benchmark.tsv"
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "te_fraction_cross_tissue", attempt=attempt).get("mem_mb", 8000),
    conda:
        "teprofile.yaml"
    container:
        sif("teprofile.yaml", "te_fraction_cross_tissue")
    params:
        script = os.path.join(ROOT_DIR, "modules", "teprofile", "bin", "te_fraction_cross_tissue.py"),
        p = tep_params,
        tissues = tissues,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("te_fraction_cross_tissue", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start te_fraction_cross_tissue at {current_time}")
            tissue_design = {}
            for t in params.tissues:
                if not _tissue_has_design(t):
                    continue
                tissue_design[t] = comparison_design(
                    design_source, t, "pairwise",
                    tissue_samples=tissue_samples)
            if not tissue_design:
                raise ValueError(
                    "sample_params.design has no comparison for any "
                    f"configured tissue {params.tissues} — fill "
                    "sample_params.design.<tissue>.comparisons")
            with open(str(output.design), "w") as f:
                json.dump(tissue_design, f, indent=2)
            script = receipt_script(f"te_fraction_cross_tissue_{current_time}.sh")
            plain_script = os.path.join(outdir, "te_fraction_cross_tissue.sh")
            cmd = [
                "python", params.script,
                "--root", str(OUTDIR_ROOT),
                "--design", str(output.design),
                "--out", str(output.fig),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write("echo 'te_fraction_cross_tissue was completed successfully'\n")
            rule_logger.info("Executing: " + " ".join(str(x) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info("te_fraction_cross_tissue completed")
        except Exception as e:
            rule_logger.error(f"Error occurred during te_fraction_cross_tissue: {e}")
            raise RuntimeError(f"te_fraction_cross_tissue error: {e}, you can check the log file {log_path} for more details.")


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
        cross_fig = outdir + "/te_fraction_cross_tissue.png" if any(_tissue_has_design(t) for t in tissues) else [],
        cross_design = outdir + "/te_fraction_compare_design.json" if any(_tissue_has_design(t) for t in tissues) else [],
        done = outdir + "/te_profile.done",
