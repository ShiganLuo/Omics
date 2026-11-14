include: "../../common/common.smk"

import json  # noqa: E402  (not in common.smk; group_samples -> engine CLI JSON)


# ---------------------------------------------------------------------------
# mods: prepare the modification long table (m6A/CpG per fiber over FIRE
# peaks or TE loci) for the sigtest module, then postprocess sigtest's
# results (significant objects split; TE-axis top-loci BED for regionplot).
#
# The sigtest module itself is a SEPARATE module consumed via subworkflow
# `use rule` — cross-module contract is files only (modules.md 模块不得
# 引用外部脚本): our long table -> sigtest's table input; sigtest's
# significance_test.tsv -> our split/top-loci inputs (sigtest_dir config).
#
# region_axis=te switches the region axis to transposable elements:
# per-locus quantification over a RepeatMasker TE BED (genome.repeat_bed,
# col4 subfamily / col7 family / col8 class); the long table carries ALL
# annotation levels, sigtest runs once per Params.sigtest.te_levels entry
# (its feature column picks the level), and a per-significant-group
# top-loci BED (browsing shortlist for regionplot — the family-level
# p-value is the only statistical assertion; loci are ranked, not tested).
#
# Upstream contract (subworkflow passes):
#   indir            = elements dir (common/4_elements) — {sample}.m6a.bed.gz,
#                      .cpg.bed.gz, .msp.bed.gz, .fire_peaks.bed
#   te_bed           = sorted-or-not TE BED (te mode; sorted copy cached in tmp)
# Group-level module: no {sample_id} wildcard, logdir gets {logdir}/group.
# ---------------------------------------------------------------------------

indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/mods")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")


# 指令串联所需的派生名(output/log/threads 指令无法调用函数);
# 配置形状由 schema(run.py)与 node.py 保证,模块只校验输入文件。
region_axis = config.get("Params", {}).get("fiberseq_mods", {}).get("region_axis") or {}
mod_type = config.get("Params", {}).get("fiberseq_mods", {}).get("mod_type", "m6a")
te_axis = region_axis.get("default", "fire_peaks") == "te"
sigtest_contrasts = [str(c) for c in (config.get("contrasts") or [])]


def _mods_longtable_inputs(wildcards):
    """Inputs for the active region axis (DAG-time file declaration).

    Grouping arrives via Params (group_samples derived from the meta
    design by node.py) — no design file input.
    """
    files = {
        "mods": expand(indir + "/{sample_id}/{sample_id}." + mod_type + ".bed.gz",
                       sample_id=samples),
        "msp": expand(indir + "/{sample_id}/{sample_id}.msp.bed.gz",
                      sample_id=samples),
    }
    if te_axis:
        te_bed = config.get("te_bed") or config.get("genome", {}).get("repeat_bed")
        files["te_bed"] = _require(te_bed, "genome.repeat_bed", "fiberseq_mods")
    else:
        files["peaks"] = expand(indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
                                sample_id=samples)
    return files


rule fiberseq_mods_longtable:
    """Aggregate per-fiber modification calls over the region axis (long table).

    fire_peaks axis: per-peak values, feature = peak locus.
    te axis: per-locus values with subfamily/family/class annotation
    columns (sigtest's feature column picks the tested level), plus the
    per-loci table consumed by fiberseq_mods_top_loci.
    """
    input:
        unpack(_mods_longtable_inputs),
    output:
        longtable = (outdir + f"/te_{mod_type}_long.tsv" if te_axis
                     else outdir + f"/{region_axis['options']['fire_peaks'].get('metric', f'{mod_type}_per_fiber')}_long.tsv"),
        loci = ([outdir + f"/te_{mod_type}_loci.tsv"] if te_axis else []),
    log:
        logdir + (f"/mods/fiberseq_mods_longtable_te_{mod_type}_long.tsv.log" if te_axis
                 else f"/mods/fiberseq_mods_longtable_{region_axis['options']['fire_peaks'].get('metric', f'{mod_type}_per_fiber')}_long.tsv.log")
    threads: rule_res("fiberseq_mods_te_longtable" if te_axis else "fiberseq_mods_longtable").get("threads", 1)
    resources:
        mem_mb = rule_res("fiberseq_mods_te_longtable" if te_axis else "fiberseq_mods_longtable").get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_mods_longtable")
    params:
        script = (ROOT_DIR + "/modules/fiberseq/mods/bin/te_modification.py" if te_axis
                  else ROOT_DIR + "/modules/fiberseq/mods/bin/extract_modification_longtable.py"),
        metric = region_axis["options"]["fire_peaks"].get("metric", f"{mod_type}_per_fiber"),
        # <- Params.fiberseq_mods.group_samples (node.py 从 meta 设计派生)
        group_samples = config.get("Params", {}).get("fiberseq_mods", {}).get("group_samples", {}),
        is_te = te_axis,
        mod_type = mod_type,
        min_fibers = int(region_axis["options"]["te"].get("min_fibers", 3)),
        loci_name = (f"te_{mod_type}_loci.tsv" if te_axis else ""),
        tmp_dir = config.get("tmp_dir") or os.path.join(outdir, "tmp"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_mods_longtable", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start mods long table at {current_time}")
            sample_outdir = os.path.dirname(str(output.longtable))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_mods_longtable_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "fiberseq_mods_longtable.sh")
            if params.is_te:
                cmd = [
                    "python", params.script,
                    "--mods", ",".join(
                        f"{s}:{indir}/{s}/{s}.{params.mod_type}.bed.gz" for s in samples),
                    "--msp", ",".join(
                        f"{indir}/{s}/{s}.msp.bed.gz" for s in samples),
                    "--te", str(input.te_bed),
                    "--group-samples", json.dumps(params.group_samples),
                    "--mod-type", params.mod_type,
                    "--min-fibers", str(params.min_fibers),
                    "--out-long", str(output.longtable),
                    "--out-loci", os.path.join(outdir, params.loci_name),
                    "--tmp-dir", params.tmp_dir,
                    "-l", log_path,
                ]
            else:
                cmd = [
                    "python", params.script,
                    "--mods", ",".join(
                        f"{s}:{indir}/{s}/{s}.{params.mod_type}.bed.gz" for s in samples),
                    "--peaks", ",".join(
                        f"{indir}/{s}/{s}.fire_peaks.bed" for s in samples),
                    "--msp", ",".join(
                        f"{indir}/{s}/{s}.msp.bed.gz" for s in samples),
                    "--group-samples", json.dumps(params.group_samples),
                    "--metric", params.metric,
                    "--out", str(output.longtable),
                    "-l", log_path,
                ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "fiberseq_mods_longtable at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_mods_longtable: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_mods_longtable: {e}, you can check the log file {log_path} for more details.\n")


def _sigtest_tsv(wildcards):
    """The sigtest module's result table for one contrast.

    The subworkflow wires sigtest's outdir to <mods>/sigtest (file-only
    cross-module contract); config.get at point of use, no global.
    """
    sigtest_dir = config.get("sigtest_dir") or os.path.join(outdir, "sigtest")
    return os.path.join(sigtest_dir, wildcards.contrast, "significance_test.tsv")


rule fiberseq_mods_split:
    """Split sigtest results into significant objects (TSV + regionplot BED).

    Consumes the sigtest module's significance_test.tsv; the split engine
    and the BED writer live in THIS module's bin/ (self-contained).
    """
    input:
        sig = _sigtest_tsv,
    output:
        sig_tsv = outdir + "/{contrast}/significant_objects.tsv",
        sig_bed = outdir + "/{contrast}/significant_objects.bed",
        rest_tsv = outdir + "/{contrast}/nonsignificant_objects.tsv",
    log:
        logdir + "/mods/split_{contrast}.log"
    threads: rule_res("fiberseq_mods_split").get("threads", 1)
    resources:
        mem_mb = rule_res("fiberseq_mods_split").get("mem_mb", 4000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_mods_split")
    params:
        split_script = ROOT_DIR + "/modules/fiberseq/mods/bin/split_significant.py",
        # <- Params.fiberseq_mods.region_axis.options.fire_peaks.top_n
        top_n = int(region_axis["options"]["fire_peaks"].get("top_n", 0)),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_mods_split", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Splitting sigtest results for contrast {wildcards.contrast}")
            sample_outdir = os.path.dirname(str(output.sig_tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_mods_split_{wildcards.contrast}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"fiberseq_mods_split_{wildcards.contrast}.sh")
            cmd = [
                "python", params.split_script,
                "-i", str(input.sig),
                "--out-tsv", str(output.sig_tsv),
                "--out-bed", str(output.sig_bed),
                "--out-rest-tsv", str(output.rest_tsv),
                "--top-n", str(params.top_n),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "fiberseq_mods_split for {wildcards.contrast} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_mods_split for {wildcards.contrast}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_mods_split for {wildcards.contrast}: {e}, you can check the log file {log_path} for more details.\n")


rule fiberseq_mods_top_loci:
    """Rank top loci within significant TE subfamilies/families (TE axis only).

    A regionplot browsing shortlist — the subfamily-level p-value from the
    sigtest module is the only statistical assertion; loci are ranked by
    |group mean diff|, not tested.
    """
    input:
        sig = _sigtest_tsv,
        loci = outdir + f"/te_{mod_type}_loci.tsv" if te_axis else [],
    output:
        tsv = [outdir + "/{contrast}/top_loci.tsv"] if te_axis else [],
        bed = [outdir + "/{contrast}/top_loci.bed"] if te_axis else [],
    log:
        logdir + "/mods/top_loci_{contrast}.log"
    threads: rule_res("fiberseq_mods_top_loci").get("threads", 1)
    resources:
        mem_mb = rule_res("fiberseq_mods_top_loci").get("mem_mb", 4000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_mods_top_loci")
    params:
        top_loci_script = ROOT_DIR + "/modules/fiberseq/mods/bin/te_top_loci.py",
        # <- Params.fiberseq_mods.group_samples (node.py 从 meta 设计派生)
        group_samples = config.get("Params", {}).get("fiberseq_mods", {}).get("group_samples", {}),
        grp_by_contrast = config.get("Params", {}).get("fiberseq_mods", {}).get("groups_by_contrast", {}),
        groups = config.get("Params", {}).get("fiberseq_mods", {}).get("groups", []),
        mod_type = mod_type,
        min_fibers = int(region_axis["options"]["te"].get("min_fibers", 3)),
        is_te = te_axis,
    run:
        if not params.is_te:
            # fire_peaks axis declares empty outputs; nothing to rank
            for out in list(output.tsv) + list(output.bed):
                open(str(out), "w").close()
        else:
            log_path = str(log)
            open(log_path, "w").close()
            rule_logger = setup_logger("fiberseq_mods_top_loci", log_file=log_path)
            try:
                current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
                rule_logger.info(f"Ranking top loci for contrast {wildcards.contrast}")
                sample_outdir = os.path.dirname(str(output.tsv[0]))
                os.makedirs(sample_outdir, exist_ok=True)
                script = os.path.join(RECEIPT_DIR, f"fiberseq_mods_top_loci_{wildcards.contrast}_{current_time}.sh")
                plain_script = os.path.join(sample_outdir, f"fiberseq_mods_top_loci_{wildcards.contrast}.sh")
                # the two compared groups for this contrast (node.py
                # groups_by_contrast -> global groups -> all groups), their
                # sample lists picked from group_samples
                grp_pair = (params.grp_by_contrast.get(str(wildcards.contrast))
                            or params.groups or list(params.group_samples))[:2]
                grp_spec = json.dumps({g: params.group_samples[g] for g in grp_pair})
                cmd = [
                    "python", params.top_loci_script,
                    "--sig", str(input.sig),
                    "--loci", str(input.loci[0]),
                    "--groups", grp_spec,
                    "--value-suffix", f"{params.mod_type}_per_fiber",
                    "--top-n", str(params.loci_top_n),
                    "--min-fibers", str(params.min_fibers),
                    "--out-tsv", str(output.tsv[0]),
                    "--out-bed", str(output.bed[0]),
                ]
                with open(script, "w") as f:
                    f.write("#!/bin/bash\n")
                    f.write("set -euo pipefail\n")
                    f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                    f.write(f'echo "fiberseq_mods_top_loci for {wildcards.contrast} completed successfully"\n')
                shutil.copyfile(script, plain_script)
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            except Exception as e:
                rule_logger.error(f"Error occurred during fiberseq_mods_top_loci for {wildcards.contrast}: {e}\n")
                raise RuntimeError(f"Error occurred during fiberseq_mods_top_loci for {wildcards.contrast}: {e}, you can check the log file {log_path} for more details.\n")


rule fiberseq_mods_result:
    """Panel: aggregate mods outputs as a dependency point for use rule.

    significance_test.tsv itself belongs to the sigtest module (its own
    result panel aggregates it); this panel covers the mods-owned outputs.
    """
    input:
        longtable = (outdir + f"/te_{mod_type}_long.tsv" if te_axis
                     else outdir + f"/{region_axis['options']['fire_peaks'].get('metric', f'{mod_type}_per_fiber')}_long.tsv"),
        sig_tsv = expand(outdir + "/{contrast}/significant_objects.tsv",
                         contrast=sigtest_contrasts),
        sig_bed = expand(outdir + "/{contrast}/significant_objects.bed",
                         contrast=sigtest_contrasts),
        top_loci = (expand(outdir + "/{contrast}/top_loci.bed",
                           contrast=sigtest_contrasts) if te_axis else []),
