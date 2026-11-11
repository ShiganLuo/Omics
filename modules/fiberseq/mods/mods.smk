include: "../../common/common.smk"

import os
import time

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
# col4 subfamily / col7 family / col8 class), coverage-weighted collapse
# to Params.mods.te_level, sigtest at that level, and a per-significant-
# group top-loci BED (browsing shortlist for regionplot — the family-level
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
p = config.get("Params", {}).get("fiberseq_mods", {})

_mod_type = p.get("mod_type", "m6a")            # m6a | cpg
_metric = p.get("metric", f"{_mod_type}_per_fiber")
# region axis: fire_peaks (default) | te (RepeatMasker loci collapsed per level)
_region_axis = p.get("region_axis", "fire_peaks")
_te_level = p.get("te_level", "subfamily")      # subfamily (col4) | family (col7)
_min_fibers = int(p.get("min_fibers", 3))
_te_bed = config.get("te_bed") or config.get("genome", {}).get("repeat_bed")
_is_te = _region_axis == "te"
# axis-dependent long-table path: TE mode namespaces outputs under te_<level>
_long_name = f"te_{_mod_type}_{_te_level}_long.tsv" if _is_te else f"{_metric}_long.tsv"
_loci_name = f"te_{_mod_type}_loci.tsv" if _is_te else None
# TE mode reuses the longtable rule name's resources but scales up (4M-loci sweeps)
_longtable_rule_key = "fiberseq_mods_te_longtable" if _is_te else "fiberseq_mods_longtable"
_tmp_dir = config.get("tmp_dir") or os.path.join(outdir, "tmp")

sigtest_contrasts = [str(c) for c in (config.get("contrasts") or [])]


def _group_samples_spec() -> str:
    """'gA:s1+s2,gB:s3' from config group_samples.

    node.py derives it from the meta design (group_pairs) — the meta is
    the single source of grouping for this module.
    """
    gs = dict(p.get("group_samples") or {})
    if not gs:
        raise ValueError(
            "fiberseq_mods needs Params.fiberseq_mods.group_samples: fill "
            "the meta design/group columns (node.py derives them) and rerun.")
    return ",".join(f"{g}:{'+'.join(ss)}" for g, ss in gs.items())


def _mods_longtable_inputs(wildcards):
    """Inputs for the active region axis (DAG-time file declaration).

    Grouping arrives via Params (group_samples derived from the meta
    design by node.py) — no design file input.
    """
    files = {
        "mods": expand(indir + "/{sample_id}/{sample_id}." + _mod_type + ".bed.gz",
                       sample_id=samples),
        "msp": expand(indir + "/{sample_id}/{sample_id}.msp.bed.gz",
                      sample_id=samples),
    }
    if _is_te:
        files["te_bed"] = _require(_te_bed, "genome.repeat_bed", "fiberseq_mods")
    else:
        files["peaks"] = expand(indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
                                sample_id=samples)
    return files


rule fiberseq_mods_longtable:
    """Aggregate per-fiber modification calls over the region axis (long table).

    fire_peaks axis: per-peak values, feature = peak locus.
    te axis: per-locus values collapsed coverage-weighted to te_level,
    plus the per-loci table consumed by fiberseq_mods_top_loci.
    """
    input:
        unpack(_mods_longtable_inputs),
    output:
        longtable = outdir + f"/{_long_name}",
        loci = ([outdir + "/" + _loci_name] if _is_te else []),
    log:
        logdir + f"/mods/fiberseq_mods_longtable_{_long_name}.log"
    threads: rule_res(_longtable_rule_key).get("threads", 1)
    resources:
        mem_mb = rule_res(_longtable_rule_key).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_mods_longtable")
    params:
        script = (ROOT_DIR + "/modules/fiberseq/mods/bin/te_modification.py" if _is_te
                  else ROOT_DIR + "/modules/fiberseq/mods/bin/extract_modification_longtable.py"),
        metric = _metric,
        group_spec = _group_samples_spec(),
        is_te = _is_te,
        mod_type = _mod_type,
        te_level = _te_level,
        min_fibers = _min_fibers,
        loci_name = _loci_name or "",
        tmp_dir = _tmp_dir,
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
                    "--group-samples", params.group_spec,
                    "--te-level", params.te_level,
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
                    "--group-samples", params.group_spec,
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


# Directory where the sigtest module's outputs land (subworkflow wires
# sigtest's outdir to <mods>/sigtest and its input table to our long table;
# cross-module contract is FILES ONLY — no script references either way).
SIGTEST_DIR = config.get("sigtest_dir") or os.path.join(outdir, "sigtest")


def _sigtest_tsv(wildcards):
    """The sigtest module's result table for one contrast."""
    return os.path.join(SIGTEST_DIR, wildcards.contrast, "significance_test.tsv")


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
        p = p,
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
                "--top-n", str(params.p.get("top_n", 0)),
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
        loci = outdir + "/" + _loci_name if _is_te else [],
    output:
        tsv = [outdir + "/{contrast}/top_loci.tsv"] if _is_te else [],
        bed = [outdir + "/{contrast}/top_loci.bed"] if _is_te else [],
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
        group_spec = _group_samples_spec(),
        p = p,
        mod_type = _mod_type,
        min_fibers = _min_fibers,
        is_te = _is_te,
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
                # group spec restricted to the two compared groups
                grp_pair = (params.p.get("groups_by_contrast", {}).get(str(wildcards.contrast))
                            or params.p.get("groups") or [
                                tok.split(":", 1)[0] for tok in params.group_spec.split(",")])
                by_group = {
                    tok.split(":", 1)[0]: tok.split(":", 1)[1]
                    for tok in params.group_spec.split(",")}
                grp_spec = ",".join(f"{g}:{by_group[g]}" for g in grp_pair[:2])
                cmd = [
                    "python", params.top_loci_script,
                    "--sig", str(input.sig),
                    "--loci", str(input.loci[0]),
                    "--groups", grp_spec,
                    "--value-suffix", f"{params.mod_type}_per_fiber",
                    "--top-n", str(params.p.get("loci_top_n", 5)),
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
        longtable = outdir + f"/{_long_name}",
        sig_tsv = expand(outdir + "/{contrast}/significant_objects.tsv",
                         contrast=sigtest_contrasts),
        sig_bed = expand(outdir + "/{contrast}/significant_objects.bed",
                         contrast=sigtest_contrasts),
        top_loci = (expand(outdir + "/{contrast}/top_loci.bed",
                           contrast=sigtest_contrasts) if _is_te else []),
