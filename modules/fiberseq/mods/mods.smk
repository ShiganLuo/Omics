include: "../../common/common.smk"

import os
import time

# ---------------------------------------------------------------------------
# mods: aggregate ft_extract modification calls (m6A/CpG) over FIRE peaks
# into a sigtest long table, test per-region group differences, and split
# significant regions into a regionplot-ready BED.
#
# Upstream contract (subworkflow passes):
#   indir            = elements dir (common/4_elements) — {sample}.m6a.bed.gz,
#                      .cpg.bed.gz, .msp.bed.gz, .fire_peaks.bed
#   quant_dir        = fiberseq_analysis/quant (producer config's outdir)
# Group-level module: no {sample_id} wildcard, logdir gets {logdir}/group.
# ---------------------------------------------------------------------------

indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/mods")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
p = config.get("Params", {}).get("fiberseq_mods", {})

# design TSV (sample/group[/contrast]) — same file diff uses
_design = p.get("design_tsv")
_mod_type = p.get("mod_type", "m6a")            # m6a | cpg
_metric = p.get("metric", f"{_mod_type}_per_fiber")

sigtest_contrasts = [str(c) for c in (config.get("contrasts") or [])]


def _mods_design(wildcards):
    """Required design TSV; missing key fails at DAG time naming it."""
    return _require(_design, "Params.mods.design_tsv", "fiberseq_mods")


rule fiberseq_mods_longtable:
    """Aggregate per-fiber modification calls over FIRE peaks (long table)."""
    input:
        design = _mods_design,
        mods = expand(indir + "/{sample_id}/{sample_id}." + _mod_type + ".bed.gz",
                      sample_id=samples),
        peaks = expand(indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
                       sample_id=samples),
        msp = expand(indir + "/{sample_id}/{sample_id}.msp.bed.gz",
                     sample_id=samples),
    output:
        longtable = outdir + f"/{_metric}_long.tsv",
    log:
        logdir + f"/mods/fiberseq_mods_longtable_{_metric}.log"
    threads: rule_res("fiberseq_mods_longtable").get("threads", 1)
    resources:
        mem_mb = rule_res("fiberseq_mods_longtable").get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_mods_longtable")
    params:
        script = ROOT_DIR + "/modules/fiberseq/mods/bin/extract_modification_longtable.py",
        metric = _metric,
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
            cmd = [
                "python", params.script,
                "--mods", ",".join(
                    f"{s}:{indir}/{s}/{s}.{_mod_type}.bed.gz" for s in samples),
                "--peaks", ",".join(
                    f"{indir}/{s}/{s}.fire_peaks.bed" for s in samples),
                "--msp", ",".join(
                    f"{indir}/{s}/{s}.msp.bed.gz" for s in samples),
                "--design", str(input.design),
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


rule fiberseq_mods_sigtest:
    """Per-region significance test on the mods long table (sigtest module)."""
    input:
        longtable = outdir + f"/{_metric}_long.tsv",
    output:
        tsv = outdir + "/{contrast}/significance_test.tsv",
        png = outdir + "/{contrast}/significance_bars.png",
        sig_tsv = outdir + "/{contrast}/significant_objects.tsv",
        sig_bed = outdir + "/{contrast}/significant_objects.bed",
        rest_tsv = outdir + "/{contrast}/nonsignificant_objects.tsv",
    log:
        logdir + "/mods/{contrast}.log"
    threads: rule_res("fiberseq_mods_sigtest").get("threads", 1)
    resources:
        mem_mb = rule_res("fiberseq_mods_sigtest").get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_mods_sigtest")
    params:
        sigtest_script = ROOT_DIR + "/modules/sigtest/bin/sigtest.py",
        split_script = ROOT_DIR + "/modules/fiberseq/mods/bin/split_significant.py",
        p = p,
        metric = _metric,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_mods_sigtest", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start mods sigtest for contrast {wildcards.contrast} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_mods_sigtest_{wildcards.contrast}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"fiberseq_mods_sigtest_{wildcards.contrast}.sh")
            sig_cmd = [
                "python", params.sigtest_script,
                "-i", str(input.longtable),
                "-o", sample_outdir,
                "-f", "feature",
                "--group-col", "group",
                "-v", "value",
                "-m", params.p.get("method", "auto"),
                "-c", params.p.get("correction", "bh"),
                "-a", str(params.p.get("alpha", 0.05)),
                "-y", params.p.get("ylabel", params.metric),
                "-l", log_path,
            ]
            groups = params.p.get("groups") or []
            if len(groups) == 2:
                sig_cmd += ["-g", str(groups[0]), str(groups[1])]
            if params.p.get("no_plot"):
                sig_cmd += ["--no-plot"]
            split_cmd = [
                "python", params.split_script,
                "-i", str(output.tsv),
                "--out-tsv", str(output.sig_tsv),
                "--out-bed", str(output.sig_bed),
                "--out-rest-tsv", str(output.rest_tsv),
                "--top-n", str(params.p.get("top_n", 0)),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in sig_cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in split_cmd) + "\n")
                f.write(f'echo "fiberseq_mods_sigtest for {wildcards.contrast} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_mods_sigtest for {wildcards.contrast}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_mods_sigtest for {wildcards.contrast}: {e}, you can check the log file {log_path} for more details.\n")


rule fiberseq_mods_result:
    """Panel: aggregate mods outputs as a dependency point for use rule."""
    input:
        longtable = outdir + f"/{_metric}_long.tsv",
        tsv = expand(outdir + "/{contrast}/significance_test.tsv",
                     contrast=sigtest_contrasts),
        sig_bed = expand(outdir + "/{contrast}/significant_objects.bed",
                         contrast=sigtest_contrasts),
