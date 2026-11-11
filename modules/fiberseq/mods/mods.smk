include: "../../common/common.smk"

import os
import time

# ---------------------------------------------------------------------------
# mods: aggregate ft_extract modification calls (m6A/CpG) over FIRE peaks
# into a sigtest long table, test per-region group differences, and split
# significant regions into a regionplot-ready BED.
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

# design TSV (sample/group[/contrast]) — same file diff uses
_design = p.get("design_tsv")
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


def _mods_design(wildcards):
    """Required design TSV; missing key fails at DAG time naming it."""
    return _require(_design, "Params.mods.design_tsv", "fiberseq_mods")


def _group_samples_spec() -> str:
    """'gA:s1+s2,gB:s3' from config group_samples (node.py writes it from
    the meta design); a hand-filled design_tsv still parses as fallback."""
    gs = dict(p.get("group_samples") or {})
    if not gs and p.get("design_tsv"):
        sample_group: dict = {}
        with open(p["design_tsv"]) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            low = [h.strip().lower() for h in hdr]
            si, gi = low.index("sample"), low.index("group")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(si, gi) and f[si].strip():
                    sample_group[f[si].strip()] = f[gi].strip()
        gs = {}
        for s, g in sample_group.items():
            gs.setdefault(g, []).append(s)
    if not gs:
        raise ValueError(
            "fiberseq_mods needs sample grouping: fill the meta design/group "
            "columns (drives Params.fiberseq_mods.group_samples) or set "
            "Params.fiberseq_mods.design_tsv.")
    return ",".join(f"{g}:{'+'.join(ss)}" for g, ss in gs.items())


def _mods_longtable_inputs(wildcards):
    """Inputs for the active region axis (DAG-time file declaration).

    Grouping arrives via Params (group_samples derived from the meta
    design by node.py) — a design file input only when hand-authored
    design_tsv is set.
    """
    files = {
        "mods": expand(indir + "/{sample_id}/{sample_id}." + _mod_type + ".bed.gz",
                       sample_id=samples),
        "msp": expand(indir + "/{sample_id}/{sample_id}.msp.bed.gz",
                      sample_id=samples),
    }
    if p.get("design_tsv"):
        files["design"] = _mods_design(wildcards)
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
            if _is_te:
                cmd = [
                    "python", params.script,
                    "--mods", ",".join(
                        f"{s}:{indir}/{s}/{s}.{_mod_type}.bed.gz" for s in samples),
                    "--msp", ",".join(
                        f"{indir}/{s}/{s}.msp.bed.gz" for s in samples),
                    "--te", str(input.te_bed),
                    "--group-samples", params.group_spec,
                    "--te-level", _te_level,
                    "--min-fibers", str(_min_fibers),
                    "--out-long", str(output.longtable),
                    "--out-loci", os.path.join(outdir, _loci_name),
                    "--tmp-dir", _tmp_dir,
                    "-l", log_path,
                ]
            else:
                cmd = [
                    "python", params.script,
                    "--mods", ",".join(
                        f"{s}:{indir}/{s}/{s}.{_mod_type}.bed.gz" for s in samples),
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


rule fiberseq_mods_sigtest:
    """Per-feature significance test on the mods long table (sigtest module).

    TE axis: features are subfamily/family names; significant groups then
    get a per-group top-loci shortlist (regionplot browsing BED) from the
    per-loci table via te_top_loci.py.
    """
    input:
        longtable = outdir + f"/{_long_name}",
        loci = ([outdir + "/" + _loci_name] if _is_te else []),
    output:
        tsv = outdir + "/{contrast}/significance_test.tsv",
        png = outdir + "/{contrast}/significance_bars.png",
        sig_tsv = outdir + "/{contrast}/significant_objects.tsv",
        sig_bed = outdir + "/{contrast}/significant_objects.bed",
        rest_tsv = outdir + "/{contrast}/nonsignificant_objects.tsv",
        top_loci_tsv = ([outdir + "/{contrast}/top_loci.tsv"] if _is_te else []),
        top_loci_bed = ([outdir + "/{contrast}/top_loci.bed"] if _is_te else []),
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
        top_loci_script = ROOT_DIR + "/modules/fiberseq/mods/bin/te_top_loci.py",
        group_spec = _group_samples_spec(),
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
            # TE axis: rank loci inside significant groups (browsing BED —
            # group-level p is the statistical assertion, loci are ranked)
            top_loci_cmd = []
            if _is_te:
                # group spec (params.group_spec = 'gA:s1+s2,gB:s3'), restricted
                # to the two groups sigtest compared (Params groups or first two)
                grp_pair = params.p.get("groups") or [
                    tok.split(":", 1)[0] for tok in params.group_spec.split(",")]
                by_group = {
                    tok.split(":", 1)[0]: tok.split(":", 1)[1]
                    for tok in params.group_spec.split(",")}
                grp_spec = ",".join(f"{g}:{by_group[g]}" for g in grp_pair[:2])
                top_loci_cmd = [
                    "python", params.top_loci_script,
                    "--sig", str(output.tsv),
                    "--loci", str(input.loci[0]),
                    "--groups", grp_spec,
                    "--value-suffix", f"{_mod_type}_per_fiber",
                    "--top-n", str(params.p.get("loci_top_n", 5)),
                    "--min-fibers", str(_min_fibers),
                    "--out-tsv", str(output.top_loci_tsv),
                    "--out-bed", str(output.top_loci_bed),
                ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in sig_cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in split_cmd) + "\n")
                if top_loci_cmd:
                    f.write(" ".join(shlex.quote(str(x)) for x in top_loci_cmd) + "\n")
                f.write(f'echo "fiberseq_mods_sigtest for {wildcards.contrast} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_mods_sigtest for {wildcards.contrast}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_mods_sigtest for {wildcards.contrast}: {e}, you can check the log file {log_path} for more details.\n")


rule fiberseq_mods_result:
    """Panel: aggregate mods outputs as a dependency point for use rule."""
    input:
        longtable = outdir + f"/{_long_name}",
        tsv = expand(outdir + "/{contrast}/significance_test.tsv",
                     contrast=sigtest_contrasts),
        sig_bed = expand(outdir + "/{contrast}/significant_objects.bed",
                         contrast=sigtest_contrasts),
        top_loci = (expand(outdir + "/{contrast}/top_loci.bed",
                           contrast=sigtest_contrasts) if _is_te else []),
