
include: "../common/common.smk"

import json

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
comparisons = config.get("comparisons", {})
comparison_names = list(comparisons)
vcf_dir = config.get("vcf_dir", "")
vcf_suffix = config.get("vcf_suffix", ".sv.vcf.gz")

def _sample_vcfs(sample_ids):
    """Resolve sample IDs to their SV VCF paths under vcf_dir."""
    return [f"{vcf_dir}/{s}/{s}{vcf_suffix}" for s in sample_ids]

_formats_sv_exp_specific = config.get("Params", {}).get("sv", {}).get("sv_exp_specific", {}).get("formats") or ["png"]
if isinstance(_formats_sv_exp_specific, str):
    _formats_sv_exp_specific = [_formats_sv_exp_specific]
_formats_sv_exp_circos = config.get("Params", {}).get("sv", {}).get("sv_exp_circos", {}).get("formats") or ["png"]
if isinstance(_formats_sv_exp_circos, str):
    _formats_sv_exp_circos = [_formats_sv_exp_circos]
_formats_sv_gene_model = config.get("Params", {}).get("sv", {}).get("sv_gene_model", {}).get("formats") or ["png"]
if isinstance(_formats_sv_gene_model, str):
    _formats_sv_gene_model = [_formats_sv_gene_model]
_formats_sv_diff_analysis = config.get("Params", {}).get("sv", {}).get("sv_diff_analysis", {}).get("formats") or ["png"]
if isinstance(_formats_sv_diff_analysis, str):
    _formats_sv_diff_analysis = [_formats_sv_diff_analysis]
_formats_sv_oncoprint = config.get("Params", {}).get("sv", {}).get("sv_oncoprint", {}).get("formats") or ["png"]
if isinstance(_formats_sv_oncoprint, str):
    _formats_sv_oncoprint = [_formats_sv_oncoprint]
_formats_sv_mei = config.get("Params", {}).get("sv", {}).get("sv_mei", {}).get("formats") or ["png"]
if isinstance(_formats_sv_mei, str):
    _formats_sv_mei = [_formats_sv_mei]

# genes 同样 rule 级(Params.sv.sv_gene_model.genes),仅 sv_gene_model 消费。
_genes_sv_gene_model = config.get("Params", {}).get("sv", {}).get("sv_gene_model", {}).get("genes", [])

# ========================
# Per-comparison rules
# ========================

rule sv_exp_specific:
    """Run experiment-specific SV analysis: merge, extract, annotate, plot, TE enrichment."""
    input:
        control_vcfs = lambda wc: _sample_vcfs(comparisons[wc.comparison]["control_samples"]),
        experiment_vcfs = lambda wc: _sample_vcfs(comparisons[wc.comparison]["experimental_samples"]),
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "sv"),
    output:
        merged_vcf = outdir + "/{comparison}/{comparison}_merged.vcf",
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        annotated_tab = outdir + "/{comparison}/{comparison}_annotated.tab",
        repeatmasker_dir = directory(outdir + "/{comparison}/repeatmasker"),
        statistics_dir = directory(outdir + "/{comparison}/statistics"),
    log:
        logdir + "/sv/{comparison}/sv_exp_specific.log"
    benchmark: logdir + "/sv/{comparison}/sv_exp_specific.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_exp_specific", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/exp_specific.py"),
        vep_cache = config.get("Params", {}).get("sv", {}).get("exp_specific", {}).get("vep_cache", "~/.vep"),
        species = config.get("Params", {}).get("sv", {}).get("exp_specific", {}).get("species", "mus_musculus"),
        assembly = config.get("Params", {}).get("sv", {}).get("exp_specific", {}).get("assembly", "GRCm39"),
        custom_gtf = config.get("genome", {}).get("custom_gtf", ""),
        custom_fasta = config.get("genome", {}).get("custom_fasta", ""),
        te_lib = config.get("genome", {}).get("te_lib", ""),
        dist = config.get("Params", {}).get("sv", {}).get("exp_specific", {}).get("dist", 500),
        formats = _formats_sv_exp_specific,   # ③ 转发
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_exp_specific")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_exp_specific",log_file = log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            outprefix = f"{outdir}/{wildcards.comparison}/{wildcards.comparison}"
            script_path = os.path.join(RECEIPT_DIR, f"sv_exp_specific_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.comparison}/sv_exp_specific.sh")
            rule_logger.info(f"Running sv_exp_specific for comparison: {wildcards.comparison}")
            fmt_args = []
            for fmt in params.formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-c", *list(input.control_vcfs),
                "-e", *list(input.experiment_vcfs),
                "-o", outprefix,
                "-d", str(params.dist),
                "--vep_cache", params.vep_cache,
                "--species", params.species,
                "--assembly", params.assembly,
                "--annotate_format", "tab",
            ] + fmt_args
            if params.custom_gtf:
                cmd += ["--custom_gtf", params.custom_gtf]
                if params.custom_fasta:
                    cmd += ["--custom_fasta", params.custom_fasta]
            if params.te_lib:
                cmd += ["--te_lib", params.te_lib]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write("echo 'sv_exp_specific completed successfully.'\n")
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_exp_specific failed for {wildcards.comparison}: {e}")
            raise RuntimeError(f"sv_exp_specific failed for {wildcards.comparison}: {e}")
rule sv_exp_circos:
    """Generate Circos plot for experiment-specific SVs."""
    input:
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "sv"),
    output:
        circos_dir = directory(outdir + "/{comparison}/circos"),
    log:
        logdir + "/sv/{comparison}/sv_exp_circos.log"
    benchmark: logdir + "/sv/{comparison}/sv_exp_circos.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_exp_circos", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_circos.py"),
        genome = config.get("Params", {}).get("sv", {}).get("circos", {}).get("genome", "mm39"),
        cytoband = config.get("Params", {}).get("sv", {}).get("circos", {}).get("cytoband", ""),
        ins_bin_size = config.get("Params", {}).get("sv", {}).get("circos", {}).get("ins_bin_size", 100000),
        formats = _formats_sv_exp_circos,   # ③ 转发
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_exp_circos")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_exp_circos", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            circos_outdir = f"{outdir}/{wildcards.comparison}/circos"
            script_path = os.path.join(RECEIPT_DIR, f"sv_exp_circos_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.comparison}/sv_exp_circos.sh")
            rule_logger.info(f"Running sv_exp_circos for comparison: {wildcards.comparison}")
            fmt_args = []
            for fmt in params.formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "--vcf", input.specific_vcf,
                "--fasta", input.fasta,
                "--outdir", circos_outdir,
                "--genome", params.genome,
                "--ins_bin_size", str(params.ins_bin_size),
            ] + fmt_args
            if params.cytoband:
                cmd += ["--cytoband", params.cytoband]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write("echo 'sv_exp_circos completed successfully.'\n")
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_exp_circos failed for {wildcards.comparison}: {e}")
            raise RuntimeError(f"sv_exp_circos failed for {wildcards.comparison}: {e}")

rule sv_gene_model:
    """Generate gene model plots with SV overlay for specified genes."""
    input:
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "sv"),
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        gtf = lambda wildcards: _require(config.get("genome", {}).get("gtf"), "genome.gtf", "sv"),
    output:
        gene_dir = directory(outdir + "/{comparison}/gene_model"),
    log:
        logdir + "/sv/{comparison}/sv_gene_model.log"
    benchmark: logdir + "/sv/{comparison}/sv_gene_model.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_gene_model", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/utils/gene_model.py"),
        threads = config.get("Params", {}).get("sv", {}).get("gene_model", {}).get("threads", 1),
        genes = _genes_sv_gene_model,   # ③ 转发
        formats = _formats_sv_gene_model,   # ③ 转发
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_gene_model")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_gene_model", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            gene_outdir = f"{outdir}/{wildcards.comparison}/gene_model"
            script_path = os.path.join(RECEIPT_DIR, f"sv_gene_model_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.comparison}/sv_gene_model.sh")
            rule_logger.info(f"Running sv_gene_model for comparison: {wildcards.comparison}")
            gene_args = []
            for g in params.genes:
                gene_args += ["-g", g]
            fmt_args = []
            for fmt in params.formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-t", input.gtf,
            ] + gene_args + [
                "-m", input.specific_vcf,
                "-o", gene_outdir,
                "-j", str(params.threads),
            ] + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write("echo 'sv_gene_model completed successfully.'\n")
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_gene_model failed for {wildcards.comparison}: {e}")
            raise RuntimeError(f"sv_gene_model failed for {wildcards.comparison}: {e}")

        # ========================
        # Aggregation rules (cross-comparison)
        # ========================

rule sv_diff_analysis:
    """Run SV differential analysis across all comparison groups."""
    input:
        specific_vcfs = expand(outdir + "/{comparison}/{comparison}_only.vcf", comparison=comparison_names),
    output:
        diff_dir = directory(outdir + "/sv_diff_analysis"),
    log:
        logdir + "/sv/sv_diff_analysis.log"
    benchmark: logdir + "/sv/sv_diff_analysis.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_diff_analysis", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/pbsv_sv_diff_analysis.py"),
        large_sv_threshold = config.get("Params", {}).get("sv", {}).get("diff_analysis", {}).get("large_sv_threshold", 10000),
        formats = _formats_sv_diff_analysis,   # ③ 转发
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_diff_analysis")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_diff_analysis", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(RECEIPT_DIR, f"sv_diff_analysis_{current_time}.sh")
            plain_script = os.path.join(outdir, f"sv_diff_analysis.sh")
            rule_logger.info(f"Running sv_diff_analysis for all comparisons")
            group_args = []
            for name in comparison_names:
                vcf_path = f"{outdir}/{name}/{name}_only.vcf"
                group_args += ["-g", f"{name}:{vcf_path}"]
            fmt_args = []
            for fmt in params.formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-o", str(output.diff_dir),
                "-s", str(params.large_sv_threshold),
            ] + group_args + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write("echo 'sv_diff_analysis completed successfully.'\n")
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_diff_analysis failed: {e}")
            raise RuntimeError(f"sv_diff_analysis failed: {e}")

rule sv_oncoprint:
    """Generate OncoPrint plot from annotated SVs, COSMIC genes, and DESeq2 results."""
    input:
        annotated_tabs = expand(outdir + "/{comp}/{comp}_annotated.tab", comp=comparison_names),
        cosmic_file = lambda wildcards: _require(config.get("Params", {}).get("sv", {}).get("oncoprint", {}).get("cosmic_file"), "Params.sv.oncoprint.cosmic_file", "sv"),
    output:
        oncoprint_matrix = outdir + "/sv_oncoprint/oncoprint_matrix.csv",
        oncoprint_dir = directory(outdir + "/sv_oncoprint"),
    log:
        logdir + "/sv/oncoprint.log"
    benchmark: logdir + "/sv/oncoprint.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_oncoprint", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_OncoPrint.py"),
        deseq2_files = config.get("Params", {}).get("sv", {}).get("oncoprint", {}).get("deseq2_files", {}),
        formats = _formats_sv_oncoprint,   # ③ 转发
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_oncoprint")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_oncoprint", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            oncoprint_outdir = f"{outdir}/sv_oncoprint"
            script_path = os.path.join(RECEIPT_DIR, f"sv_oncoprint_{current_time}.sh")
            plain_script = os.path.join(outdir, f"sv_oncoprint.sh")
            rule_logger.info(f"Running sv_oncoprint for all comparisons")
            group_args = []
            for name in comparison_names:
                tab_path = f"{outdir}/{name}/{name}_annotated.tab"
                group_args += ["-g", f"{name}:{tab_path}"]
            deseq2_args = []
            for cond, path in params.deseq2_files.items():
                deseq2_args += ["-d", f"{cond}:{path}"]
            fmt_args = []
            for fmt in params.formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "--cosmic_file", input.cosmic_file,
                "-o", oncoprint_outdir,
            ] + group_args + deseq2_args + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write("echo 'sv_oncoprint completed successfully.'\n")
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_oncoprint failed: {e}")
            raise RuntimeError(f"sv_oncoprint failed: {e}")


def sv_mei_inputs(wildcards):
    """sv_mei file inputs — declared so the DAG validates them up front.

    Classifier source: Dfam FamDB h5 partition(s) (Params.sv.mei.dfam_h5,
    preferred) or the legacy te_lib consensus FASTA (genome.te_lib).
    Optional reference files join the inputs only when configured, so a
    missing file fails at DAG build instead of mid-execution.
    """
    mei = dict(config.get("Params", {}).get("sv", {}).get("mei", {}) or {})
    dfam_h5 = mei.get("dfam_h5") or []
    mei["dfam_h5"] = [dfam_h5] if isinstance(dfam_h5, str) else list(dfam_h5)
    te_lib = config.get("genome", {}).get("te_lib") or ""
    classifier = mei["dfam_h5"] or ([te_lib] if te_lib and os.path.exists(te_lib) else [])
    if not classifier:
        raise ValueError("sv_mei: set Params.sv.mei.dfam_h5 or genome.te_lib")
    files = {
        "merged_vcf": outdir + f"/{wildcards.comparison}/{wildcards.comparison}_merged.vcf",
        "classifier": classifier,
    }
    genome = config.get("genome", {})
    for name, key in (("fasta", "fasta"), ("gtf", "gtf"), ("te_bed", "repeat_bed")):
        if genome.get(key):
            files[name] = genome[key]
    return files


rule sv_mei:
    """MEI (mobile element insertion) analysis per comparison.

    Consumes the SURVIVOR-merged VCF from sv_exp_specific; classifies INS
    alleles by RepeatMasker (Dfam famdb export or the te_lib consensus),
    detects poly(A)/TSD fingerprints, annotates genomic context, and
    renders the 6-figure suite via the single entry
    modules/sv/bin/mei_analysis.py. TE insertion dating
    (te_age_analysis.py) is a separate standalone analysis, not run here.
    """
    input:
        unpack(sv_mei_inputs),
    output:
        mei_done = outdir + "/{comparison}/mei/mei.done",
    log:
        logdir + "/sv/{comparison}/sv_mei.log"
    benchmark: logdir + "/sv/{comparison}/sv_mei.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_mei", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    params:
        analysis = os.path.join(ROOT_DIR, "modules/sv/bin/mei_analysis.py"),
        mei = config.get("Params", {}).get("sv", {}).get("mei", {}) or {},
        te_lib = config.get("genome", {}).get("te_lib") or "",
        threads = lambda wildcards: rule_res("sv_mei").get("threads", 4),
        formats = _formats_sv_mei,   # ③ 转发
        ctrl_name = lambda wildcards: (comparisons.get(wildcards.comparison, {}) or {}).get("control_group_name") or "",
        exp_name = lambda wildcards: (comparisons.get(wildcards.comparison, {}) or {}).get("experimental_group_name") or "",
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_mei")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_mei", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            mei_out = os.path.join(str(outdir), wildcards.comparison, "mei")
            os.makedirs(mei_out, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"sv_mei_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(str(outdir), wildcards.comparison, "sv_mei.sh")
            mei = dict(params.mei)
            dfam_h5 = mei.get("dfam_h5") or []
            mei["dfam_h5"] = [dfam_h5] if isinstance(dfam_h5, str) else list(dfam_h5)
            cmd = [
                "python", params.analysis,
                "-i", str(input.merged_vcf),
                "-o", mei_out,
                "--min-len", str(mei.get("min_len", 100)),
                "--full-length-l1", str(mei.get("full_length", 6000)),
                "--full-length-family", mei.get("full_length_family", "L1"),
                "--threads", str(params.threads),
                "--title", wildcards.comparison,
            ]
            if mei["dfam_h5"]:
                for h5 in mei["dfam_h5"]:
                    cmd += ["--dfam-h5", h5]
                cmd += ["--rm-species", shlex.quote(mei.get("rm_species") or ""),
                        "--dfam-mode", mei.get("dfam_mode", "export")]
            elif params.te_lib:
                cmd += ["--te-lib", params.te_lib]
            if mei.get("force_rerun"):
                cmd += ["--force-rerun"]
            # TE bucket list (Params.sv.mei.te_buckets), JSON transport
            cmd += ["--te-buckets",
                    shlex.quote(json.dumps(mei.get("te_buckets"), ensure_ascii=False))]
            for flag, val in (("--ctrl-name", params.ctrl_name),
                              ("--exp-name", params.exp_name),
                              ("--fasta", getattr(input, "fasta", "")),
                              ("--gtf", getattr(input, "gtf", "")),
                              ("--te-bed", getattr(input, "te_bed", ""))):
                if val:
                    cmd += [flag, shlex.quote(str(val))]
            # image formats: one -f per requested format (mei figures +
            # report are written per format by mei_analysis -> mei_plot)
            for fmt in params.formats:
                cmd += ["-f", str(fmt)]
            # figures run INSIDE the entry (utils/mei_plot.py chained)
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "sv_mei completed successfully for {wildcards.comparison}"\n')
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            # completion marker: mei/ is NOT a directory() output, so the
            # script-level caches and artifacts inside survive reruns
            with open(str(output.mei_done), "w") as f:
                f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"sv_mei failed for {wildcards.comparison}: {e}")
            raise RuntimeError(f"sv_mei failed for {wildcards.comparison}: {e}")


def sv_tra_te_inputs(wildcards):
    """sv_tra_te file inputs — validated up front.

    Needs the SURVIVOR-merged VCF (same as sv_mei), the reference FASTA
    (samtools faidx TE-copy extraction), the TE BED (official col7
    family_id / col8 class_id used verbatim when present) and optionally
    its .fai (real chromosome lengths for the panorama figure).
    """
    files = {
        "merged_vcf": outdir + f"/{wildcards.comparison}/{wildcards.comparison}_merged.vcf",
    }
    genome = config.get("genome", {})
    for name, key in (("fasta", "fasta"), ("te_bed", "repeat_bed"),
                      ("fai", "fai")):
        if genome.get(key):
            files[name] = genome[key]
    if "fasta" not in files:
        raise ValueError("sv_tra_te: genome.fasta is required "
                         "(TE-copy sequence extraction via samtools faidx)")
    if "te_bed" not in files:
        raise ValueError("sv_tra_te: genome.repeat_bed is required "
                         "(breakpoint TE annotation)")
    return files


rule sv_tra_te:
    """TE-associated TRA analysis per comparison (T0-T3 evidence funnel).

    Consumes the SURVIVOR-merged VCF; extracts inter-chromosomal TRA
    events (mirror dedup), tests breakpoint-TE enrichment against a
    permutation null, pairs official RepeatMasker family ids against a
    shuffle expectation, and validates NAHR-grade homology by blastn of
    the two TE copies. Methods: modules/sv/README.md section 3.3 (TRA/TE-TE funnel).
    (bin/utils/tra_te_analysis.py). Figures (tra_te_figures) are part of
    that entry's stage 4.
    """
    input:
        unpack(sv_tra_te_inputs),
    output:
        tra_te_done = outdir + "/{comparison}/tra_te/tra_te.done",
    log:
        logdir + "/sv/{comparison}/sv_tra_te.log"
    benchmark: logdir + "/sv/{comparison}/sv_tra_te.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_tra_te", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    params:
        analysis = os.path.join(ROOT_DIR, "modules/sv/bin/utils/tra_te_analysis.py"),
        tra_te = config.get("Params", {}).get("sv", {}).get("tra_te", {}) or {},
        threads = lambda wildcards: rule_res("sv_tra_te").get("threads", 4),
        samtools_sif = sif("../samtools/samtools.yaml"),
        blast_sif = sif("sv.yaml"),
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_tra_te")
    run:
        import shlex
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_tra_te", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            tra_out = os.path.join(str(outdir), wildcards.comparison, "tra_te")
            os.makedirs(tra_out, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"sv_tra_te_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(str(outdir), wildcards.comparison, "sv_tra_te.sh")
            tt = dict(params.tra_te)
            cmd = [
                "python", params.analysis,
                "--vcf", str(input.merged_vcf),
                "--te-bed", str(input.te_bed),
                "--fasta", str(input.fasta),
                "-o", tra_out,
                "--threads", str(params.threads),
                "--samtools-sif", params.samtools_sif,
                "--blast-sif", params.blast_sif,
                "--n-perm", str(tt.get("n_perm", 1000)),
            ]
            fai = getattr(input, "fai", "")
            if fai:
                cmd += ["--fai", str(fai)]
            if tt.get("seed") is not None:
                cmd += ["--seed", str(tt.get("seed", 42))]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
                f.write(f'echo "sv_tra_te completed successfully for {wildcards.comparison}"\n')
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            # completion marker (tra_te/ is not a directory() output)
            with open(str(output.tra_te_done), "w") as f:
                f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"sv_tra_te failed for {wildcards.comparison}: {e}")
            raise RuntimeError(f"sv_tra_te failed for {wildcards.comparison}: {e}")

rule sv_result:
    """Aggregate all SV analysis outputs as dependency endpoint."""
    input:
        annotated_tabs = expand(outdir + "/{comp}/{comp}_annotated.tab", comp=comparison_names),
        repeatmasker_dirs = expand(outdir + "/{comp}/repeatmasker", comp=comparison_names),
        statistics_dirs = expand(outdir + "/{comp}/statistics", comp=comparison_names),
        enrichment_dirs = expand(outdir + "/{comp}/enrichment", comp=comparison_names),
        circos_dirs = expand(outdir + "/{comp}/circos", comp=comparison_names),
        gene_model_dirs = expand(outdir + "/{comp}/gene_model", comp=comparison_names),
        mei_done = expand(outdir + "/{comp}/mei/mei.done", comp=comparison_names),
        tra_te_done = expand(outdir + "/{comp}/tra_te/tra_te.done", comp=comparison_names),
        diff_dir = outdir + "/sv_diff_analysis",
        oncoprint_dir = outdir + "/sv_oncoprint",
