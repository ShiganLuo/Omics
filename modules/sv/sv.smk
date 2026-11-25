
include: "../common/common.smk"

import json  # noqa: E402  (common.smk 未导入而本模块需要的包)

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
logdir_combine = config.get("logdir_combine", logdir)
comparisons = config.get("comparisons", {})
comparison_names = list(comparisons)
vcf_dir = config.get("vcf_dir", "")
vcf_suffix = config.get("vcf_suffix", ".sv.vcf.gz")

# 指令派生名: parse 期格式列表(output 指令无法调用函数;str 收成单元素 list)
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
# gene_model 线程数的 Params 兜底(resources.<rule>.threads 优先,见 rule_res)
_threads_sv_gene_model = config.get("Params", {}).get("sv", {}).get("gene_model", {}).get("threads", 1)


def sv_exp_specific_inputs(wildcards):
    """sv_exp_specific file inputs — declared so the DAG validates them up front.

    Per-sample SV VCFs resolve under vcf_dir as
    ``{vcf_dir}/{sample_id}/{sample_id}{vcf_suffix}``.
    """
    c = comparisons[wildcards.comparison]
    return {
        "control_vcfs": [f"{vcf_dir}/{s}/{s}{vcf_suffix}" for s in c["control_samples"]],
        "experiment_vcfs": [f"{vcf_dir}/{s}/{s}{vcf_suffix}" for s in c["experimental_samples"]],
        "fasta": _require(config.get("genome", {}).get("fasta"), "genome.fasta", "sv"),
    }


rule sv_exp_specific:
    """Run experiment-specific SV analysis: merge, extract, annotate, plot, TE enrichment."""
    input:
        unpack(sv_exp_specific_inputs),
    output:
        merged_vcf = outdir + "/{comparison}/{comparison}_merged.vcf",
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        annotated_tab = outdir + "/{comparison}/{comparison}_annotated.tab",
        repeatmasker_done = outdir + "/{comparison}/repeatmasker/repeatmasker.done",
        statistics_done = outdir + "/{comparison}/statistics/statistics.done",
    log:
        logdir + "/{comparison}/sv_exp_specific.log"
    threads: rule_res("sv_exp_specific").get("threads", 1)
    benchmark: logdir + "/{comparison}/sv_exp_specific.benchmark.tsv",
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
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_exp_specific", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sv_exp_specific for comparison {wildcards.comparison} at {current_time}")
            sample_outdir = os.path.dirname(str(output.merged_vcf))
            os.makedirs(sample_outdir, exist_ok=True)
            outprefix = os.path.join(sample_outdir, wildcards.comparison)
            script = os.path.join(RECEIPT_DIR, f"sv_exp_specific_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"sv_exp_specific_{wildcards.comparison}.sh")
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
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "sv_exp_specific for {wildcards.comparison} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                # completion markers: repeatmasker/ and statistics/ are NOT
                # directory() outputs, so the engine's artifacts inside survive
                # reruns; the markers track the dirs in the DAG instead
                with open(str(output.repeatmasker_done), "w") as f:
                    f.write("ok\n")
                with open(str(output.statistics_done), "w") as f:
                    f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"Error occurred during sv_exp_specific for comparison {wildcards.comparison}: {e}\n")
            raise RuntimeError(f"Error occurred during sv_exp_specific for comparison {wildcards.comparison}: {e}, you can check the log file {log_path} for more details.\n")

rule sv_exp_circos:
    """Generate Circos plot for experiment-specific SVs."""
    input:
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "sv"),
    output:
        circos_done = outdir + "/{comparison}/circos/circos.done",
    log:
        logdir + "/{comparison}/sv_exp_circos.log"
    threads: rule_res("sv_exp_circos").get("threads", 1)
    benchmark: logdir + "/{comparison}/sv_exp_circos.benchmark.tsv",
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
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_exp_circos", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sv_exp_circos for comparison {wildcards.comparison} at {current_time}")
            circos_outdir = os.path.dirname(str(output.circos_done))
            os.makedirs(circos_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"sv_exp_circos_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(circos_outdir), f"sv_exp_circos_{wildcards.comparison}.sh")
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
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "sv_exp_circos for {wildcards.comparison} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                with open(str(output.circos_done), "w") as f:
                    f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"Error occurred during sv_exp_circos for comparison {wildcards.comparison}: {e}\n")
            raise RuntimeError(f"Error occurred during sv_exp_circos for comparison {wildcards.comparison}: {e}, you can check the log file {log_path} for more details.\n")

rule sv_gene_model:
    """Generate gene model plots with SV overlay for specified genes."""
    input:
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "sv"),
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        gtf = lambda wildcards: _require(config.get("genome", {}).get("gtf"), "genome.gtf", "sv"),
    output:
        gene_model_done = outdir + "/{comparison}/gene_model/gene_model.done",
    log:
        logdir + "/{comparison}/sv_gene_model.log"
    threads: rule_res("sv_gene_model").get("threads", _threads_sv_gene_model)
    benchmark: logdir + "/{comparison}/sv_gene_model.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_gene_model", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/utils/gene_model.py"),
        genes = _genes_sv_gene_model,   # ③ 转发
        formats = _formats_sv_gene_model,   # ③ 转发
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_gene_model")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_gene_model", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sv_gene_model for comparison {wildcards.comparison} at {current_time}")
            gene_outdir = os.path.dirname(str(output.gene_model_done))
            os.makedirs(gene_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"sv_gene_model_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(gene_outdir), f"sv_gene_model_{wildcards.comparison}.sh")
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
                "-j", str(threads),
            ] + fmt_args
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "sv_gene_model for {wildcards.comparison} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                with open(str(output.gene_model_done), "w") as f:
                    f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"Error occurred during sv_gene_model for comparison {wildcards.comparison}: {e}\n")
            raise RuntimeError(f"Error occurred during sv_gene_model for comparison {wildcards.comparison}: {e}, you can check the log file {log_path} for more details.\n")

rule sv_diff_analysis:
    """Run SV differential analysis across all comparison groups."""
    input:
        specific_vcfs = expand(outdir + "/{comparison}/{comparison}_only.vcf", comparison=comparison_names),
    output:
        diff_done = outdir + "/sv_diff_analysis/sv_diff_analysis.done",
    log:
        logdir_combine + "/sv/sv_diff_analysis.log"
    threads: rule_res("sv_diff_analysis").get("threads", 1)
    benchmark: logdir_combine + "/sv/sv_diff_analysis.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_diff_analysis", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/pbsv_sv_diff_analysis.py"),
        large_sv_threshold = config.get("Params", {}).get("sv", {}).get("diff_analysis", {}).get("large_sv_threshold", 10000),
        comparison_names = comparison_names,
        formats = _formats_sv_diff_analysis,   # ③ 转发
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_diff_analysis")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_diff_analysis", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sv_diff_analysis at {current_time}")
            diff_outdir = os.path.dirname(str(output.diff_done))
            os.makedirs(diff_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"sv_diff_analysis_{current_time}.sh")
            plain_script = os.path.join(diff_outdir, "sv_diff_analysis.sh")
            fmt_args = []
            for fmt in params.formats:
                fmt_args += ["-f", fmt]
            groups = {name: str(vcf) for name, vcf in zip(params.comparison_names, input.specific_vcfs)}
            cmd = [
                "python", params.script,
                "-o", diff_outdir,
                "-s", str(params.large_sv_threshold),
                "--groups-json", json.dumps(groups, ensure_ascii=False),
            ] + fmt_args
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write('echo "sv_diff_analysis was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                with open(str(output.diff_done), "w") as f:
                    f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"Error occurred during sv_diff_analysis: {e}\n")
            raise RuntimeError(f"Error occurred during sv_diff_analysis: {e}, you can check the log file {log_path} for more details.\n")

rule sv_oncoprint:
    """Generate OncoPrint plot from annotated SVs, COSMIC genes, and DESeq2 results."""
    input:
        annotated_tabs = expand(outdir + "/{comp}/{comp}_annotated.tab", comp=comparison_names),
        cosmic_file = lambda wildcards: _require(config.get("Params", {}).get("sv", {}).get("oncoprint", {}).get("cosmic_file"), "Params.sv.oncoprint.cosmic_file", "sv"),
    output:
        oncoprint_matrix = outdir + "/sv_oncoprint/oncoprint_matrix.csv",
        oncoprint_done = outdir + "/sv_oncoprint/oncoprint.done",
    log:
        logdir_combine + "/sv/sv_oncoprint.log"
    threads: rule_res("sv_oncoprint").get("threads", 1)
    benchmark: logdir_combine + "/sv/sv_oncoprint.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_oncoprint", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_OncoPrint.py"),
        deseq2_files = config.get("Params", {}).get("sv", {}).get("oncoprint", {}).get("deseq2_files", {}),
        comparison_names = comparison_names,
        formats = _formats_sv_oncoprint,   # ③ 转发
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_oncoprint")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_oncoprint", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sv_oncoprint at {current_time}")
            oncoprint_outdir = os.path.dirname(str(output.oncoprint_matrix))
            os.makedirs(oncoprint_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"sv_oncoprint_{current_time}.sh")
            plain_script = os.path.join(oncoprint_outdir, "sv_oncoprint.sh")
            fmt_args = []
            for fmt in params.formats:
                fmt_args += ["-f", fmt]
            groups = {name: str(tab) for name, tab in zip(params.comparison_names, input.annotated_tabs)}
            cmd = [
                "python", params.script,
                "--cosmic_file", input.cosmic_file,
                "-o", oncoprint_outdir,
                "--groups-json", json.dumps(groups, ensure_ascii=False),
                "--deseq2-json", json.dumps(params.deseq2_files, ensure_ascii=False),
            ] + fmt_args
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write('echo "sv_oncoprint was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                with open(str(output.oncoprint_done), "w") as f:
                    f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"Error occurred during sv_oncoprint: {e}\n")
            raise RuntimeError(f"Error occurred during sv_oncoprint: {e}, you can check the log file {log_path} for more details.\n")


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
        logdir + "/{comparison}/sv_mei.log"
    threads: rule_res("sv_mei").get("threads", 4)
    benchmark: logdir + "/{comparison}/sv_mei.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_mei", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    params:
        analysis = os.path.join(ROOT_DIR, "modules/sv/bin/mei_analysis.py"),
        mei = config.get("Params", {}).get("sv", {}).get("mei", {}) or {},
        te_lib = config.get("genome", {}).get("te_lib") or "",
        formats = _formats_sv_mei,   # ③ 转发
        ctrl_name = lambda wildcards: (comparisons.get(wildcards.comparison, {}) or {}).get("control_group_name") or "",
        exp_name = lambda wildcards: (comparisons.get(wildcards.comparison, {}) or {}).get("experimental_group_name") or "",
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_mei")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_mei", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sv_mei for comparison {wildcards.comparison} at {current_time}")
            mei_out = os.path.dirname(str(output.mei_done))
            os.makedirs(mei_out, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"sv_mei_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(mei_out), f"sv_mei_{wildcards.comparison}.sh")
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
                "--threads", str(threads),
                "--title", wildcards.comparison,
            ]
            if mei["dfam_h5"]:
                for h5 in mei["dfam_h5"]:
                    cmd += ["--dfam-h5", h5]
                cmd += ["--rm-species", mei.get("rm_species") or "",
                        "--dfam-mode", mei.get("dfam_mode", "export")]
            elif params.te_lib:
                cmd += ["--te-lib", params.te_lib]
            if mei.get("force_rerun"):
                cmd += ["--force-rerun"]
            # TE bucket list (Params.sv.mei.te_buckets), JSON transport
            cmd += ["--te-buckets",
                    json.dumps(mei.get("te_buckets"), ensure_ascii=False)]
            for flag, val in (("--ctrl-name", params.ctrl_name),
                              ("--exp-name", params.exp_name),
                              ("--fasta", getattr(input, "fasta", "")),
                              ("--gtf", getattr(input, "gtf", "")),
                              ("--te-bed", getattr(input, "te_bed", ""))):
                if val:
                    cmd += [flag, str(val)]
            # image formats: one -f per requested format (mei figures + report
            # are written per format by mei_analysis -> mei_plot)
            for fmt in params.formats:
                cmd += ["-f", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "sv_mei for {wildcards.comparison} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                # completion marker: mei/ is NOT a directory() output, so the
                # script-level caches and artifacts inside survive reruns
                with open(str(output.mei_done), "w") as f:
                    f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"Error occurred during sv_mei for comparison {wildcards.comparison}: {e}\n")
            raise RuntimeError(f"Error occurred during sv_mei for comparison {wildcards.comparison}: {e}, you can check the log file {log_path} for more details.\n")


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
        logdir + "/{comparison}/sv_tra_te.log"
    threads: rule_res("sv_tra_te").get("threads", 4)
    benchmark: logdir + "/{comparison}/sv_tra_te.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_tra_te", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    params:
        analysis = os.path.join(ROOT_DIR, "modules/sv/bin/utils/tra_te_analysis.py"),
        tra_te = config.get("Params", {}).get("sv", {}).get("tra_te", {}) or {},
        samtools_sif = sif("../samtools/samtools.yaml"),
        blast_sif = sif("sv.yaml"),
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml", "sv_tra_te")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_tra_te", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sv_tra_te for comparison {wildcards.comparison} at {current_time}")
            tra_out = os.path.dirname(str(output.tra_te_done))
            os.makedirs(tra_out, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"sv_tra_te_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(tra_out), f"sv_tra_te_{wildcards.comparison}.sh")
            tt = dict(params.tra_te)
            cmd = [
                "python", params.analysis,
                "--vcf", str(input.merged_vcf),
                "--te-bed", str(input.te_bed),
                "--fasta", str(input.fasta),
                "-o", tra_out,
                "--threads", str(threads),
                "--samtools-sif", params.samtools_sif,
                "--blast-sif", params.blast_sif,
                "--n-perm", str(tt.get("n_perm", 1000)),
            ]
            fai = getattr(input, "fai", "")
            if fai:
                cmd += ["--fai", str(fai)]
            if tt.get("seed") is not None:
                cmd += ["--seed", str(tt.get("seed", 42))]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "sv_tra_te for {wildcards.comparison} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
                # completion marker (tra_te/ is not a directory() output)
                with open(str(output.tra_te_done), "w") as f:
                    f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"Error occurred during sv_tra_te for comparison {wildcards.comparison}: {e}\n")
            raise RuntimeError(f"Error occurred during sv_tra_te for comparison {wildcards.comparison}: {e}, you can check the log file {log_path} for more details.\n")

rule sv_result:
    """Aggregate all SV analysis outputs as dependency endpoint."""
    input:
        annotated_tabs = expand(outdir + "/{comp}/{comp}_annotated.tab", comp=comparison_names),
        repeatmasker_done = expand(outdir + "/{comp}/repeatmasker/repeatmasker.done", comp=comparison_names),
        statistics_done = expand(outdir + "/{comp}/statistics/statistics.done", comp=comparison_names),
        circos_done = expand(outdir + "/{comp}/circos/circos.done", comp=comparison_names),
        gene_model_done = expand(outdir + "/{comp}/gene_model/gene_model.done", comp=comparison_names),
        mei_done = expand(outdir + "/{comp}/mei/mei.done", comp=comparison_names),
        tra_te_done = expand(outdir + "/{comp}/tra_te/tra_te.done", comp=comparison_names),
        diff_done = outdir + "/sv_diff_analysis/sv_diff_analysis.done",
        oncoprint_done = outdir + "/sv_oncoprint/oncoprint.done",
