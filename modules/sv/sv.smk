from snakemake.logging import logger
include: "../common/common.smk"
import time
import os

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# comparisons: DESeq2-style group_pairs dict keyed by comparison name
# ("{control_group_name}_vs_{experimental_group_name}"), each entry:
#   {control_group_name, experimental_group_name,
#    control_samples: [...], experimental_samples: [...]}
# SV VCF paths are resolved per sample from vcf_dir/vcf_suffix (pbsv layout);
# no file paths live in the config.
comparisons = config.get("comparisons", {})
comparison_names = list(comparisons)
vcf_dir = config.get("vcf_dir", "")
vcf_suffix = config.get("vcf_suffix", ".sv.vcf.gz")

def _sample_vcfs(sample_ids):
    """Resolve sample IDs to their SV VCF paths under vcf_dir."""
    return [f"{vcf_dir}/{s}/{s}{vcf_suffix}" for s in sample_ids]

# Global parameters
image_formats = config.get("Params", {}).get("image_formats", ["png"])
genes = config.get("Params", {}).get("genes", [])

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
        # TE enrichment + per-SV stats are produced by the same script run:
        # run_exp_enricher writes repeatmasker/ (bed2fa + RepeatMasker +
        # subfamily enrichment test), run_vcf_analysis writes plot/ and
        # table/. Declaring them keeps "output drives the workflow" intact —
        # previously they were undeclared side products that node.py never
        # registered, so nothing requested them.
        repeatmasker_dir = directory(outdir + "/{comparison}/repeatmasker"),
        plot_dir = directory(outdir + "/{comparison}/plot"),
        table_dir = directory(outdir + "/{comparison}/table"),
    log:
        logdir + "/{comparison}/sv_exp_specific.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/exp_specific.py"),
        vep_cache = config.get("Params", {}).get("exp_specific", {}).get("vep_cache", "~/.vep"),
        species = config.get("Params", {}).get("exp_specific", {}).get("species", "mus_musculus"),
        assembly = config.get("Params", {}).get("exp_specific", {}).get("assembly", "GRCm39"),
        custom_gtf = config.get("Params", {}).get("exp_specific", {}).get("custom_gtf", ""),
        custom_fasta = config.get("Params", {}).get("exp_specific", {}).get("custom_fasta", ""),
        te_lib = config.get("Params", {}).get("exp_specific", {}).get("te_lib", ""),
        dist = config.get("Params", {}).get("exp_specific", {}).get("dist", 500),
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_exp_specific",log_file = log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            outprefix = f"{outdir}/{wildcards.comparison}/{wildcards.comparison}"
            script_path = os.path.join(outdir, f"{wildcards.comparison}/sv_exp_specific_{current_time}.sh")
            rule_logger.info(f"Running sv_exp_specific for comparison: {wildcards.comparison}")
            fmt_args = []
            for fmt in image_formats:
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
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_exp_specific failed for {wildcards.comparison}: {e}")
            raise RuntimeError(f"sv_exp_specific failed for {wildcards.comparison}: {e}")

rule sv_exp_enrichment:
    """Run GO/KEGG enrichment and hotspot analysis on annotated SVs."""
    input:
        annotated_tab = outdir + "/{comparison}/{comparison}_annotated.tab",
        gtf = lambda wildcards: config.get("genome", {}).get("gtf") or [],
        fai = lambda wildcards: config.get("genome", {}).get("fai") or [],
    output:
        enrichment_dir = directory(outdir + "/{comparison}/enrichment"),
    log:
        logdir + "/{comparison}/sv_exp_enrichment.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_enrichment.py"),
        # exp_specific.species carries the assembly name (e.g.
        # T2T_mhaESC_v1.5); map it onto an Enrichr organism for gseapy
        organism = (
            "mouse" if "mouse" in str(config.get("Params", {}).get("exp_specific", {}).get("species", "mouse")).lower()
            or "mhESC".lower() in str(config.get("Params", {}).get("exp_specific", {}).get("species", "")).lower()
            or "GRCm".lower() in str(config.get("Params", {}).get("exp_specific", {}).get("species", "")).lower()
            else "human"
        ),
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_exp_enrichment", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            enrichment_outdir = f"{outdir}/{wildcards.comparison}/enrichment"
            script_path = os.path.join(outdir, f"{wildcards.comparison}/sv_exp_enrichment_{current_time}.sh")
            rule_logger.info(f"Running sv_exp_enrichment for comparison: {wildcards.comparison}")
            fmt_args = []
            for fmt in image_formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-a", input.annotated_tab,
                "-o", enrichment_outdir,
                "--organism", params.organism,
            ] + fmt_args
            # gtf/fai are optional (hotspot analysis); only pass when configured
            if input.gtf:
                cmd += ["-g", input.gtf]
            if input.fai:
                cmd += ["-fi", input.fai]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write("echo 'sv_exp_enrichment completed successfully.'\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_exp_enrichment failed for {wildcards.comparison}: {e}")
            raise RuntimeError(f"sv_exp_enrichment failed for {wildcards.comparison}: {e}")

rule sv_exp_circos:
    """Generate Circos plot for experiment-specific SVs."""
    input:
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "sv"),
    output:
        circos_dir = directory(outdir + "/{comparison}/circos"),
    log:
        logdir + "/{comparison}/sv_exp_circos.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_circos.py"),
        genome = config.get("Params", {}).get("circos", {}).get("genome", "mm39"),
        cytoband = config.get("Params", {}).get("circos", {}).get("cytoband", ""),
        ins_bin_size = config.get("Params", {}).get("circos", {}).get("ins_bin_size", 100000),
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_exp_circos", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            circos_outdir = f"{outdir}/{wildcards.comparison}/circos"
            script_path = os.path.join(outdir, f"{wildcards.comparison}/sv_exp_circos_{current_time}.sh")
            rule_logger.info(f"Running sv_exp_circos for comparison: {wildcards.comparison}")
            fmt_args = []
            for fmt in image_formats:
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
            shell(f"bash {script_path} >> {log_path} 2>&1")
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
        logdir + "/{comparison}/sv_gene_model.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/utils/gene_model.py"),
        threads = config.get("Params", {}).get("gene_model", {}).get("threads", 1),
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_gene_model", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            gene_outdir = f"{outdir}/{wildcards.comparison}/gene_model"
            script_path = os.path.join(outdir, f"{wildcards.comparison}/sv_gene_model_{current_time}.sh")
            rule_logger.info(f"Running sv_gene_model for comparison: {wildcards.comparison}")
            gene_args = []
            for g in genes:
                gene_args += ["-g", g]
            fmt_args = []
            for fmt in image_formats:
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
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_gene_model failed for {wildcards.comparison}: {e}")
            raise RuntimeError(f"sv_gene_model failed for {wildcards.comparison}: {e}")

# ========================
# Aggregation rules (cross-comparison)
# ========================

rule sv_diff_analysis:
    """Run SV differential analysis across all comparison groups."""
    input:
        specific_vcfs = expand(outdir + "/{comp}/{comp}_only.vcf", comp=comparison_names),
    output:
        diff_dir = directory(outdir + "/sv_diff_analysis"),
    log:
        logdir + "/all/sv_diff_analysis.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/pbsv_sv_diff_analysis.py"),
        large_sv_threshold = config.get("Params", {}).get("diff_analysis", {}).get("large_sv_threshold", 10000),
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_diff_analysis", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            diff_outdir = f"{outdir}/sv_diff_analysis"
            script_path = os.path.join(outdir, f"sv_diff_analysis_{current_time}.sh")
            rule_logger.info(f"Running sv_diff_analysis for all comparisons")
            group_args = []
            for name in comparison_names:
                vcf_path = f"{outdir}/{name}/{name}_only.vcf"
                group_args += ["-g", f"{name}:{vcf_path}"]
            fmt_args = []
            for fmt in image_formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-o", diff_outdir,
                "-s", str(params.large_sv_threshold),
            ] + group_args + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write("echo 'sv_diff_analysis completed successfully.'\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_diff_analysis failed: {e}")
            raise RuntimeError(f"sv_diff_analysis failed: {e}")

rule sv_oncoprint:
    """Generate OncoPrint plot from annotated SVs, COSMIC genes, and DESeq2 results."""
    input:
        annotated_tabs = expand(outdir + "/{comp}/{comp}_annotated.tab", comp=comparison_names),
        cosmic_file = lambda wildcards: _require(config.get("Params", {}).get("oncoprint", {}).get("cosmic_file"), "Params.sv.oncoprint.cosmic_file", "sv"),
    output:
        oncoprint_matrix = outdir + "/sv_oncoprint/oncoprint_matrix.csv",
        oncoprint_dir = directory(outdir + "/sv_oncoprint"),
    log:
        logdir + "/all/sv_oncoprint.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_OncoPrint.py"),
        deseq2_files = config.get("Params", {}).get("oncoprint", {}).get("deseq2_files", {}),
    conda:
        "sv.yaml"
    container:
        sif("sv.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sv_oncoprint", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            oncoprint_outdir = f"{outdir}/sv_oncoprint"
            script_path = os.path.join(outdir, f"sv_oncoprint_{current_time}.sh")
            rule_logger.info(f"Running sv_oncoprint for all comparisons")
            group_args = []
            for name in comparison_names:
                tab_path = f"{outdir}/{name}/{name}_annotated.tab"
                group_args += ["-g", f"{name}:{tab_path}"]
            deseq2_args = []
            for cond, path in params.deseq2_files.items():
                deseq2_args += ["-d", f"{cond}:{path}"]
            fmt_args = []
            for fmt in image_formats:
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
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sv_oncoprint failed: {e}")
            raise RuntimeError(f"sv_oncoprint failed: {e}")

# ========================
# Result aggregation rule
# ========================

rule sv_result:
    """Aggregate all SV analysis outputs as dependency endpoint."""
    input:
        annotated_tabs = expand(outdir + "/{comp}/{comp}_annotated.tab", comp=comparison_names),
        repeatmasker_dirs = expand(outdir + "/{comp}/repeatmasker", comp=comparison_names),
        plot_dirs = expand(outdir + "/{comp}/plot", comp=comparison_names),
        table_dirs = expand(outdir + "/{comp}/table", comp=comparison_names),
        enrichment_dirs = expand(outdir + "/{comp}/enrichment", comp=comparison_names),
        circos_dirs = expand(outdir + "/{comp}/circos", comp=comparison_names),
        gene_model_dirs = expand(outdir + "/{comp}/gene_model", comp=comparison_names),
        diff_dir = outdir + "/sv_diff_analysis",
        oncoprint_dir = outdir + "/sv_oncoprint",
