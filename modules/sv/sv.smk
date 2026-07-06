from snakemake.logging import logger
include: "../common/common.smk"
import time
import os

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# comparisons: list of dicts, each with control_vcf, experiment_vcf, name
comparisons = config.get("comparisons", [])
comparison_names = [c["name"] for c in comparisons]

# Global parameters
image_formats = config.get("Params", {}).get("image_formats", ["png"])
genes = config.get("Params", {}).get("genes", [])

# ========================
# Per-comparison rules
# ========================

rule sv_exp_specific:
    """Run experiment-specific SV analysis: merge, extract, annotate, plot, TE enrichment."""
    input:
        control_vcf = lambda wc: next(c["control_vcf"] for c in comparisons if c["name"] == wc.comparison),
        experiment_vcf = lambda wc: next(c["experiment_vcf"] for c in comparisons if c["name"] == wc.comparison),
        fasta = config.get("genome", {}).get("fasta", ""),
    output:
        merged_vcf = outdir + "/{comparison}/{comparison}_merged.vcf",
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        annotated_tab = outdir + "/{comparison}/{comparison}_annotated.tab",
    log:
        logdir + "/{comparison}/sv_exp_specific.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/exp_specific.py"),
        vep_cache = config.get("Params", {}).get("exp_specific", {}).get("vep_cache", "~/.vep"),
        species = config.get("Params", {}).get("exp_specific", {}).get("species", "mus_musculus"),
        assembly = config.get("Params", {}).get("exp_specific", {}).get("assembly", "GRCm39"),
        dist = config.get("Params", {}).get("exp_specific", {}).get("dist", 500),
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            outprefix = f"{outdir}/{wildcards.comparison}/{wildcards.comparison}"
            script_path = os.path.join(outdir, f"{wildcards.comparison}/sv_exp_specific_{current_time}.sh")
            fmt_args = []
            for fmt in image_formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-c", input.control_vcf,
                "-e", input.experiment_vcf,
                "-o", outprefix,
                "-d", str(params.dist),
                "--vep_cache", params.vep_cache,
                "--species", params.species,
                "--assembly", params.assembly,
                "--annotate_format", "tab",
            ] + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"sv_exp_specific failed for {wildcards.comparison}: {e}\n")
            raise

rule sv_exp_enrichment:
    """Run GO/KEGG enrichment and hotspot analysis on annotated SVs."""
    input:
        annotated_tab = outdir + "/{comparison}/{comparison}_annotated.tab",
        gtf = config.get("genome", {}).get("gtf", ""),
        fai = config.get("genome", {}).get("fai", ""),
    output:
        enrichment_dir = directory(outdir + "/{comparison}/enrichment"),
    log:
        logdir + "/{comparison}/sv_exp_enrichment.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_enrichment.py"),
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            enrichment_outdir = f"{outdir}/{wildcards.comparison}/enrichment"
            script_path = os.path.join(outdir, f"{wildcards.comparison}/sv_exp_enrichment_{current_time}.sh")
            fmt_args = []
            for fmt in image_formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-a", input.annotated_tab,
                "-o", enrichment_outdir,
                "-g", input.gtf,
                "-fi", input.fai,
            ] + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"sv_exp_enrichment failed for {wildcards.comparison}: {e}\n")
            raise

rule sv_exp_circos:
    """Generate Circos plot for experiment-specific SVs."""
    input:
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        fasta = config.get("genome", {}).get("fasta", ""),
    output:
        circos_dir = directory(outdir + "/{comparison}/circos"),
    log:
        logdir + "/{comparison}/sv_exp_circos.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_circos.py"),
        genome = config.get("Params", {}).get("circos", {}).get("genome", "mm39"),
        cytoband = config.get("Params", {}).get("circos", {}).get("cytoband", ""),
        ins_bin_size = config.get("Params", {}).get("circos", {}).get("ins_bin_size", 100000),
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            circos_outdir = f"{outdir}/{wildcards.comparison}/circos"
            script_path = os.path.join(outdir, f"{wildcards.comparison}/sv_exp_circos_{current_time}.sh")
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
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"sv_exp_circos failed for {wildcards.comparison}: {e}\n")
            raise

rule sv_gene_model:
    """Generate gene model plots with SV overlay for specified genes."""
    input:
        fasta = config.get("genome", {}).get("fasta", ""),
        specific_vcf = outdir + "/{comparison}/{comparison}_only.vcf",
        gtf = config.get("genome", {}).get("gtf", ""),
    output:
        gene_dir = directory(outdir + "/{comparison}/gene_model"),
    log:
        logdir + "/{comparison}/sv_gene_model.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/utils/gene_model.py"),
        threads = config.get("Params", {}).get("gene_model", {}).get("threads", 1),
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            gene_outdir = f"{outdir}/{wildcards.comparison}/gene_model"
            script_path = os.path.join(outdir, f"{wildcards.comparison}/sv_gene_model_{current_time}.sh")
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
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"sv_gene_model failed for {wildcards.comparison}: {e}\n")
            raise

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
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            diff_outdir = f"{outdir}/sv_diff_analysis"
            script_path = os.path.join(outdir, f"sv_diff_analysis_{current_time}.sh")
            group_args = []
            for c in comparisons:
                vcf_path = f"{outdir}/{c['name']}/{c['name']}_only.vcf"
                group_args += ["-g", f"{c['name']}:{vcf_path}"]
            fmt_args = []
            for fmt in image_formats:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-o", diff_outdir,
                "-s", str(params.large_sv_threshold),
            ] + group_args + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"sv_diff_analysis failed: {e}\n")
            raise

rule sv_oncoprint:
    """Generate OncoPrint plot from annotated SVs, COSMIC genes, and DESeq2 results."""
    input:
        annotated_tabs = expand(outdir + "/{comp}/{comp}_annotated.tab", comp=comparison_names),
        cosmic_file = config.get("Params", {}).get("oncoprint", {}).get("cosmic_file", ""),
    output:
        oncoprint_matrix = outdir + "/sv_oncoprint/oncoprint_matrix.csv",
        oncoprint_dir = directory(outdir + "/sv_oncoprint"),
    log:
        logdir + "/all/sv_oncoprint.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_OncoPrint.py"),
        deseq2_files = config.get("Params", {}).get("oncoprint", {}).get("deseq2_files", {}),
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            oncoprint_outdir = f"{outdir}/sv_oncoprint"
            script_path = os.path.join(outdir, f"sv_oncoprint_{current_time}.sh")
            group_args = []
            for c in comparisons:
                tab_path = f"{outdir}/{c['name']}/{c['name']}_annotated.tab"
                group_args += ["-g", f"{c['name']}:{tab_path}"]
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
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"sv_oncoprint failed: {e}\n")
            raise

# ========================
# Result aggregation rule
# ========================

rule sv_result:
    """Aggregate all SV analysis outputs as dependency endpoint."""
    input:
        annotated_tabs = expand(outdir + "/{comp}/{comp}_annotated.tab", comp=comparison_names),
        enrichment_dirs = expand(outdir + "/{comp}/enrichment", comp=comparison_names),
        circos_dirs = expand(outdir + "/{comp}/circos", comp=comparison_names),
        gene_model_dirs = expand(outdir + "/{comp}/gene_model", comp=comparison_names),
        diff_dir = outdir + "/sv_diff_analysis",
        oncoprint_dir = outdir + "/sv_oncoprint",
