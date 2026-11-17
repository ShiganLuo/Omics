
include: "../../common/common.smk"


ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# comparisons: list of dicts, each with control_vcf, experiment_vcf, name
comparisons = config.get("comparisons", [])
comparison_names = [c["name"] for c in comparisons]

# Global parameters
# 图片格式(modules.md #9):键名 "formats"。参数精确到 rule:formats 按 rule
# 独立配置(Params.<rule_name>.formats,无 module 级共享回退;rule_param 见
# common.smk),每 rule 一个 parse 期变量。
_sv_params = config.get("Params", {})


def _rule_formats(rule_name):
    """Per-rule figure formats (parse 期解析,CLI -f 转发共用)。"""
    fmts = rule_param(_sv_params, rule_name, "formats", ["png"])
    return [fmts] if isinstance(fmts, str) else fmts


_formats_sv_exp_specific = _rule_formats("sv_exp_specific")
_formats_sv_exp_enrichment = _rule_formats("sv_exp_enrichment")
_formats_sv_exp_circos = _rule_formats("sv_exp_circos")
_formats_sv_gene_model = _rule_formats("sv_gene_model")
_formats_sv_diff_analysis = _rule_formats("sv_diff_analysis")
_formats_sv_oncoprint = _rule_formats("sv_oncoprint")

# genes 同样 rule 级(Params.sv_gene_model.genes),仅 sv_gene_model 消费。
_genes_sv_gene_model = rule_param(_sv_params, "sv_gene_model", "genes", [])

# ========================
# Per-comparison rules
# ========================

rule sv_exp_specific:
    """Run experiment-specific SV analysis: merge, extract, annotate, plot, TE enrichment."""
    input:
        control_vcf = lambda wc: next(c["control_vcf"] for c in comparisons if c["name"] == wc.comparison),
        experiment_vcf = lambda wc: next(c["experiment_vcf"] for c in comparisons if c["name"] == wc.comparison),
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"],
    output:
        merged_vcf = outdir + "/{genome}/{comparison}/{comparison}_merged.vcf",
        specific_vcf = outdir + "/{genome}/{comparison}/{comparison}_only.vcf",
        annotated_tab = outdir + "/{genome}/{comparison}/{comparison}_annotated.tab",
    log:
        logdir + "/{genome}/{comparison}/sv_exp_specific.log"
    benchmark: logdir + "/{genome}/{comparison}/sv_exp_specific.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_exp_specific", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/exp_specific.py"),
        vep_cache = config.get("Params", {}).get("exp_specific", {}).get("vep_cache", "~/.vep"),
        species = config.get("Params", {}).get("exp_specific", {}).get("species", "mus_musculus"),
        assembly = config.get("Params", {}).get("exp_specific", {}).get("assembly", "GRCm39"),
        dist = config.get("Params", {}).get("exp_specific", {}).get("dist", 500),
    conda:
        "../sv.yaml"
    container:
        sif("../sv.yaml", "sv_exp_specific")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            outprefix = f"{outdir}/{wildcards.genome}/{wildcards.comparison}/{wildcards.comparison}"
            script_path = os.path.join(RECEIPT_DIR, f"sv_exp_specific_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.comparison}/sv_exp_specific.sh")
            fmt_args = []
            for fmt in _formats_sv_exp_specific:
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
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule sv_exp_enrichment:
    """Run GO/KEGG enrichment and hotspot analysis on annotated SVs."""
    input:
        annotated_tab = outdir + "/{genome}/{comparison}/{comparison}_annotated.tab",
        gtf = lambda wildcards: config["genome"][wildcards.genome].get("gtf", ""),
        fai = lambda wildcards: config["genome"][wildcards.genome].get("fai", ""),
    output:
        enrichment_dir = directory(outdir + "/{genome}/{comparison}/enrichment"),
    log:
        logdir + "/{genome}/{comparison}/sv_exp_enrichment.log"
    benchmark: logdir + "/{genome}/{comparison}/sv_exp_enrichment.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_exp_enrichment", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_enrichment.py"),
    conda:
        "../sv.yaml"
    container:
        sif("../sv.yaml", "sv_exp_enrichment")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            enrichment_outdir = f"{outdir}/{wildcards.genome}/{wildcards.comparison}/enrichment"
            script_path = os.path.join(RECEIPT_DIR, f"sv_exp_enrichment_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.comparison}/sv_exp_enrichment.sh")
            fmt_args = []
            for fmt in _formats_sv_exp_enrichment:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-a", input.annotated_tab,
                "-o", enrichment_outdir,
                "-g", input.gtf,
                "-fi", input.fai,
            ] + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule sv_exp_circos:
    """Generate Circos plot for experiment-specific SVs."""
    input:
        specific_vcf = outdir + "/{genome}/{comparison}/{comparison}_only.vcf",
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"],
    output:
        circos_dir = directory(outdir + "/{genome}/{comparison}/circos"),
    log:
        logdir + "/{genome}/{comparison}/sv_exp_circos.log"
    benchmark: logdir + "/{genome}/{comparison}/sv_exp_circos.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_exp_circos", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_circos.py"),
        genome = config.get("Params", {}).get("circos", {}).get("genome", "mm39"),
        cytoband = config.get("Params", {}).get("circos", {}).get("cytoband", ""),
        ins_bin_size = config.get("Params", {}).get("circos", {}).get("ins_bin_size", 100000),
    conda:
        "../sv.yaml"
    container:
        sif("../sv.yaml", "sv_exp_circos")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            circos_outdir = f"{outdir}/{wildcards.genome}/{wildcards.comparison}/circos"
            script_path = os.path.join(RECEIPT_DIR, f"sv_exp_circos_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.comparison}/sv_exp_circos.sh")
            fmt_args = []
            for fmt in _formats_sv_exp_circos:
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
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule sv_gene_model:
    """Generate gene model plots with SV overlay for specified genes."""
    input:
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"],
        specific_vcf = outdir + "/{genome}/{comparison}/{comparison}_only.vcf",
        gtf = lambda wildcards: config["genome"][wildcards.genome].get("gtf", ""),
    output:
        gene_dir = directory(outdir + "/{genome}/{comparison}/gene_model"),
    log:
        logdir + "/{genome}/{comparison}/sv_gene_model.log"
    benchmark: logdir + "/{genome}/{comparison}/sv_gene_model.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_gene_model", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/utils/gene_model.py"),
        threads = config.get("Params", {}).get("gene_model", {}).get("threads", 1),
    conda:
        "../sv.yaml"
    container:
        sif("../sv.yaml", "sv_gene_model")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            gene_outdir = f"{outdir}/{wildcards.genome}/{wildcards.comparison}/gene_model"
            script_path = os.path.join(RECEIPT_DIR, f"sv_gene_model_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.comparison}/sv_gene_model.sh")
            gene_args = []
            for g in _genes_sv_gene_model:
                gene_args += ["-g", g]
            fmt_args = []
            for fmt in _formats_sv_gene_model:
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
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

        # ========================
        # Aggregation rules (cross-comparison)
        # ========================

rule sv_diff_analysis:
    """Run SV differential analysis across all comparison groups."""
    input:
        specific_vcfs = expand(outdir + "/{{genome}}/{comp}/{comp}_only.vcf", comp=comparison_names),
    output:
        diff_dir = directory(outdir + "/{genome}/sv_diff_analysis"),
    log:
        logdir + "/{genome}/all/sv_diff_analysis.log"
    benchmark: logdir + "/{genome}/all/sv_diff_analysis.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_diff_analysis", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/pbsv_sv_diff_analysis.py"),
        large_sv_threshold = config.get("Params", {}).get("diff_analysis", {}).get("large_sv_threshold", 10000),
    conda:
        "../sv.yaml"
    container:
        sif("../sv.yaml", "sv_diff_analysis")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            diff_outdir = f"{outdir}/{wildcards.genome}/sv_diff_analysis"
            script_path = os.path.join(RECEIPT_DIR, f"sv_diff_analysis_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.genome}/sv_diff_analysis.sh")
            group_args = []
            for c in comparisons:
                vcf_path = f"{outdir}/{wildcards.genome}/{c['name']}/{c['name']}_only.vcf"
                group_args += ["-g", f"{c['name']}:{vcf_path}"]
            fmt_args = []
            for fmt in _formats_sv_diff_analysis:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "-o", diff_outdir,
                "-s", str(params.large_sv_threshold),
            ] + group_args + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule sv_oncoprint:
    """Generate OncoPrint plot from annotated SVs, COSMIC genes, and DESeq2 results."""
    input:
        annotated_tabs = expand(outdir + "/{{genome}}/{comp}/{comp}_annotated.tab", comp=comparison_names),
        cosmic_file = config.get("Params", {}).get("oncoprint", {}).get("cosmic_file", ""),
    output:
        oncoprint_matrix = outdir + "/{genome}/sv_oncoprint/oncoprint_matrix.csv",
        oncoprint_dir = directory(outdir + "/{genome}/sv_oncoprint"),
    log:
        logdir + "/{genome}/all/sv_oncoprint.log"
    benchmark: logdir + "/{genome}/all/sv_oncoprint.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_oncoprint", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/sv/bin/run_OncoPrint.py"),
        deseq2_files = config.get("Params", {}).get("oncoprint", {}).get("deseq2_files", {}),
    conda:
        "../sv.yaml"
    container:
        sif("../sv.yaml", "sv_oncoprint")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            oncoprint_outdir = f"{outdir}/{wildcards.genome}/sv_oncoprint"
            script_path = os.path.join(RECEIPT_DIR, f"sv_oncoprint_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.genome}/sv_oncoprint.sh")
            group_args = []
            for c in comparisons:
                tab_path = f"{outdir}/{wildcards.genome}/{c['name']}/{c['name']}_annotated.tab"
                group_args += ["-g", f"{c['name']}:{tab_path}"]
            deseq2_args = []
            for cond, path in params.deseq2_files.items():
                deseq2_args += ["-d", f"{cond}:{path}"]
            fmt_args = []
            for fmt in _formats_sv_oncoprint:
                fmt_args += ["-f", fmt]
            cmd = [
                "python", params.script,
                "--cosmic_file", input.cosmic_file,
                "-o", oncoprint_outdir,
            ] + group_args + deseq2_args + fmt_args
            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

        # ========================
        # Result aggregation rule
        # ========================

rule sv_result:
    """Aggregate all SV analysis outputs as dependency endpoint."""
    input:
        annotated_tabs = expand(outdir + "/{{genome}}/{comp}/{comp}_annotated.tab", comp=comparison_names),
        enrichment_dirs = expand(outdir + "/{{genome}}/{comp}/enrichment", comp=comparison_names),
        circos_dirs = expand(outdir + "/{{genome}}/{comp}/circos", comp=comparison_names),
        gene_model_dirs = expand(outdir + "/{{genome}}/{comp}/gene_model", comp=comparison_names),
        diff_dir = outdir + "/{genome}/sv_diff_analysis",
        oncoprint_dir = outdir + "/{genome}/sv_oncoprint",

    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sv_result", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
