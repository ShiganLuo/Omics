include: "../../common/common.smk"

indir = config.get("indir", "data")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
genome_paired_samples = config.get("genome_paired_samples", {})
genome_single_samples = config.get("genome_single_samples", {})


_formats = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("formats") or ["png"]
if isinstance(_formats, str):
    _formats = [_formats]

_RASTER_FMT_ORDER = ("png", "jpg", "jpeg", "tif", "tiff")


def _producer_fmt(formats):
    """Raster-first extension among one producer's configured formats."""
    fmts = formats or ["png"]
    if isinstance(fmts, dict):
        fmts = fmts.get("default") or ["png"]
    if isinstance(fmts, str):
        fmts = [fmts]
    for fmt in _RASTER_FMT_ORDER:
        if fmt in fmts:
            return fmt
    return fmts[0]


_st_formats = config.get("Params", {}).get("StringTie", {}).get("generate_report", {}).get("formats") or ["png"]
if isinstance(_st_formats, str):
    _st_formats = [_st_formats]
_deseq2_formats = config.get("Params", {}).get("DESeq2", {}).get("generate_report", {}).get("formats") or ["png"]
if isinstance(_deseq2_formats, str):
    _deseq2_formats = [_deseq2_formats]
_arriba_formats = config.get("Params", {}).get("arriba", {}).get("generate_report", {}).get("formats") or ["png"]
if isinstance(_arriba_formats, str):
    _arriba_formats = [_arriba_formats]


REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "RNAseq_report", "bin", "generate_report.py")
FUSION_FIG_BASES = [
    "fig1_per_sample_counts",
    "fig2_fusion_type_distribution",
    "fig3_type_by_sample",
    "fig4_reading_frame",
    "fig5_recurrent_heatmap",
    "fig6_inframe_support",
]

def get_input_for_generate_report(wildcards):
    logger.debug(f"[get_input_for_generate_report] called with wildcards: {wildcards}")
    in_dict = {}
    gsva_fmt = config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("generate_report", {}).get("gsva_formats") or ["png"]
    if isinstance(gsva_fmt, str):
        gsva_fmt = [gsva_fmt]
    go_kegg_fmt = config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("generate_report", {}).get("go_kegg_formats") or ["png"]
    if isinstance(go_kegg_fmt, str):
        go_kegg_fmt = [go_kegg_fmt]
    gsea_fmt = config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("generate_report", {}).get("gsea_formats") or ["jpeg"]
    if isinstance(gsea_fmt, str):
        gsea_fmt = [gsea_fmt]
    in_dict["te_sample_summary"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_sample_summary.tsv"
    in_dict["te_group_summary"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_group_summary.tsv"
    in_dict["te_type_counts"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_counts.tsv"
    in_dict["te_group_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_group_stacked." + _producer_fmt(_st_formats)
    in_dict["te_type_top_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_top." + _producer_fmt(_st_formats)
    in_dict["te_type_group_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_by_group." + _producer_fmt(_st_formats)
    in_dict["tecount_matrix"] = indir + f"/counts/TEtranscripts/{wildcards.genome}/{wildcards.genome}_TEcount.tsv"
    in_dict["tecount_name_matrix"] = indir + f"/counts/TEtranscripts/{wildcards.genome}/{wildcards.genome}_TEcount_name.tsv"
    in_dict["featurecounts_matrix"] = indir + f"/counts/featureCounts/{wildcards.genome}/{wildcards.genome}_featureCounts.tsv"
    in_dict["gsva_heatmap"] = indir + f"/function/{wildcards.genome}/gsva/gsva_heatmap." + _producer_fmt(gsva_fmt)
    in_dict["gsva_scores"] = indir + f"/function/{wildcards.genome}/gsva/gsva_scores.tsv"
    in_dict["fusion_summary"] = indir + f"/fusion/{wildcards.genome}/arriba_report/per_sample_summary.tsv"
    in_dict["recurrent_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/recurrent_fusions.tsv"
    in_dict["high_medium_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/high_medium_confidence_fusions.tsv"
    in_dict["inframe_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/inframe_fusions.tsv"
    in_dict["fusion_figs"] = [indir + f"/fusion/{wildcards.genome}/arriba_report/figures/{base}." + _producer_fmt(_arriba_formats) for base in FUSION_FIG_BASES]
    contrasts = list(config.get("comparisons", {}).get(wildcards.genome, {}).keys())
    if not contrasts:
        logger.warning(f"No contrasts found for genome {wildcards.genome} in comparisons configuration.")
    in_dict["contrast_group"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/group.tsv" for contrast in contrasts]
    in_dict["contrast_pca"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/PCA/{contrast}.cpmPCA." + _producer_fmt(_deseq2_formats) for contrast in contrasts]
    in_dict["contrast_gene_volcano"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/volcano/{contrast}.TEcount_Gene_volcano." + _producer_fmt(_deseq2_formats) for contrast in contrasts]
    in_dict["contrast_gene_heatmap"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/heatmap/{contrast}.TEcount_Gene_updown." + _producer_fmt(_deseq2_formats) for contrast in contrasts]
    in_dict["contrast_gene_updown"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/upDown/{contrast}.TEcount_Gene_updown.tsv" for contrast in contrasts]
    in_dict["contrast_te_updown"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/upDown/{contrast}.TEcount_TE_updown.tsv" for contrast in contrasts]
    in_dict["contrast_gene_te_updown"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/upDown/{contrast}.TEcount_Gene_TE_updown.tsv" for contrast in contrasts]
    in_dict["contrast_gene_name"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/{contrast}.TEcount_Gene.name.tsv" for contrast in contrasts]
    in_dict["contrast_te_name"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/{contrast}.TEcount_TE.name.tsv" for contrast in contrasts]
    in_dict["contrast_gene_te_name"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/TEcount/{contrast}.TEcount_Gene_TE.name.tsv" for contrast in contrasts]
    in_dict["func_go_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_back_to_back." + _producer_fmt(go_kegg_fmt) for contrast in contrasts]
    in_dict["func_kegg_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_back_to_back." + _producer_fmt(go_kegg_fmt) for contrast in contrasts]
    in_dict["func_go_up"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_up.csv" for contrast in contrasts]
    in_dict["func_go_down"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_down.csv" for contrast in contrasts]
    in_dict["func_kegg_up"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_up.csv" for contrast in contrasts]
    in_dict["func_kegg_down"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_down.csv" for contrast in contrasts]
    in_dict["func_up_genes"] = [indir + f"/function/{wildcards.genome}/{contrast}/up_genes.txt" for contrast in contrasts]
    in_dict["func_down_genes"] = [indir + f"/function/{wildcards.genome}/{contrast}/down_genes.txt" for contrast in contrasts]
    in_dict["func_gsea_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/GSEA/TEcount_Gene_GSEA." + _producer_fmt(gsea_fmt) for contrast in contrasts]
    in_dict["func_gsea_csv"] = [indir + f"/function/{wildcards.genome}/{contrast}/GSEA/TEcount_Gene_GSEA.csv" for contrast in contrasts]
    return in_dict

rule generate_report:
    input:
        unpack(get_input_for_generate_report)
    output:
        report = outdir + "/{genome}/RNAseq_report.pptx",
        file_inventory = outdir + "/{genome}/RNAseq_report_files.xlsx",
        figures_manifest = outdir + "/{genome}/figures/manifest.tsv",
    log:
        logdir + "/RNAseq_report/{genome}/RNAseq_report.log"
    threads: rule_res("generate_report").get("threads", 1)
    benchmark: logdir + "/RNAseq_report/{genome}/RNAseq_report.benchmark.tsv",
    resources:
        mem_mb=rule_res("generate_report").get("mem_mb", 8000),
    conda:
        "../RNAseq_report.yaml"
    container:
        sif("../RNAseq_report.yaml", "generate_report")
    params:
        samples = lambda wildcards: genome_paired_samples.get(wildcards.genome, []) + genome_single_samples.get(wildcards.genome, []),
        paired_samples = lambda wildcards: genome_paired_samples.get(wildcards.genome, []),
        single_samples = lambda wildcards: genome_single_samples.get(wildcards.genome, []),
        contrasts = lambda wildcards: list(config.get("comparisons", {}).get(wildcards.genome, {}).keys()),
        title = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("title") or "RNA-seq Analysis Report",
        subtitle = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("date") or "",
        lang = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("lang") or "zh",
        formats = _formats,
        img_dir = lambda wildcards: outdir + f"/{wildcards.genome}/figures",
        script = REPORT_SCRIPT,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("RNAseq_generate_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start RNAseq report generation for genome {wildcards.genome} at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)
            os.makedirs(params.img_dir, exist_ok=True)
            script_path = receipt_script(f"RNAseq_report_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(report_dir, f"RNAseq_report.sh")
            genome_dir = f"{outdir}/{wildcards.genome}"
            for subdir in ["transcripts", "fusion", "diff_expression", "function"]:
                src = f"{indir}/{subdir}/{wildcards.genome}"
                dst = f"{genome_dir}/{subdir}"
                if SCRIPT_ONLY:
                    script_only_touch(output)
                else:
                    shell(f"mkdir -p $(dirname {dst}) && ln -sfn {src} {dst}")
            for count_tool in ["TEtranscripts", "featureCounts"]:
                src = f"{indir}/counts/{count_tool}/{wildcards.genome}"
                dst = f"{genome_dir}/counts/{count_tool}"
                if SCRIPT_ONLY:
                    script_only_touch(output)
                else:
                    shell(f"mkdir -p $(dirname {dst}) && ln -sfn {src} {dst}")
            cmd = [
                "python3", params.script,
                "--analysis-dir", genome_dir,
                "--output", str(output.report),
                "--file-inventory", str(output.file_inventory),
                "--title", params.title,
                "--subtitle", params.subtitle,
                "--pipeline", params.pipeline,
                "--genome", wildcards.genome,
                "--date", params.date,
                "--lang", params.lang,
                "--img-dir", params.img_dir,
            ]
            if params.samples:
                cmd.extend(["--samples", *params.samples])
            if params.paired_samples:
                cmd.extend(["--paired-samples", *params.paired_samples])
            if params.single_samples:
                cmd.extend(["--single-samples", *params.single_samples])
            if params.contrasts:
                cmd.extend(["--contrasts", *params.contrasts])
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join([shlex.quote(str(x)) for x in cmd]) + "\n")
                f.write("echo 'RNAseq report generation completed successfully.'\n")
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"RNAseq report generation failed: {e}")
            raise RuntimeError(f"RNAseq report generation failed: {e}, please check the log file: {log_path} for details.")

rule report_result:
    input:
        report = outdir + "/{genome}/RNAseq_report.pptx",
        file_inventory = outdir + "/{genome}/RNAseq_report_files.xlsx",
        figures_manifest = outdir + "/{genome}/figures/manifest.tsv",
