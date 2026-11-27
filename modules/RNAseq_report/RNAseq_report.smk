include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
contrasts = config.get("contrasts", [])

genome_label = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("genome") or "all"

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
_go_kegg_formats = config.get("Params", {}).get("function", {}).get("generate_report", {}).get("go_kegg_formats") or ["png"]
if isinstance(_go_kegg_formats, str):
    _go_kegg_formats = [_go_kegg_formats]
_gsea_formats = config.get("Params", {}).get("function", {}).get("generate_report", {}).get("gsea_formats") or ["jpeg"]
if isinstance(_gsea_formats, str):
    _gsea_formats = [_gsea_formats]
_gsva_formats = config.get("Params", {}).get("function", {}).get("generate_report", {}).get("gsva_formats") or ["png"]
if isinstance(_gsva_formats, str):
    _gsva_formats = [_gsva_formats]

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "RNAseq_report", "bin", "generate_report.py")
FUSION_FIG_BASES = [
    "fig1_per_sample_counts",
    "fig2_fusion_type_distribution",
    "fig3_type_by_sample",
    "fig4_reading_frame",
    "fig5_recurrent_heatmap",
    "fig6_inframe_support",
]

rule generate_report:
    input:
        te_sample_summary = outdir + "/transcripts/TE_chimeric/TE_chimeric_sample_summary.tsv",
        te_group_summary = outdir + "/transcripts/TE_chimeric/TE_chimeric_group_summary.tsv",
        te_type_counts = outdir + "/transcripts/TE_chimeric/TE_chimeric_te_type_counts.tsv",
        te_group_plot = outdir + "/transcripts/TE_chimeric/TE_chimeric_group_stacked." + _producer_fmt(_st_formats),
        te_type_top_plot = outdir + "/transcripts/TE_chimeric/TE_chimeric_te_type_top." + _producer_fmt(_st_formats),
        te_type_group_plot = outdir + "/transcripts/TE_chimeric/TE_chimeric_te_type_by_group." + _producer_fmt(_st_formats),
        tecount_matrix = outdir + f"/counts/TEtranscripts/{genome_label}_TEcount.tsv",
        tecount_name_matrix = outdir + f"/counts/TEtranscripts/{genome_label}_TEcount_name.tsv",
        featurecounts_matrix = outdir + f"/counts/featureCounts/{genome_label}_featureCounts.tsv",
        gsva_heatmap = outdir + "/function/gsva/gsva_heatmap." + _producer_fmt(_gsva_formats),
        gsva_scores = outdir + "/function/gsva/gsva_scores.tsv",
        fusion_summary = outdir + "/fusion/arriba_report/per_sample_summary.tsv",
        recurrent_fusions = outdir + "/fusion/arriba_report/recurrent_fusions.tsv",
        high_medium_fusions = outdir + "/fusion/arriba_report/high_medium_confidence_fusions.tsv",
        inframe_fusions = outdir + "/fusion/arriba_report/inframe_fusions.tsv",
        fusion_figs = expand(outdir + "/fusion/arriba_report/figures/{base}." + _producer_fmt(_arriba_formats),
                             base=FUSION_FIG_BASES),
        contrast_group = expand(outdir + "/diff_expression/{contrast}/group.tsv", contrast=contrasts),
        contrast_pca = expand(outdir + "/diff_expression/{contrast}/PCA/{contrast}.cpmPCA." + _producer_fmt(_deseq2_formats), contrast=contrasts),
        contrast_gene_volcano = expand(outdir + "/diff_expression/{contrast}/volcano/{contrast}.TEcount_Gene_volcano." + _producer_fmt(_deseq2_formats), contrast=contrasts),
        contrast_gene_heatmap = expand(outdir + "/diff_expression/{contrast}/heatmap/{contrast}.TEcount_Gene_updown." + _producer_fmt(_deseq2_formats), contrast=contrasts),
        contrast_gene_updown = expand(outdir + "/diff_expression/{contrast}/upDown/{contrast}.TEcount_Gene_updown.tsv", contrast=contrasts),
        contrast_te_updown = expand(outdir + "/diff_expression/{contrast}/upDown/{contrast}.TEcount_TE_updown.tsv", contrast=contrasts),
        contrast_gene_te_updown = expand(outdir + "/diff_expression/{contrast}/upDown/{contrast}.TEcount_Gene_TE_updown.tsv", contrast=contrasts),
        contrast_gene_name = expand(outdir + "/diff_expression/{contrast}/TEcount/{contrast}.TEcount_Gene.name.tsv", contrast=contrasts),
        contrast_te_name = expand(outdir + "/diff_expression/{contrast}/TEcount/{contrast}.TEcount_TE.name.tsv", contrast=contrasts),
        contrast_gene_te_name =  expand(outdir + "/diff_expression/{contrast}/TEcount/{contrast}.TEcount_Gene_TE.name.tsv", contrast=contrasts),
        func_go_plot = expand(outdir + "/function/{contrast}/go_back_to_back." + _producer_fmt(_go_kegg_formats), contrast=contrasts),
        func_kegg_plot = expand(outdir + "/function/{contrast}/kegg_back_to_back." + _producer_fmt(_go_kegg_formats), contrast=contrasts),
        func_go_up = expand(outdir + "/function/{contrast}/go_up.csv", contrast=contrasts),
        func_go_down = expand(outdir + "/function/{contrast}/go_down.csv", contrast=contrasts),
        func_kegg_up = expand(outdir + "/function/{contrast}/kegg_up.csv", contrast=contrasts),
        func_kegg_down = expand(outdir + "/function/{contrast}/kegg_down.csv", contrast=contrasts),
        func_up_genes = expand(outdir + "/function/{contrast}/up_genes.txt", contrast=contrasts),
        func_down_genes = expand(outdir + "/function/{contrast}/down_genes.txt", contrast=contrasts),
        func_gsea_plot = expand(outdir + "/function/{contrast}/GSEA/TEcount_Gene_GSEA." + _producer_fmt(_gsea_formats), contrast=contrasts),
        func_gsea_csv = expand(outdir + "/function/{contrast}/GSEA/TEcount_Gene_GSEA.csv", contrast=contrasts),
    output:
        report = outdir + "/RNAseq_report.pptx",
        file_inventory = outdir + "/RNAseq_report_files.xlsx",
        figures_manifest = outdir + "/figures/manifest.tsv",
    log:
        logdir + "/RNAseq_report.log"
    threads: rule_res("generate_report").get("threads", 1)
    benchmark: logdir + "/RNAseq_report.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "generate_report", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "RNAseq_report.yaml"
    container:
        sif("RNAseq_report.yaml", "generate_report")
    params:
        samples = samples,
        paired_samples = paired_samples,
        single_samples = single_samples,
        contrasts = contrasts,
        title = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("title") or "RNA-seq Analysis Report",
        subtitle = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("date") or "",
        lang = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("lang") or "zh",
        formats = _formats,
        img_dir = outdir + "/figures",
        script = REPORT_SCRIPT,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("RNAseq_generate_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start RNAseq report generation at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)
            os.makedirs(params.img_dir, exist_ok=True)
            script_path = receipt_script(f"RNAseq_report_{current_time}.sh")
            plain_script = os.path.join(report_dir, f"RNAseq_report.sh")
            cmd = [
                "python3", params.script,
                "--analysis-dir", outdir,
                "--output", str(output.report),
                "--file-inventory", str(output.file_inventory),
                "--title", params.title,
                "--subtitle", params.subtitle,
                "--pipeline", params.pipeline,
                "--genome", params.genome,
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
        report = outdir + "/RNAseq_report.pptx",
        file_inventory = outdir + "/RNAseq_report_files.xlsx",
        figures_manifest = outdir + "/figures/manifest.tsv",
