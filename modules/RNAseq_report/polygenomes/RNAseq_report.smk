include: "../../common/common.smk"
indir = config.get("indir", "data")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
genome_paired_samples = config.get("genome_paired_samples", {})
genome_single_samples = config.get("genome_single_samples", {})

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "RNAseq_report", "bin", "generate_report.py")
FUSION_FIGS = [
    "fig1_per_sample_counts.png",
    "fig2_fusion_type_distribution.png",
    "fig3_type_by_sample.png",
    "fig4_reading_frame.png",
    "fig5_recurrent_heatmap.png",
    "fig6_inframe_support.png",
]

def get_input_for_generate_report(wildcards):
    logger.info(f"[get_input_for_generate_report] called with wildcards: {wildcards}")
    in_dict = {}
    in_dict["te_sample_summary"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_sample_summary.tsv"
    in_dict["te_group_summary"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_group_summary.tsv"
    in_dict["te_type_counts"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_counts.tsv"
    in_dict["te_group_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_group_stacked.png"
    in_dict["te_type_top_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_top.png"
    in_dict["te_type_group_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_by_group.png"
    in_dict["tecount_matrix"] = indir + f"/counts/{wildcards.genome}/TEcount/all_TEcount.tsv"
    in_dict["fusion_summary"] = indir + f"/fusion/{wildcards.genome}/arriba_report/per_sample_summary.tsv"
    in_dict["recurrent_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/recurrent_fusions.tsv"
    in_dict["high_medium_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/high_medium_confidence_fusions.tsv"
    in_dict["inframe_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/inframe_fusions.tsv"
    in_dict["fusion_figs"] = [indir + f"/fusion/{wildcards.genome}/arriba_report/figures/{fig}" for fig in FUSION_FIGS]
    contrasts = list(config.get("Params", {}).get("DESeq2", {}).get("group_pairs", {}).get(wildcards.genome, {}).keys())
    if not contrasts:
        logger.warning(f"No contrasts found for genome {wildcards.genome} in group_pairs configuration.")
    in_dict["contrast_group"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/group.tsv" for contrast in contrasts]
    in_dict["contrast_pca"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/PCA/{contrast}.cpmPCA.png" for contrast in contrasts]
    in_dict["contrast_gene_volcano"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/volcano/{contrast}.TEcount_Gene_volcano.png" for contrast in contrasts]
    in_dict["contrast_gene_heatmap"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/heatmap/{contrast}.TEcount_Gene_updown.png" for contrast in contrasts]
    in_dict["contrast_gene_updown"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/upDown/{contrast}.TEcount_Gene_updown.tsv" for contrast in contrasts]
    in_dict["contrast_te_updown"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/upDown/{contrast}.TEcount_TE_updown.tsv" for contrast in contrasts]
    in_dict["contrast_gene_te_updown"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/upDown/{contrast}.TEcount_Gene_TE_updown.tsv" for contrast in contrasts]
    in_dict["contrast_gene_name"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/{contrast}.TEcount_Gene.name.tsv" for contrast in contrasts]
    in_dict["contrast_te_name"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/{contrast}.TEcount_TE.name.tsv" for contrast in contrasts]
    in_dict["contrast_gene_te_name"] = [indir + f"/diff_expression/{wildcards.genome}/{contrast}/{contrast}.TEcount_Gene_TE.name.tsv" for contrast in contrasts]
    in_dict["func_go_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_back_to_back.png" for contrast in contrasts]
    in_dict["func_kegg_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_back_to_back.png" for contrast in contrasts]
    in_dict["func_go_up"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_up.csv" for contrast in contrasts]
    in_dict["func_go_down"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_down.csv" for contrast in contrasts]
    in_dict["func_kegg_up"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_up.csv" for contrast in contrasts]
    in_dict["func_kegg_down"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_down.csv" for contrast in contrasts]
    in_dict["func_up_genes"] = [indir + f"/function/{wildcards.genome}/{contrast}/up_genes.txt" for contrast in contrasts]
    in_dict["func_down_genes"] = [indir + f"/function/{wildcards.genome}/{contrast}/down_genes.txt" for contrast in contrasts]
    in_dict["func_gsea_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/GSEA/TEcount_Gene_GSEA.jpeg" for contrast in contrasts]
    in_dict["func_gsea_csv"] = [indir + f"/function/{wildcards.genome}/{contrast}/GSEA/TEcount_Gene_GSEA.csv" for contrast in contrasts]
    return in_dict

rule generate_report:
    input:
        unpack(get_input_for_generate_report)
    output:
        report = outdir + "/{genome}/RNAseq_report.pptx",
        file_inventory = outdir + "/{genome}/RNAseq_report_files.xlsx",
    log:
        logdir + "/{genome}/RNAseq_report.log"
    threads: 1
    conda:
        "../RNAseq_report.yaml"
    container:
        sif("../RNAseq_report.yaml")
    params:
        samples = lambda wildcards: genome_paired_samples.get(wildcards.genome, []) + genome_single_samples.get(wildcards.genome, []),
        paired_samples = lambda wildcards: genome_paired_samples.get(wildcards.genome, []),
        single_samples = lambda wildcards: genome_single_samples.get(wildcards.genome, []),
        contrasts = lambda wildcards: list(config.get("Params", {}).get("DESeq2", {}).get("group_pairs", {}).get(wildcards.genome, {}).keys()),
        title = config.get("Params", {}).get("report", {}).get("title") or "RNA-seq Analysis Report",
        subtitle = config.get("Params", {}).get("report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("date") or "",
        lang = config.get("Params", {}).get("report", {}).get("lang") or "zh",
        img_dir = lambda wildcards: outdir + f"/{wildcards.genome}/ppt_results",
        script = REPORT_SCRIPT,
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("RNAseq_generate_report", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start RNAseq report generation for genome {wildcards.genome} at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)
            os.makedirs(params.img_dir, exist_ok=True)
            script_path = os.path.join(report_dir, f"RNAseq_report_{current_time}.sh")
            # Create symlink structure for generate_report.py compatibility
            # Script expects: analysis_dir/transcripts/..., analysis_dir/fusion/...
            # Actual paths: {outdir}/transcripts/{genome}/..., {outdir}/fusion/{genome}/...
            genome_dir = f"{outdir}/{wildcards.genome}"
            for subdir in ["transcripts", "fusion", "diff_expression", "counts", "function"]:
                src = f"{indir}/{subdir}/{wildcards.genome}"
                dst = f"{genome_dir}/{subdir}"
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
            with open(script_path, "w") as fh:
                fh.write("#!/bin/bash\n")
                fh.write("set -euo pipefail\n")
                fh.write(" ".join([shlex.quote(str(x)) for x in cmd]) + "\n")
            shell(f"bash {script_path} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as fh:
                fh.write(f"RNAseq report generation failed: {e}\n")
            logger.error(f"RNAseq report generation failed: {e}")
            raise


rule report_result:
    input:
        report = outdir + "/{genome}/RNAseq_report.pptx",
        file_inventory = outdir + "/{genome}/RNAseq_report_files.xlsx",
