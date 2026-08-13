include: "../common/common.smk"
import shlex
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
contrasts = config.get("contrasts", [])

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "RNAseq_report", "bin", "generate_report.py")
FUSION_FIGS = [
    "fig1_per_sample_counts.png",
    "fig2_fusion_type_distribution.png",
    "fig3_type_by_sample.png",
    "fig4_reading_frame.png",
    "fig5_recurrent_heatmap.png",
    "fig6_inframe_support.png",
]

rule generate_report:
    input:
        te_sample_summary = outdir + "/transcripts/TE_chimeric/TE_chimeric_sample_summary.tsv",
        te_group_summary = outdir + "/transcripts/TE_chimeric/TE_chimeric_group_summary.tsv",
        te_type_counts = outdir + "/transcripts/TE_chimeric/TE_chimeric_te_type_counts.tsv",
        te_group_plot = outdir + "/transcripts/TE_chimeric/TE_chimeric_group_stacked.png",
        te_type_top_plot = outdir + "/transcripts/TE_chimeric/TE_chimeric_te_type_top.png",
        te_type_group_plot = outdir + "/transcripts/TE_chimeric/TE_chimeric_te_type_by_group.png",
        tecount_matrix = outdir + "/results/counts/TEcount/all_TEcount.tsv",
        fusion_summary = outdir + "/fusion/arriba_report/per_sample_summary.tsv",
        recurrent_fusions = outdir + "/fusion/arriba_report/recurrent_fusions.tsv",
        high_medium_fusions = outdir + "/fusion/arriba_report/high_medium_confidence_fusions.tsv",
        inframe_fusions = outdir + "/fusion/arriba_report/inframe_fusions.tsv",
        fusion_figs = expand(outdir + "/fusion/arriba_report/figures/{fig}", fig=FUSION_FIGS),
        contrast_group = expand(outdir + "/diff_expression/{contrast}/group.tsv", contrast=contrasts),
        contrast_pca = expand(outdir + "/diff_expression/{contrast}/PCA/{contrast}.cpmPCA.png", contrast=contrasts),
        contrast_gene_volcano = expand(outdir + "/diff_expression/{contrast}/volcano/{contrast}.TEcount_Gene_volcano.png", contrast=contrasts),
        contrast_gene_heatmap = expand(outdir + "/diff_expression/{contrast}/heatmap/{contrast}.TEcount_Gene_updown.png", contrast=contrasts),
        contrast_gene_updown = expand(outdir + "/diff_expression/{contrast}/upDown/{contrast}.TEcount_Gene_updown.tsv", contrast=contrasts),
        contrast_te_updown = expand(outdir + "/diff_expression/{contrast}/upDown/{contrast}.TEcount_TE_updown.tsv", contrast=contrasts),
        contrast_gene_te_updown = expand(outdir + "/diff_expression/{contrast}/upDown/{contrast}.TEcount_Gene_TE_updown.tsv", contrast=contrasts),
        contrast_gene_name = expand(outdir + "/diff_expression/{contrast}/{contrast}.TEcount_Gene.name.tsv", contrast=contrasts),
        contrast_te_name = expand(outdir + "/diff_expression/{contrast}/{contrast}.TEcount_TE.name.tsv", contrast=contrasts),
        contrast_gene_te_name =  expand(outdir + "/diff_expression/{contrast}/{contrast}.TEcount_Gene_TE.name.tsv", contrast=contrasts),
        func_go_plot = expand(outdir + "/function/{contrast}/go_back_to_back.png", contrast=contrasts),
        func_kegg_plot = expand(outdir + "/function/{contrast}/kegg_back_to_back.png", contrast=contrasts),
        func_go_up = expand(outdir + "/function/{contrast}/go_up.csv", contrast=contrasts),
        func_go_down = expand(outdir + "/function/{contrast}/go_down.csv", contrast=contrasts),
        func_kegg_up = expand(outdir + "/function/{contrast}/kegg_up.csv", contrast=contrasts),
        func_kegg_down = expand(outdir + "/function/{contrast}/kegg_down.csv", contrast=contrasts),
        func_up_genes = expand(outdir + "/function/{contrast}/up_genes.txt", contrast=contrasts),
        func_down_genes = expand(outdir + "/function/{contrast}/down_genes.txt", contrast=contrasts),
        func_gsea_plot = expand(outdir + "/function/{contrast}/GSEA/TEcount_Gene_GSEA.jpeg", contrast=contrasts),
        func_gsea_csv = expand(outdir + "/function/{contrast}/GSEA/TEcount_Gene_GSEA.csv", contrast=contrasts),
    output:
        report = outdir + "/RNAseq_report.pptx",
        file_inventory = outdir + "/RNAseq_report_files.xlsx",
    log:
        logdir + "/RNAseq_report.log"
    threads: 1
    conda:
        "RNAseq_report.yaml"
    container:
        sif("RNAseq_report.yaml")
    params:
        samples = samples,
        paired_samples = paired_samples,
        single_samples = single_samples,
        contrasts = contrasts,
        title = config.get("Params", {}).get("report", {}).get("title") or "RNA-seq Analysis Report",
        subtitle = config.get("Params", {}).get("report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("date") or "",
        lang = config.get("Params", {}).get("report", {}).get("lang") or "zh",
        img_dir = outdir + "/ppt_results",
        script = REPORT_SCRIPT,
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("RNAseq_generate_report", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start RNAseq report generation at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)
            os.makedirs(params.img_dir, exist_ok=True)
            script_path = os.path.join(report_dir, f"RNAseq_report_{current_time}.sh")
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
        report = outdir + "/RNAseq_report.pptx",
        file_inventory = outdir + "/RNAseq_report_files.xlsx",
