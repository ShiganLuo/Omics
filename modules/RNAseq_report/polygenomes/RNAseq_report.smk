include: "../../common/common.smk"

indir = config.get("indir", "data")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
genome_paired_samples = config.get("genome_paired_samples", {})
genome_single_samples = config.get("genome_single_samples", {})

# 图片格式（report 模块 figures/<name>.<fmt> 多格式输出）① parse 期读入 +
# 规范化（str 收成单元素 list；合法值由 schema 声明）
_formats = config.get("Params", {}).get("report", {}).get("formats") or ["png"]
if isinstance(_formats, str):
    _formats = [_formats]

# ── 上游产出图 raster-first 变体解析(modules.md #9 formats 规范)──
# report input 消费上游生产者的图;生产者按 Params.<producer>.formats 每格式
# 写一个 `<base>.<fmt>`。raster-first(png → jpg → jpeg → tif → tiff,都没有
# 则取配置列表首个)选定一个具体扩展名,保证 DAG 只要求生产者真正会写出的
# 那一个变体。默认配置下解析结果与历史硬编码扩展名逐字节一致。function 的
# formats 在 Params.function.<genome>.<section> 下,按 wildcards.genome 在
# get_input_for_generate_report 里逐 genome 解析(_function_fmt)。
_RASTER_FMT_ORDER = ("png", "jpg", "jpeg", "tif", "tiff")


def _producer_fmt(formats):
    """Raster-first extension among one producer's configured formats."""
    fmts = formats or ["png"]
    if isinstance(fmts, dict):
        # schema-metadata 占位形状 {"type": ..., "default": [...]} — 取 default
        fmts = fmts.get("default") or ["png"]
    if isinstance(fmts, str):
        fmts = [fmts]
    for fmt in _RASTER_FMT_ORDER:
        if fmt in fmts:
            return fmt
    return fmts[0]


_st_formats = (config.get("Params", {}).get("StringTie") or {}).get("formats") or ["png"]
_deseq2_formats = (config.get("Params", {}).get("DESeq2") or {}).get("formats") or ["png"]
_arriba_formats = (config.get("Params", {}).get("arriba") or {}).get("formats") or ["png"]


def _function_fmt(genome, section, default):
    """Per-genome Params.function.<genome>.<section>.formats raster-first pick."""
    gp = (config.get("Params", {}).get("function") or {}).get(genome) or {}
    return _producer_fmt((gp.get(section) or {}).get("formats") or default)


REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "RNAseq_report", "bin", "generate_report.py")
# arriba fusion figure base names (no extension) — 变体扩展名由
# Params.arriba.formats raster-first 解析(见 _producer_fmt)
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
    in_dict["te_sample_summary"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_sample_summary.tsv"
    in_dict["te_group_summary"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_group_summary.tsv"
    in_dict["te_type_counts"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_counts.tsv"
    in_dict["te_group_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_group_stacked." + _producer_fmt(_st_formats)
    in_dict["te_type_top_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_top." + _producer_fmt(_st_formats)
    in_dict["te_type_group_plot"] = indir + f"/transcripts/{wildcards.genome}/TE_chimeric/TE_chimeric_te_type_by_group." + _producer_fmt(_st_formats)
    in_dict["tecount_matrix"] = indir + f"/counts/TEtranscripts/{wildcards.genome}/{wildcards.genome}_TEcount.tsv"
    in_dict["tecount_name_matrix"] = indir + f"/counts/TEtranscripts/{wildcards.genome}/{wildcards.genome}_TEcount_name.tsv"
    in_dict["featurecounts_matrix"] = indir + f"/counts/featureCounts/{wildcards.genome}/{wildcards.genome}_featureCounts.tsv"
    in_dict["gsva_heatmap"] = indir + f"/function/{wildcards.genome}/gsva/gsva_heatmap." + _function_fmt(wildcards.genome, "gsva", ["png"])
    in_dict["gsva_scores"] = indir + f"/function/{wildcards.genome}/gsva/gsva_scores.tsv"
    in_dict["fusion_summary"] = indir + f"/fusion/{wildcards.genome}/arriba_report/per_sample_summary.tsv"
    in_dict["recurrent_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/recurrent_fusions.tsv"
    in_dict["high_medium_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/high_medium_confidence_fusions.tsv"
    in_dict["inframe_fusions"] = indir + f"/fusion/{wildcards.genome}/arriba_report/inframe_fusions.tsv"
    in_dict["fusion_figs"] = [indir + f"/fusion/{wildcards.genome}/arriba_report/figures/{base}." + _producer_fmt(_arriba_formats) for base in FUSION_FIG_BASES]
    contrasts = list(config.get("Params", {}).get("DESeq2", {}).get("group_pairs", {}).get(wildcards.genome, {}).keys())
    if not contrasts:
        logger.warning(f"No contrasts found for genome {wildcards.genome} in group_pairs configuration.")
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
    in_dict["func_go_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_back_to_back." + _function_fmt(wildcards.genome, "go_kegg", ["png"]) for contrast in contrasts]
    in_dict["func_kegg_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_back_to_back." + _function_fmt(wildcards.genome, "go_kegg", ["png"]) for contrast in contrasts]
    in_dict["func_go_up"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_up.csv" for contrast in contrasts]
    in_dict["func_go_down"] = [indir + f"/function/{wildcards.genome}/{contrast}/go_down.csv" for contrast in contrasts]
    in_dict["func_kegg_up"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_up.csv" for contrast in contrasts]
    in_dict["func_kegg_down"] = [indir + f"/function/{wildcards.genome}/{contrast}/kegg_down.csv" for contrast in contrasts]
    in_dict["func_up_genes"] = [indir + f"/function/{wildcards.genome}/{contrast}/up_genes.txt" for contrast in contrasts]
    in_dict["func_down_genes"] = [indir + f"/function/{wildcards.genome}/{contrast}/down_genes.txt" for contrast in contrasts]
    in_dict["func_gsea_plot"] = [indir + f"/function/{wildcards.genome}/{contrast}/GSEA/TEcount_Gene_GSEA." + _function_fmt(wildcards.genome, "gsea", ["jpeg"]) for contrast in contrasts]
    in_dict["func_gsea_csv"] = [indir + f"/function/{wildcards.genome}/{contrast}/GSEA/TEcount_Gene_GSEA.csv" for contrast in contrasts]
    return in_dict

rule generate_report:
    input:
        unpack(get_input_for_generate_report)
    output:
        report = outdir + "/{genome}/RNAseq_report.pptx",
        file_inventory = outdir + "/{genome}/RNAseq_report_files.xlsx",
        # 按格式展开的 figures/<name>.<fmt> 由脚本落盘，DAG 用文件清单
        # manifest.tsv 锚定（no directory() outputs）
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
        contrasts = lambda wildcards: list(config.get("Params", {}).get("DESeq2", {}).get("group_pairs", {}).get(wildcards.genome, {}).keys()),
        title = config.get("Params", {}).get("report", {}).get("title") or "RNA-seq Analysis Report",
        subtitle = config.get("Params", {}).get("report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("date") or "",
        lang = config.get("Params", {}).get("report", {}).get("lang") or "zh",
        formats = _formats,   # ③ 转发
        img_dir = lambda wildcards: outdir + f"/{wildcards.genome}/figures",
        script = REPORT_SCRIPT,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("RNAseq_generate_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start RNAseq report generation for genome {wildcards.genome} at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)
            os.makedirs(params.img_dir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"RNAseq_report_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(report_dir, f"RNAseq_report.sh")
            # Create symlink structure for generate_report.py compatibility
            # Script expects: analysis_dir/transcripts/..., analysis_dir/fusion/...
            # Actual paths: {outdir}/transcripts/{genome}/..., {outdir}/fusion/{genome}/...
            genome_dir = f"{outdir}/{wildcards.genome}"
            for subdir in ["transcripts", "fusion", "diff_expression", "function"]:
                src = f"{indir}/{subdir}/{wildcards.genome}"
                dst = f"{genome_dir}/{subdir}"
                shell(f"mkdir -p $(dirname {dst}) && ln -sfn {src} {dst}")
            # counts: TEtranscripts and featureCounts have separate subdirs
            for count_tool in ["TEtranscripts", "featureCounts"]:
                src = f"{indir}/counts/{count_tool}/{wildcards.genome}"
                dst = f"{genome_dir}/counts/{count_tool}"
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
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"RNAseq report generation failed: {e}")
            raise RuntimeError(f"RNAseq report generation failed: {e}, please check the log file: {log_path} for details.")

rule report_result:
    input:
        report = outdir + "/{genome}/RNAseq_report.pptx",
        file_inventory = outdir + "/{genome}/RNAseq_report_files.xlsx",
        figures_manifest = outdir + "/{genome}/figures/manifest.tsv",

    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "report_result", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
