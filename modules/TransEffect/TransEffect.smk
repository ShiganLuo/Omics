include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
sample_group_dict = config.get("sample_group_dict", {})
drclass_dir = config.get("drclass_dir", "")
counts_matrix = config.get("counts_matrix", "")
cnv_dir = config.get("cnv_dir", "")
gmt = config.get("gmt")
gtf = config.get("gtf")
Procedure = config.get("Procedure", {})
teffect_params = config.get("Params", {}).get("TransEffect", {}).get("teffect_analyze", {})

# 输出决定流程走向(modules.md): program_overlap 及其图只在配置了 gmt 时声明,
# node.py 同条件注册;未配置时 P5 分析整体跳过,不产生悬空输出。
_has_gmt = bool(gmt)

TEFFECT_TABLES = [
    "rna_support_sites.tsv",
    "vaf_ratios.tsv",
    "variant_gene_delta.tsv",
    "cnv_dosage_summary.tsv",
    "cnv_segments.tsv",
    "teffect_report.md",
]
if _has_gmt:
    TEFFECT_TABLES.append("program_overlap.tsv")

TEFFECT_FIGS = [
    "fig_rna_support_vaf.png",
    "fig_variant_gene_delta.png",
    "fig_cnv_overview.png",
]
if _has_gmt:
    TEFFECT_FIGS.append("fig_program_overlap.png")


def teffect_analyze_inputs(wildcards):
    """声明全部消费文件;可选文件(gmt/gtf)配置了才进 input。"""
    files = {
        "counts": counts_matrix,
    }
    for s in samples:
        files[f"drclass_dna_{s}"] = f"{drclass_dir}/{s}/{s}.dna_anchored.tsv"
        files[f"anchors_{s}"] = f"{drclass_dir}/{s}/anchors.tsv"
        files[f"cnv_{s}"] = f"{cnv_dir}/{s}/{s}.call.cns"
    if gmt:
        files["gmt"] = gmt
    if gtf:
        files["gtf"] = gtf
    return files


rule teffect_analyze:
    """转录效应分析表: RNA 支持位点 VAF / RNA-DNA 比值 / 变异基因 Δ' / CNV 剂量 / P5 交集."""
    input:
        unpack(teffect_analyze_inputs),
    output:
        rna_sites = outdir + "/rna_support_sites.tsv",
        vaf = outdir + "/vaf_ratios.tsv",
        delta = outdir + "/variant_gene_delta.tsv",
        cnv_summary = outdir + "/cnv_dosage_summary.tsv",
        cnv_segments = outdir + "/cnv_segments.tsv",
        report = outdir + "/teffect_report.md",
        program = (outdir + "/program_overlap.tsv") if config.get("gmt") else [],
    log:
        logdir + "/teffect_analyze.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "TransEffect", "bin", "transeffect.py"),
        py_proc=Procedure.get("python3") or "python3",
        sample_group_args=[f"{s}={sample_group_dict.get(s, 'G')}" for s in samples],
        dna_af_min=teffect_params.get("dna_af_min", 0.1),
        dna_dp_min=teffect_params.get("dna_dp_min", 10),
        min_cpm=teffect_params.get("min_cpm", 1.0),
        min_expressed=teffect_params.get("min_expressed", 2),
        null_reps=teffect_params.get("null_reps", 20),
        null_plot_reps=teffect_params.get("null_plot_reps", 3),
        seed=teffect_params.get("seed", 42),
        sex_chroms=teffect_params.get("sex_chroms"),
        gene_sets=teffect_params.get("gene_sets"),
    benchmark: logdir + "/teffect_analyze.benchmark.tsv",
    resources:
        mem_mb=rule_res("teffect_analyze").get("mem_mb", 8000),
    conda:
        "TransEffect.yaml"
    container:
        sif("TransEffect.yaml", "teffect_analyze")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("teffect_analyze", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start teffect analyze at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"teffect_analyze_{current_time}.sh")
            plain_script = os.path.join(outdir, "teffect_analyze.sh")

            cmd = [
                params.py_proc, params.py, "analyze",
                "--samples", ",".join(samples),
                "--drclass-dir", str(drclass_dir),
                "--counts", str(input.counts),
                "--cnv-dir", str(cnv_dir),
                "--outdir", str(outdir),
                "--dna-af-min", str(params.dna_af_min),
                "--dna-dp-min", str(params.dna_dp_min),
                "--min-cpm", str(params.min_cpm),
                "--min-expressed", str(params.min_expressed),
                "--null-reps", str(params.null_reps),
                "--null-plot-reps", str(params.null_plot_reps),
                "--seed", str(params.seed),
            ]
            for sg in params.sample_group_args:
                cmd += ["--sample-group", sg]
            if params.sex_chroms:
                cmd += ["--sex-chroms", str(params.sex_chroms)]
            if params.gene_sets:
                cmd += ["--gene-sets", str(params.gene_sets)]
            if getattr(input, "gmt", ""):
                cmd += ["--gmt", str(input.gmt)]
            if getattr(input, "gtf", ""):
                cmd += ["--gtf", str(input.gtf)]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write('echo "teffect analyze was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during teffect analyze: {e}\n")
            raise RuntimeError(f"Error occurred during teffect analyze: {e}, you can check the log file {log_path} for more details.\n")


rule teffect_figures:
    """读 analyze 产物画图(重绘不重算)."""
    input:
        tables=expand(outdir + "/{name}", name=TEFFECT_TABLES),
    output:
        f_rna_sites = outdir + "/fig_rna_support_vaf.png",
        f_delta = outdir + "/fig_variant_gene_delta.png",
        f_cnv = outdir + "/fig_cnv_overview.png",
        f_program = (outdir + "/fig_program_overlap.png") if config.get("gmt") else [],
    log:
        logdir + "/teffect_figures.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "TransEffect", "bin", "transeffect.py"),
        py_proc=Procedure.get("python3") or "python3",
    benchmark: logdir + "/teffect_figures.benchmark.tsv",
    resources:
        mem_mb=rule_res("teffect_figures").get("mem_mb", 4000),
    conda:
        "TransEffect.yaml"
    container:
        sif("TransEffect.yaml", "teffect_figures")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("teffect_figures", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start teffect figures at {current_time}")
            script_path = os.path.join(RECEIPT_DIR, f"teffect_figures_{current_time}.sh")
            plain_script = os.path.join(outdir, "teffect_figures.sh")

            cmd = [
                params.py_proc, params.py, "figures",
                "--outdir", str(outdir),
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write('echo "teffect figures was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during teffect figures: {e}\n")
            raise RuntimeError(f"Error occurred during teffect figures: {e}, you can check the log file {log_path} for more details.\n")
