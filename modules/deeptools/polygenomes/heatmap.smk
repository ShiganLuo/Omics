include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("Params", {}).get("computeMatrix", {}).get("samples", None) or config.get("samples", [])
bigwig_dir = config.get("bigwig_dir", "")
sample_ip_input_map = config.get("sample_ip_input_map", {})

# --- Config ---
regions_cfg = config.get("Params", {}).get("computeMatrix", {}).get("regions", "tss")
MODULE_DIR = os.path.join(config.get("ROOT_DIR", "."), "modules", "deeptools_heatmap")
HEATMAP_SCRIPT = os.path.join(MODULE_DIR, "bin", "run_heatmap.py")

# --- Image formats (forwarded to plotHeatmap --plotFileFormat) ---
# ① parse 期读入 + 规范化（str 收成单元素 list；合法值由 schema 声明）
# 参数精确到 rule:formats 按 rule 独立配置(Params.plotHeatmap.<rule_name>.formats,
# 无 module 级共享回退,config 全量路径内联读取)。
# 每 rule 一个 parse 期变量;heatmap / heatmap_gene 各自消费;
# deeptools_heatmap_result 聚合 heatmap 输出,展开同一 rule 键。
_formats_heatmap = config.get("Params", {}).get("plotHeatmap", {}).get("heatmap", {}).get("formats", ["png"])
if isinstance(_formats_heatmap, str):
    _formats_heatmap = [_formats_heatmap]
_formats_heatmap_gene = config.get("Params", {}).get("plotHeatmap", {}).get("heatmap_gene", {}).get("formats", ["png"])
if isinstance(_formats_heatmap_gene, str):
    _formats_heatmap_gene = [_formats_heatmap_gene]
_formats_deeptools_heatmap_result = config.get("Params", {}).get("plotHeatmap", {}).get("heatmap", {}).get("formats", ["png"])
if isinstance(_formats_deeptools_heatmap_result, str):
    _formats_deeptools_heatmap_result = [_formats_deeptools_heatmap_result]

def _is_gene_regions():
    return isinstance(regions_cfg, dict) and "genes" in regions_cfg

def _is_te_regions():
    """TE regions use te_gtf and should NOT be merged (each locus independent)."""
    return isinstance(regions_cfg, dict) and regions_cfg.get("gtf") == "te"

def _get_gene_gtf(wildcards):
    if isinstance(regions_cfg, dict) and regions_cfg.get("gtf") == "te":
        return config['genome'][wildcards.genome]['te_gtf']
    return config['genome'][wildcards.genome]['gtf']

def _get_match_by():
    return regions_cfg.get("match_by", "gene_name") if isinstance(regions_cfg, dict) else "gene_name"

def _get_gene_names():
    return regions_cfg.get("genes", []) if isinstance(regions_cfg, dict) else []

def get_bigwig(wildcards):
    return os.path.join(bigwig_dir, wildcards.sample_id, wildcards.sample_id + ".bigwig")

def get_input_bigwig(wildcards):
    inp = sample_ip_input_map.get(wildcards.sample_id)
    return os.path.join(bigwig_dir, inp, inp + ".bigwig") if inp else []

# Heatmap output suffix based on regions mode (tss / peaks / genes)
def _heatmap_suffix():
    if regions_cfg == "tss":
        return "tss"
    elif regions_cfg == "peaks":
        return "peaks"
    elif _is_gene_regions():
        return "genes"
    return "heatmap"

HEATMAP_SUFFIX = _heatmap_suffix()

# ============================================================
# bigwigCompare: IP/Input ratio bigwig
# ============================================================

rule bigwig_ratio:
    """Generate IP/Input ratio bigwig using deeptools bigwigCompare."""
    input:
        ip_bigwig = get_bigwig,
        input_bigwig = get_input_bigwig,
    output:
        ratio_bigwig = outdir + "/{genome}/{sample_id}/{sample_id}_IP_over_Input.bigwig",
    log: logdir + "/{genome}/{sample_id}/bigwig_ratio.log"
    threads: rule_res("bigwig_ratio").get("threads", 4)
    benchmark: logdir + "/{genome}/{sample_id}/bigwig_ratio.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "bigwig_ratio", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda: "../deeptools_heatmap.yaml"
    container: sif("../deeptools_heatmap.yaml", "bigwig_ratio")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("bigwig_ratio", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start bigwig_ratio for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.ratio_bigwig))
            os.makedirs(sample_outdir, exist_ok=True)
            cmd = [
                "bigwigCompare",
                "-b1", input.ip_bigwig,
                "-b2", input.input_bigwig,
                "--operation", "ratio",
                "--pseudocount", "1",
                "-o", output.ratio_bigwig,
                "-p", str(threads),
            ]
            script = os.path.join(RECEIPT_DIR, f"bigwig_ratio_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"bigwig_ratio_{wildcards.sample_id}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write(f'echo "bigwig_ratio for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during bigwig_ratio for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule bigwig_ratio_result:
    """Aggregation rule for bigwig_ratio."""
    input:
        ratio_bigwig = outdir + "/{genome}/{sample_id}/{sample_id}_IP_over_Input.bigwig",
def get_ratio_bigwig(wildcards):
    return os.path.join(outdir, wildcards.genome, wildcards.sample_id, wildcards.sample_id + "_IP_over_Input.bigwig")

# ============================================================
# Heatmap (computeMatrix + plotHeatmap unified)
# ============================================================

def _build_heatmap_cmd(wildcards, input, out_base, threads, title, gene_names=None):
    """Build run_heatmap.py command for a given sample.

    ``out_base`` is the figure base path (no extension); run_heatmap.py writes
    one ``<out_base>.<fmt>`` file per ``--format`` (forwarded to deeptools
    plotHeatmap ``--plotFileFormat``).
    """
    mode, ref_point = config.get("Params", {}).get("computeMatrix", {}).get("mode", "reference-point"), config.get("Params", {}).get("computeMatrix", {}).get("referencePoint", "center")
    if regions_cfg == "tss":
        mode, ref_point = "reference-point", "TSS"
    elif regions_cfg == "peaks":
        mode, ref_point = "reference-point", "center"
    elif _is_gene_regions():
        # genes mode: force scale-regions to show TSS and TES labels
        mode, ref_point = "scale-regions", "TSS"

    tss_bed = config['genome'][wildcards.genome].get('tss_bed')
    cmd = [
        "python3", HEATMAP_SCRIPT,
        "--ip-bigwig", input.ratio_bigwig,
        "--output", out_base,
        "--mode", mode, "--reference-point", ref_point,
        "--before", str(config.get("Params", {}).get("computeMatrix", {}).get("before", 3000)),
        "--after", str(config.get("Params", {}).get("computeMatrix", {}).get("after", 3000)),
        "--upstream", str(config.get("Params", {}).get("computeMatrix", {}).get("upstream", 3000)),
        "--downstream", str(config.get("Params", {}).get("computeMatrix", {}).get("downstream", 3000)),
        "--body-length", str(config.get("Params", {}).get("computeMatrix", {}).get("bodyLength", 5000)),
        "--bin-size", str(config.get("Params", {}).get("computeMatrix", {}).get("binSize", 10)),
        "--sort-using", config.get("Params", {}).get("computeMatrix", {}).get("sortUsing", "mean"),
        "--top-n", str(config.get("Params", {}).get("computeMatrix", {}).get("top_n", 0)),
        "--threads", str(threads),
        "--title", title,
        "--color-map", config.get("Params", {}).get("plotHeatmap", {}).get("colorMap", "YlOrRd"),
        "--height", str(config.get("Params", {}).get("plotHeatmap", {}).get("heatmapHeight", 15)),
        "--width", str(config.get("Params", {}).get("plotHeatmap", {}).get("heatmapWidth", 8)),
        "--what-to-show", config.get("Params", {}).get("plotHeatmap", {}).get("whatToShow", "plot, heatmap and colorbar"),
    ]

    # Region source + mode selection based on regions_cfg
    if regions_cfg == "tss" and tss_bed and os.path.isfile(tss_bed):
        cmd += ["--regions", tss_bed]
    elif regions_cfg == "tss":
        cmd += ["--gtf", input.gtf, "--region-mode", "tss",
                "--tss-flank", str(config.get("Params", {}).get("computeMatrix", {}).get("tss_flank", 1000))]
    elif regions_cfg == "peaks":
        # Peaks have no directionality — use center reference point
        cmd += ["--regions", input.regions]
    elif _is_gene_regions():
        cmd += ["--gtf", input.gtf, "--region-mode", "genes",
                "--match-by", _get_match_by()]
        names = gene_names if gene_names is not None else _get_gene_names()
        for name in names:
            cmd += ["--gene-names", name]
        # TE subfamilies: don't merge loci, each locus gets its own row with correct strand
        # Gene list: merge exons into gene body
        if not _is_te_regions():
            cmd += ["--merge"]

    return cmd

def _get_gtf_for_heatmap(wildcards):
    """Return GTF path when mode needs it, empty list otherwise."""
    tss_bed = config['genome'][wildcards.genome].get('tss_bed')
    if regions_cfg == "tss":
        return [] if (tss_bed and os.path.isfile(tss_bed)) else (config['genome'][wildcards.genome].get('gtf') or [])
    if _is_gene_regions():
        return _get_gene_gtf(wildcards) or []
    return []

def _get_regions_for_heatmap(wildcards):
    """Return pre-made BED only for peaks mode."""
    if regions_cfg == "peaks":
        return os.path.join(indir, wildcards.sample_id, wildcards.sample_id + "_peaks.narrowPeak")
    # tss/genes modes generate BED from GTF inside run_heatmap.py
    tss_bed = config['genome'][wildcards.genome].get('tss_bed')
    if regions_cfg == "tss" and tss_bed and os.path.isfile(tss_bed):
        return tss_bed
    return []

rule heatmap:
    """TSS / peaks / gene-list mode — computeMatrix + plotHeatmap in one step."""
    input:
        ratio_bigwig = get_ratio_bigwig,
        gtf = _get_gtf_for_heatmap,
        regions = _get_regions_for_heatmap,
    output:
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        plots = [outdir + "/{genome}/{sample_id}/{sample_id}_" + HEATMAP_SUFFIX + "_heatmap." + fmt for fmt in _formats_heatmap],
    log: logdir + "/{genome}/{sample_id}/heatmap.log"
    threads: rule_res("heatmap").get("threads", 4)
    benchmark: logdir + "/{genome}/{sample_id}/heatmap.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "heatmap", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    params:
        formats = _formats_heatmap,   # ③ 转发
        out_base = outdir + "/{genome}/{sample_id}/{sample_id}_" + HEATMAP_SUFFIX + "_heatmap",
    conda: "../deeptools_heatmap.yaml"
    container: sif("../deeptools_heatmap.yaml", "heatmap")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("heatmap", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start heatmap for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(params.out_base)
            os.makedirs(sample_outdir, exist_ok=True)
            inp = sample_ip_input_map.get(wildcards.sample_id, "")
            title = f"{wildcards.sample_id} vs {inp}" if inp else wildcards.sample_id
            cmd = _build_heatmap_cmd(wildcards, input, params.out_base, threads, title)
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            script = os.path.join(RECEIPT_DIR, f"heatmap_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"heatmap_{wildcards.sample_id}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write(f'echo "heatmap for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during heatmap for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule heatmap_gene:
    """Per gene/TE name — computeMatrix + plotHeatmap in one step."""
    input:
        ratio_bigwig = get_ratio_bigwig,
        gtf = _get_gene_gtf,
    output:
        plots = [outdir + "/{genome}/{sample_id}/{sample_id}_{gene_name}_heatmap." + fmt for fmt in _formats_heatmap_gene],
    log: logdir + "/{genome}/{sample_id}/heatmap_{gene_name}.log"
    threads: rule_res("heatmap_gene").get("threads", 4)
    benchmark: logdir + "/{genome}/{sample_id}/heatmap_{gene_name}.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "heatmap_gene", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    params:
        formats = _formats_heatmap_gene,   # ③ 转发
        out_base = outdir + "/{genome}/{sample_id}/{sample_id}_{gene_name}_heatmap",
    conda: "../deeptools_heatmap.yaml"
    container: sif("../deeptools_heatmap.yaml", "heatmap_gene")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("heatmap_gene", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start heatmap_gene for sample {wildcards.sample_id} {wildcards.gene_name} at {current_time}")
            sample_outdir = os.path.dirname(params.out_base)
            os.makedirs(sample_outdir, exist_ok=True)
            inp = sample_ip_input_map.get(wildcards.sample_id, "")
            title = f"{wildcards.sample_id} {wildcards.gene_name}"
            if inp:
                title += f" vs {inp}"
            cmd = _build_heatmap_cmd(wildcards, input, params.out_base, threads, title, gene_names=[wildcards.gene_name])
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            script = os.path.join(RECEIPT_DIR, f"heatmap_gene_{wildcards.sample_id}_{wildcards.gene_name}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"heatmap_gene_{wildcards.sample_id}_{wildcards.gene_name}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write(f'echo "heatmap_gene for {wildcards.sample_id} {wildcards.gene_name} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during heatmap_gene for sample {wildcards.sample_id} {wildcards.gene_name}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

        # ============================================================
        # Result
        # ============================================================

rule deeptools_heatmap_result:
    input:
        # {genome}/{sample_id} stay rule wildcards; only the extension is
        # expanded over formats — comprehension, not expand()
        plots = [outdir + "/{genome}/{sample_id}/{sample_id}_" + HEATMAP_SUFFIX + "_heatmap." + fmt for fmt in _formats_deeptools_heatmap_result],
