include: "../../common/common.smk"

indir = config.get("indir", "data")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
group_pairs = config.get("group_pairs", {})

# ① parse 期读入 + 规范化（str 收成单元素 list；合法值由 schema 声明）
# R figure formats forwarded to DESeq2.r --format (repeatable). Valid set is
# the R ggsave device set (png pdf svg ps eps tif tiff jpg jpeg); DESeq2.r
# fails fast on anything else (pgf/raw/rgba are matplotlib-only).
# 参数精确到 rule:formats 按 rule 独立配置(config 全量路径读取,
# Params.DESeq2.<rule_name>.formats)。
_formats_TEcount = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TEcount", {}).get("formats") or ["png"]
if isinstance(_formats_TEcount, str):
    _formats_TEcount = [_formats_TEcount]
_formats_TElocal = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal", {}).get("formats") or ["png"]
if isinstance(_formats_TElocal, str):
    _formats_TElocal = [_formats_TElocal]

def get_input_for_DESeq2_TEcount(wildcards):
    """Dynamically determines the input count matrix for DESeq2 TEcount based on the contrast."""
    logger.debug(f"[get_input_for_DESeq2_TEcount] called with wildcards: {wildcards}")
    in_dict = {}
    count_matrix = indir + "/{genome}/{genome}_TEcount.tsv"
    geneIDAnno = config.get('genome',{}).get('references', {}).get(wildcards.genome, {}).get('geneIDAnno')
    gtf = config.get('genome',{}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not gtf or not os.path.exists(gtf):
        logger.error(f"GTF file for genome {wildcards.genome} not found in config or does not exist")
        raise ValueError(f"GTF file for genome {wildcards.genome} not found in config or does not exist")
    gene_map = config.get('genome',{}).get('references', {}).get(wildcards.genome, {}).get('gene_map')
    if not gene_map or not os.path.exists(gene_map):
        logger.error(f"Gene annotation files for genome {wildcards.genome} not found in config or does not exist")
        raise ValueError(f"Gene annotation files for genome {wildcards.genome} not found in config or does not exist")
    in_dict['count_matrix'] = count_matrix
    in_dict['geneIDAnno'] = geneIDAnno
    in_dict['gtf'] = gtf
    in_dict['gene_map'] = gene_map
    return in_dict

rule DESeq2_TEcount:
    input:
        unpack(get_input_for_DESeq2_TEcount)
    output:
        deseq2_results = outdir + "/{genome}/{contrast}/TEcount/DESeq2.done",
        contrast_group = outdir + "/{genome}/{contrast}/TEcount/group.tsv",
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        contrast_pca = [outdir + "/{genome}/{contrast}/TEcount/PCA/{contrast}.cpmPCA." + fmt for fmt in _formats_TEcount],
        contrast_gene_volcano = [outdir + "/{genome}/{contrast}/TEcount/volcano/{contrast}.TEcount_Gene_volcano." + fmt for fmt in _formats_TEcount],
        contrast_gene_heatmap = [outdir + "/{genome}/{contrast}/TEcount/heatmap/{contrast}.TEcount_Gene_updown." + fmt for fmt in _formats_TEcount],
        contrast_gene_updown = outdir + "/{genome}/{contrast}/TEcount/upDown/{contrast}.TEcount_Gene_updown.tsv",
        contrast_te_updown = outdir + "/{genome}/{contrast}/TEcount/upDown/{contrast}.TEcount_TE_updown.tsv",
        contrast_gene_te_updown = outdir + "/{genome}/{contrast}/TEcount/upDown/{contrast}.TEcount_Gene_TE_updown.tsv",
        contrast_gene_name = outdir + "/{genome}/{contrast}/TEcount/{contrast}.TEcount_Gene.name.tsv",
        contrast_te_name = outdir + "/{genome}/{contrast}/TEcount/{contrast}.TEcount_TE.name.tsv",
        contrast_gene_te_name = outdir + "/{genome}/{contrast}/TEcount/{contrast}.TEcount_Gene_TE.name.tsv",
        contrast_gene = outdir + "/{genome}/{contrast}/TEcount/{contrast}.TEcount_Gene.tsv",
        contrast_te = outdir + "/{genome}/{contrast}/TEcount/{contrast}.TEcount_TE.tsv",
        contrast_gene_te = outdir + "/{genome}/{contrast}/TEcount/{contrast}.TEcount_Gene_TE.tsv",
    threads: rule_res("DESeq2_TEcount").get("threads", 4)
    benchmark: logdir + "/DESeq2/{genome}/{contrast}.benchmark.tsv",
    resources:
        mem_mb=lambda wildcards, input, attempt: rule_res(
            "DESeq2_TEcount", input_size=sample_count(),
            attempt=attempt).get("mem_mb", 16000),
    params:
        DESeq2_script = ROOT_DIR + "/modules/DESeq2/bin/DESeq2.r",
        write_group_script = ROOT_DIR + "/modules/DESeq2/bin/write_group_tsv.py",
        annotation_script = ROOT_DIR + "/modules/DESeq2/bin/gene_id2name.py",
        formats = _formats_TEcount,   # ③ 转发
    conda:
        "../DESeq2.yaml"
    container:
        sif("../DESeq2.yaml", "DESeq2_TEcount")
    log:
        logdir + "/DESeq2/{genome}/{contrast}.log"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("DESeq2_TEcount", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DESeq2 TEcount at {current_time}")
            sample_outdir = os.path.dirname(str(output.deseq2_results))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"DESeq2_TEcount_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"DESeq2_TEcount.sh")
            control_group_name = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name")
            experimental_group_name = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name")
            if wildcards.contrast not in group_pairs.get(wildcards.genome, {}):
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            control_samples = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_samples", [])
            experimental_samples = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_samples", [])
            if not control_samples or len(control_samples) == 0:
                raise ValueError(f"No control samples found for contrast {wildcards.contrast} in genome {wildcards.genome}.")
            if not experimental_samples or len(experimental_samples) == 0:
                raise ValueError(f"No experimental samples found for contrast {wildcards.contrast} in genome {wildcards.genome}.")
            rule_logger.info(f"Control samples: {control_samples}")
            rule_logger.info(f"Experimental samples: {experimental_samples}")
            cmd1 = [
                "python", params.write_group_script,
                "-o", f"{sample_outdir}/group.tsv",
                "-c", ",".join(control_samples),
                "-t", ",".join(experimental_samples),
                "-p", control_group_name,
                "-e", experimental_group_name
            ]
            cmd2 = [
                "Rscript", params.DESeq2_script,
                "-m", "TEcount",
                "-i", input.count_matrix,
                "-g", f"{sample_outdir}/group.tsv",
                "-p", control_group_name, experimental_group_name,
                "-f", "heatmap", "volcano", "pca",
                "-o", sample_outdir,
                "-a", input.geneIDAnno,
                "-Tcm", "all",
                "-px", wildcards.contrast
            ]
            for fmt in params.formats:
                cmd2 += ["--format", str(fmt)]
            cmd3 = [
                "python", params.annotation_script,
                "-i", 
                f"{sample_outdir}/{wildcards.contrast}.TEcount_Gene.tsv",
                f"{sample_outdir}/{wildcards.contrast}.TEcount_TE.tsv",
                f"{sample_outdir}/{wildcards.contrast}.TEcount_Gene_TE.tsv",
                "-g", input.gtf,
                "-m", input.gene_map
            ]

            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(shlex.join(cmd1) + "\n")
                f.write(shlex.join(cmd2) + "\n")
                f.write(shlex.join(cmd3) + "\n")
                f.write("echo 'DESeq2 TEcount analysis completed successfully'\n")

            shutil.copyfile(script, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during DESeq2 TEcount: {e}\n")
            raise RuntimeError(f"Error occurred during DESeq2 TEcount: {e}\n")

def get_input_for_DESeq2_TElocal(wildcards):
    """Dynamically determines the input count matrix for DESeq2 TElocal based on the contrast."""
    logger.debug(f"[get_input_for_DESeq2_TElocal] called with wildcards: {wildcards}")
    in_dict = {}
    count_matrix = indir + "/{genome}/{genome}_TElocal.tsv"
    geneIDAnno = config.get('genome',{}).get('references', {}).get(wildcards.genome, {}).get('geneIDAnno')
    gtf = config.get('genome',{}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not gtf or not os.path.exists(gtf):
        logger.error(f"GTF file for genome {wildcards.genome} not found in config or does not exist")
        raise ValueError(f"GTF file for genome {wildcards.genome} not found in config or does not exist")
    gene_map = config.get('genome',{}).get('references', {}).get(wildcards.genome, {}).get('gene_map')
    if not gene_map or not os.path.exists(gene_map):
        logger.error(f"Gene annotation files for genome {wildcards.genome} not found in config or does not exist")
        raise ValueError(f"Gene annotation files for genome {wildcards.genome} not found in config or does not exist")
    in_dict['count_matrix'] = count_matrix
    in_dict['geneIDAnno'] = geneIDAnno
    in_dict['gtf'] = gtf
    in_dict['gene_map'] = gene_map
    return in_dict

rule DESeq2_TElocal:
    input:
        unpack(get_input_for_DESeq2_TElocal)
    output:
        deseq2_results = outdir + "/{genome}/{contrast}/TElocal/DESeq2.done",
        contrast_group = outdir + "/{genome}/{contrast}/TElocal/group.tsv",
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        contrast_pca = [outdir + "/{genome}/{contrast}/TElocal/PCA/{contrast}.cpmPCA." + fmt for fmt in _formats_TElocal],
        contrast_gene_volcano = [outdir + "/{genome}/{contrast}/TElocal/volcano/{contrast}.TElocal_Gene_volcano." + fmt for fmt in _formats_TElocal],
        contrast_gene_heatmap = [outdir + "/{genome}/{contrast}/TElocal/heatmap/{contrast}.TElocal_Gene_updown." + fmt for fmt in _formats_TElocal],
        contrast_gene_updown = outdir + "/{genome}/{contrast}/TElocal/upDown/{contrast}.TElocal_Gene_updown.tsv",
        contrast_te_updown = outdir + "/{genome}/{contrast}/TElocal/upDown/{contrast}.TElocal_TE_updown.tsv",
        contrast_gene_te_updown = outdir + "/{genome}/{contrast}/TElocal/upDown/{contrast}.TElocal_Gene_TE_updown.tsv",
        contrast_gene_name = outdir + "/{genome}/{contrast}/TElocal/{contrast}.TElocal_Gene.name.tsv",
        contrast_te_name = outdir + "/{genome}/{contrast}/TElocal/{contrast}.TElocal_TE.name.tsv",
        contrast_gene_te_name = outdir + "/{genome}/{contrast}/TElocal/{contrast}.TElocal_Gene_TE.name.tsv",
        contrast_gene = outdir + "/{genome}/{contrast}/TElocal/{contrast}.TElocal_Gene.tsv",
        contrast_te = outdir + "/{genome}/{contrast}/TElocal/{contrast}.TElocal_TE.tsv",
        contrast_gene_te = outdir + "/{genome}/{contrast}/TElocal/{contrast}.TElocal_Gene_TE.tsv",
    threads: rule_res("DESeq2_TElocal").get("threads", 4)
    benchmark: logdir + "/DESeq2/{genome}/{contrast}.benchmark.tsv",
    resources:
        mem_mb=lambda wildcards, input, attempt: rule_res(
            "DESeq2_TElocal", input_size=sample_count(),
            attempt=attempt).get("mem_mb", 16000),
    params:
        DESeq2_script = ROOT_DIR + "/modules/DESeq2/bin/DESeq2.r",
        write_group_script = ROOT_DIR + "/modules/DESeq2/bin/write_group_tsv.py",
        annotation_script = ROOT_DIR + "/modules/DESeq2/bin/gene_id2name.py",
        formats = _formats_TElocal,   # ③ 转发
    conda:
        "../DESeq2.yaml"
    container:
        sif("../DESeq2.yaml", "DESeq2_TElocal")
    log:
        logdir + "/DESeq2/{genome}/{contrast}.log"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("DESeq2_TElocal", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DESeq2 TElocal at {current_time}")
            sample_outdir = os.path.dirname(str(output.deseq2_results))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"DESeq2_TElocal_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"DESeq2_TElocal.sh")
            control_group_name = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name")
            experimental_group_name = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name")
            if wildcards.contrast not in group_pairs.get(wildcards.genome, {}):
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            control_samples = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_samples", [])
            experimental_samples = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_samples", [])
            if not control_samples or len(control_samples) == 0:
                raise ValueError(f"No control samples found for contrast {wildcards.contrast} in genome {wildcards.genome}.")
            if not experimental_samples or len(experimental_samples) == 0:
                raise ValueError(f"No experimental samples found for contrast {wildcards.contrast} in genome {wildcards.genome}.")
            rule_logger.info(f"Control samples: {control_samples}")
            rule_logger.info(f"Experimental samples: {experimental_samples}")
            cmd1 = [
                "python", params.write_group_script,
                "-o", f"{sample_outdir}/group.tsv",
                "-c", ",".join(control_samples),
                "-t", ",".join(experimental_samples),
                "-p", control_group_name,
                "-e", experimental_group_name
            ]
            cmd2 = [
                "Rscript", params.DESeq2_script,
                "-m", "TElocal",
                "-i", input.count_matrix,
                "-g", f"{sample_outdir}/group.tsv",
                "-p", control_group_name, experimental_group_name,
                "-f", "heatmap", "volcano", "pca",
                "-o", sample_outdir,
                "-a", input.geneIDAnno,
                "-Tcm", "all",
                "-px", wildcards.contrast
            ]
            for fmt in params.formats:
                cmd2 += ["--format", str(fmt)]
            cmd3 = [
                "python", params.annotation_script,
                "-i", 
                f"{sample_outdir}/{wildcards.contrast}.TElocal_Gene.tsv",
                f"{sample_outdir}/{wildcards.contrast}.TElocal_TE.tsv",
                f"{sample_outdir}/{wildcards.contrast}.TElocal_Gene_TE.tsv",
                "-g", input.gtf,
                "-m", input.gene_map
            ]

            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(shlex.join(cmd1) + "\n")
                f.write(shlex.join(cmd2) + "\n")
                f.write(shlex.join(cmd3) + "\n")
                f.write("echo 'DESeq2 TElocal analysis completed successfully'\n")

            shutil.copyfile(script, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during DESeq2 TElocal: {e}\n")
            raise RuntimeError(f"Error occurred during DESeq2 TElocal: {e}\n")

# ---------------------------------------------------------------------------
# TE locus figures over TElocal results (Manhattan / locus box / subfamily bar).
# Figure params resolve from Params.DESeq2.<rule_name>.<key> only, via config
# 全量路径(参数精确到 rule,无共享回退).
# Per-contrast family patterns live in Params.DESeq2.<rule_name>.<key>[genome][contrast].
# ---------------------------------------------------------------------------

# ① parse 期读入 + 规范化（str 收成单元素 list；合法值由 schema 声明）
# Output formats for the TE locus figures (matplotlib full set; the te_*.py
# scripts fail fast on anything else). Resolution: config 全量路径
# Params.DESeq2.<rule_name>.formats only (rule-scoped; built-in ["png"]).
_manhattan_formats = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_manhattan", {}).get("formats") or ["png"]
if isinstance(_manhattan_formats, str):
    _manhattan_formats = [_manhattan_formats]
_locus_box_formats = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_locus_boxplot", {}).get("formats") or ["png"]
if isinstance(_locus_box_formats, str):
    _locus_box_formats = [_locus_box_formats]
_subfamily_bar_formats = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_subfamily_bar", {}).get("formats") or ["png"]
if isinstance(_subfamily_bar_formats, str):
    _subfamily_bar_formats = [_subfamily_bar_formats]


def get_TElocal_TEind(wildcards):
    """Resolve the TElocal locus index (.locInd) for per-locus coordinates."""
    teind = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('TEind')
    if not teind:
        raise ValueError(
            f"genome.references.{wildcards.genome}.TEind (TElocal .locInd) not found in config."
        )
    return teind


rule DESeq2_TElocal_manhattan:
    input:
        deseq = outdir + "/{genome}/{contrast}/TElocal/{contrast}.TElocal_TE.tsv",
        group = outdir + "/{genome}/{contrast}/TElocal/group.tsv",
        loci = get_TElocal_TEind,
    output:
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        manhattan = [outdir + "/{genome}/{contrast}/TElocal/manhattan/{contrast}.manhattan." + fmt for fmt in _manhattan_formats],
    threads: rule_res("DESeq2_TElocal_manhattan").get("threads", 2),
    benchmark: logdir + "/DESeq2/{genome}/{contrast}.manhattan.benchmark.tsv",
    resources:
        mem_mb=rule_res("DESeq2_TElocal_manhattan").get("mem_mb", 8000),
    params:
        figure_script = ROOT_DIR + "/modules/DESeq2/bin/te_manhattan.py",
        control_group_name = lambda wildcards: group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name"),
        experimental_group_name = lambda wildcards: group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name"),
        highlight_families = lambda wildcards: list(config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_manhattan", {}).get("highlight_families", {}).get(wildcards.genome, {}).get(wildcards.contrast, [])),
        label_top = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_manhattan", {}).get("label_top", 10),
        ylabel = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_manhattan", {}).get("ylabel", "log2(1+|delta TPM|)"),
        formats = _manhattan_formats,   # ③ 转发
    conda:
        "../tefig.yaml"
    container:
        sif("tefig.yaml", "DESeq2_TElocal_manhattan")
    log:
        logdir + "/DESeq2/{genome}/{contrast}.manhattan.log"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("DESeq2_TElocal_manhattan", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DESeq2 TElocal manhattan at {current_time}")
            fig_outdir = os.path.dirname(str(output.manhattan[0]))
            os.makedirs(fig_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"DESeq2_TElocal_manhattan_{wildcards.genome}_{wildcards.contrast}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(fig_outdir, "DESeq2_TElocal_manhattan.sh")
            if not params.control_group_name or not params.experimental_group_name:
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            cmd = [
                "python", params.figure_script,
                "--deseq", input.deseq,
                "--group", input.group,
                "--control", params.control_group_name,
                "--experiment", params.experimental_group_name,
                "--loci-tsv", input.loci,
                "--signif-col", "padj", "--signif-thresh", "0.05",
                "--value-thresh", "0", "--burst-sign", "up",
                "--y-transform", "log2", "--y-break-auto",
                "--label-top", str(params.label_top),
                "--ylabel", params.ylabel,
                "--hline", "0",
                "-o", fig_outdir,
                "--out-prefix", f"{wildcards.contrast}.manhattan",
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            for family in params.highlight_families:
                cmd += ["--highlight-family", family]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write("echo 'DESeq2 TElocal manhattan completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during DESeq2 TElocal manhattan: {e}\n")
            raise RuntimeError(f"Error occurred during DESeq2 TElocal manhattan: {e}\n")


rule DESeq2_TElocal_locus_boxplot:
    input:
        matrix = indir + "/{genome}/{genome}_TElocal.tsv",
        annotation = outdir + "/{genome}/{contrast}/TElocal/group.tsv",
        loci = get_TElocal_TEind,
    output:
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        locus_box = [outdir + "/{genome}/{contrast}/TElocal/locus_box/{contrast}.locus_box." + fmt for fmt in _locus_box_formats],
    threads: rule_res("DESeq2_TElocal_locus_boxplot").get("threads", 2),
    benchmark: logdir + "/DESeq2/{genome}/{contrast}.locus_box.benchmark.tsv",
    resources:
        mem_mb=rule_res("DESeq2_TElocal_locus_boxplot").get("mem_mb", 8000),
    params:
        figure_script = ROOT_DIR + "/modules/DESeq2/bin/te_locus_boxplot.py",
        control_group_name = lambda wildcards: group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name"),
        experimental_group_name = lambda wildcards: group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name"),
        families = lambda wildcards: list(config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_locus_boxplot", {}).get("families", {}).get(wildcards.genome, {}).get(wildcards.contrast, [])),
        match_field = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_locus_boxplot", {}).get("match_field", "subfamily"),
        match_mode = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_locus_boxplot", {}).get("match_mode", "exact"),
        normalize = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_locus_boxplot", {}).get("normalize", "log2tpm"),
        kind = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_locus_boxplot", {}).get("kind", "box"),
        formats = _locus_box_formats,   # ③ 转发
    conda:
        "../tefig.yaml"
    container:
        sif("tefig.yaml", "DESeq2_TElocal_locus_boxplot")
    log:
        logdir + "/DESeq2/{genome}/{contrast}.locus_box.log"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("DESeq2_TElocal_locus_boxplot", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DESeq2 TElocal locus boxplot at {current_time}")
            fig_outdir = os.path.dirname(str(output.locus_box[0]))
            os.makedirs(fig_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"DESeq2_TElocal_locus_boxplot_{wildcards.genome}_{wildcards.contrast}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(fig_outdir, "DESeq2_TElocal_locus_boxplot.sh")
            if not params.control_group_name or not params.experimental_group_name:
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            if not params.families:
                raise ValueError(
                    f"Params.DESeq2.locus_boxplot.families[{wildcards.genome}][{wildcards.contrast}] "
                    f"is empty; configure the TE family patterns to plot."
                )
            cmd = [
                "python", params.figure_script,
                "--matrix", input.matrix,
                "--annotation", input.annotation,
                "--match-field", params.match_field,
                "--match-mode", params.match_mode,
                "--normalize", params.normalize,
                "--kind", params.kind,
                "--group-order", f"{params.control_group_name},{params.experimental_group_name}",
                "--loci-tsv", input.loci,
                "-o", fig_outdir,
                "--out-prefix", f"{wildcards.contrast}.locus_box",
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            for family in params.families:
                cmd += ["--family", family]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write("echo 'DESeq2 TElocal locus boxplot completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during DESeq2 TElocal locus boxplot: {e}\n")
            raise RuntimeError(f"Error occurred during DESeq2 TElocal locus boxplot: {e}\n")


rule DESeq2_TElocal_subfamily_bar:
    input:
        matrix = indir + "/{genome}/{genome}_TElocal.tsv",
        annotation = outdir + "/{genome}/{contrast}/TElocal/group.tsv",
        loci = get_TElocal_TEind,
    output:
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        subfamily_bar = [outdir + "/{genome}/{contrast}/TElocal/subfamily_bar/{contrast}.subfamily_bar." + fmt for fmt in _subfamily_bar_formats],
    threads: rule_res("DESeq2_TElocal_subfamily_bar").get("threads", 2),
    benchmark: logdir + "/DESeq2/{genome}/{contrast}.subfamily_bar.benchmark.tsv",
    resources:
        mem_mb=rule_res("DESeq2_TElocal_subfamily_bar").get("mem_mb", 8000),
    params:
        figure_script = ROOT_DIR + "/modules/DESeq2/bin/te_subfamily_bar.py",
        control_group_name = lambda wildcards: group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name"),
        experimental_group_name = lambda wildcards: group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name"),
        families = lambda wildcards: list(config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_subfamily_bar", {}).get("families", {}).get(wildcards.genome, {}).get(wildcards.contrast, [])),
        match_field = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_subfamily_bar", {}).get("match_field", "subfamily"),
        match_mode = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_subfamily_bar", {}).get("match_mode", "exact"),
        value = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_subfamily_bar", {}).get("value", "total"),
        normalize = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_subfamily_bar", {}).get("normalize", "tpm"),
        yscale = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_subfamily_bar", {}).get("yscale", "linear"),
        ylabel = config.get("Params", {}).get("DESeq2", {}).get("DESeq2_TElocal_subfamily_bar", {}).get("ylabel", "TPM"),
        formats = _subfamily_bar_formats,   # ③ 转发
    conda:
        "../tefig.yaml"
    container:
        sif("tefig.yaml", "DESeq2_TElocal_subfamily_bar")
    log:
        logdir + "/DESeq2/{genome}/{contrast}.subfamily_bar.log"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("DESeq2_TElocal_subfamily_bar", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DESeq2 TElocal subfamily bar at {current_time}")
            fig_outdir = os.path.dirname(str(output.subfamily_bar[0]))
            os.makedirs(fig_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"DESeq2_TElocal_subfamily_bar_{wildcards.genome}_{wildcards.contrast}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(fig_outdir, "DESeq2_TElocal_subfamily_bar.sh")
            if not params.control_group_name or not params.experimental_group_name:
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            if not params.families:
                raise ValueError(
                    f"Params.DESeq2.subfamily_bar.families[{wildcards.genome}][{wildcards.contrast}] "
                    f"is empty; configure the TE family patterns to plot."
                )
            cmd = [
                "python", params.figure_script,
                "--matrix", input.matrix,
                "--annotation", input.annotation,
                "--match-field", params.match_field,
                "--match-mode", params.match_mode,
                "--value", params.value,
                "--normalize", params.normalize,
                "--yscale", params.yscale,
                "--group-order", f"{params.control_group_name},{params.experimental_group_name}",
                "--loci-tsv", input.loci,
                "--ylabel", params.ylabel,
                "-o", fig_outdir,
                "--out-prefix", f"{wildcards.contrast}.subfamily_bar",
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            for family in params.families:
                cmd += ["--family", family]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write("echo 'DESeq2 TElocal subfamily bar completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during DESeq2 TElocal subfamily bar: {e}\n")
            raise RuntimeError(f"Error occurred during DESeq2 TElocal subfamily bar: {e}\n")