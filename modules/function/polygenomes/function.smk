include: "../../common/common.smk"

indir = config.get("indir", "output/diff_expression")
counts_dir = config.get("counts_dir", "output/counts")
outdir = config.get("outdir", "output/function")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
group_pairs = config.get("group_pairs", {})
genome = config.get("genome", {})

# 参数精确到 rule:参数挂在消费它的 rule 名下
# (Params.function.<genome>.<rule_name>.<key>;原 go_kegg/gsea/gsva 嵌套组名
# 拍平为 rule 名,叶子名作 key;无 module 级共享,读点写 config 全量路径字面链)。
# formats 为每 rule 一个 parse 期变量——output 不能是函数,格式列表必须是
# 模块级变量(decoil.smk 同款注释),各 rule 独立展开 output,因此必须跨
# genome 一致:每个 genome 段各自声明,取值不一致时 parse 期直接报错,
# 而不是悄悄少产声明的输出文件。
# 存量例外(2026-11 返工规范任务一):跨 genome 一致性校验保留函数形态,
# 内部读取为 config 全量路径。
def _uniform_formats(rule_name, key, default):
    """Formats for one rule's key, uniform across genome param blocks.

    A rule sub-section without an explicit list inherits the one that is
    set; more than one distinct explicit list is a config error (rule
    output lists cannot vary per genome wildcard).

    Parameters
    ----------
    rule_name : str
        The consuming rule's name (e.g. function_go_kegg / function_result).
    key : str
        Parameter name inside the rule's sub-section (e.g. formats).
    default : list
        Fallback list (= the extension the rules emit today).

    Returns
    -------
    list
        The single formats list used for both output expansion and params.
    """
    found = {}
    for genome, genome_params in (config.get("Params", {}).get("function", {}) or {}).items():
        if not isinstance(genome_params, dict):
            continue
        fmts = config.get("Params", {}).get("function", {}).get(genome, {}).get(rule_name, {}).get(key)
        if fmts is None:
            continue
        if isinstance(fmts, str):
            fmts = [fmts]
        found[genome] = list(fmts)
    if not found:
        return list(default)
    distinct = {tuple(v) for v in found.values()}
    if len(distinct) > 1:
        raise ValueError(
            f"Params.function.<genome>.{rule_name}.{key} must be identical "
            f"across genomes (rule outputs are expanded at parse time): "
            f"got {found}")
    return list(next(iter(distinct)))

_go_kegg_formats = _uniform_formats("function_go_kegg", "formats", ["png"])
_gsea_formats = _uniform_formats("function_gsea", "formats", ["jpeg"])
_gsva_formats = _uniform_formats("function_gsva", "formats", ["png"])
# function_result 消费 function_go_kegg/function_gsea 的图(跨 rule 消费,
# 变体期望挂在消费 rule 名下;每 rule 一个 parse 期变量)
_result_go_fmts = _uniform_formats("function_result", "go_kegg_formats", ["png"])
_result_gsea_fmts = _uniform_formats("function_result", "gsea_formats", ["jpeg"])


def _raster_pick(formats):
    """Report/consumer variant for a figure: png first, else first raster.

    Cross-module consumers (report embedding) resolve one concrete file;
    follow the raster-first convention from modules.md #9 so the consumer
    never asks for a variant the producer did not write.

    Parameters
    ----------
    formats : list
        Producer formats list.

    Returns
    -------
    str
        The variant the consumer should reference.
    """
    if "png" in formats:
        return "png"
    for fmt in ("jpg", "jpeg", "tif", "tiff"):
        if fmt in formats:
            return fmt
    return formats[0]

rule function_go_kegg:
    input:
        deseq2_result = indir + "/{genome}/{contrast}/TEcount/{contrast}.TEcount_Gene.name.tsv",
    output:
        # ② 按格式展开:每种格式一个输出文件(base 不带扩展名)
        func_go_plot = [outdir + "/{genome}/{contrast}/go_back_to_back." + fmt
                        for fmt in _go_kegg_formats],
        func_kegg_plot = [outdir + "/{genome}/{contrast}/kegg_back_to_back." + fmt
                          for fmt in _go_kegg_formats],
        func_go_up = outdir + "/{genome}/{contrast}/go_up.csv",
        func_go_down = outdir + "/{genome}/{contrast}/go_down.csv",
        func_kegg_up = outdir + "/{genome}/{contrast}/kegg_up.csv",
        func_kegg_down = outdir + "/{genome}/{contrast}/kegg_down.csv",
        func_up_genes = outdir + "/{genome}/{contrast}/up_genes.txt",
        func_down_genes = outdir + "/{genome}/{contrast}/down_genes.txt"
    log:
        logdir + "/function/{genome}/{contrast}.go_kegg.log"
    threads: rule_res("function_go_kegg").get("threads", 1)
    benchmark: logdir + "/function/{genome}/{contrast}.go_kegg.benchmark.tsv",
    resources:
        mem_mb=rule_res("function_go_kegg").get("mem_mb", 8000),
    params:
        go_kegg_script = ROOT_DIR + "/modules/function/bin/go-kegg.r",
        species = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("function_go_kegg", {}).get("species", "mouse") or "mouse",
        lfc_cut = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("function_go_kegg", {}).get("lfc_cut", 1) or 1,
        p_cut = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("function_go_kegg", {}).get("p_cut", 0.05) or 0.05,
        top = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("function_go_kegg", {}).get("top", 10) or 10,
        formats = _go_kegg_formats,   # ③ 转发(跨 genome 一致,见 _uniform_formats)
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml", "function_go_kegg")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("function_go_kegg", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start GO/KEGG analysis for {wildcards.contrast} at {current_time}")

            if wildcards.contrast not in group_pairs.get(wildcards.genome, {}):
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            ctrl = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name")
            exp = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name")
            sample_outdir = os.path.dirname(str(output.func_go_plot[0]))

            if not os.path.exists(input.deseq2_result):
                raise FileNotFoundError(
                    f"DESeq2 result not found: {input.deseq2_result}. "
                    f"Ensure DESeq2_TEcount has completed for {wildcards.contrast}."
                )

            cmd = [
                "Rscript",
                params.go_kegg_script,
                "-i", input.deseq2_result,
                "-o", sample_outdir,
                "-s", params.species,
                "--lfc-cut", str(params.lfc_cut),
                "--p-cut", str(params.p_cut),
                "--top", str(params.top),
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]

            script = os.path.join(RECEIPT_DIR, f"go_kegg_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"go_kegg.sh")
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'GO/KEGG analysis completed successfully for {wildcards.genome} and {wildcards.contrast}'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"GO/KEGG analysis failed: {e}\n")
            raise RuntimeError(f"GO/KEGG analysis failed: {e}\n")

def get_input_for_function_gsea(wildcards):
    logger.debug(f"[get_input_for_function_gsea] called with wildcards: {wildcards}")
    in_dict = {}
    deseq2_result = indir + f"/{wildcards.genome}/{wildcards.contrast}/TEcount/{wildcards.contrast}.TEcount_Gene.tsv"
    annotation = config.get('genome',{}).get('references', {}).get(wildcards.genome, {}).get('geneIDAnno')
    if not annotation or not os.path.exists(annotation):
        raise FileNotFoundError(f"Gene ID annotation not found for genome {wildcards.genome}.")
    gmt = config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("function_gsea", {}).get("gmt", None)
    if not gmt or not os.path.exists(gmt):
        raise ValueError(f"GSEA requires a GMT file. Please set Params.function.{wildcards.genome}.function_gsea.gmt in config.")
    in_dict['deseq2_result'] = deseq2_result
    in_dict['annotation'] = annotation
    in_dict['gmt'] = gmt
    return in_dict

rule function_gsea:
    input:
        unpack(get_input_for_function_gsea)
    output:
        # ② 按格式展开:每种格式一个输出文件(base 不带扩展名)
        func_gsea_plot = [outdir + "/{genome}/{contrast}/GSEA/TEcount_Gene_GSEA." + fmt
                          for fmt in _gsea_formats],
        func_gsea_csv = outdir + "/{genome}/{contrast}/GSEA/TEcount_Gene_GSEA.csv"
    log:
        logdir + "/function/{genome}/{contrast}.gsea.log"
    threads: rule_res("function_gsea").get("threads", 1)
    benchmark: logdir + "/function/{genome}/{contrast}.gsea.benchmark.tsv",
    resources:
        mem_mb=rule_res("function_gsea").get("mem_mb", 8000),
    params:
        gsea_script = ROOT_DIR + "/modules/function/bin/gsea.r",
        formats = _gsea_formats,   # ③ 转发(跨 genome 一致,见 _uniform_formats)
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml", "function_gsea")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("function_gsea", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            if wildcards.contrast not in group_pairs.get(wildcards.genome, {}):
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            ctrl = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name")
            exp = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name")
            rule_logger.info(f"Start GSEA analysis for {wildcards.contrast} at {current_time}")

            if wildcards.contrast not in group_pairs.get(wildcards.genome, {}):
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            sample_outdir = os.path.dirname(str(output.func_gsea_plot[0]))
            script = os.path.join(RECEIPT_DIR, f"gsea_{wildcards.contrast}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"gsea.sh")
            graph_title = wildcards.contrast

            cmd = [
                "Rscript",
                params.gsea_script,
                "-m", "Gene",
                "-g", input.gmt,
                "-i", input.deseq2_result,
                "-o", os.path.dirname(sample_outdir),
                "-a", input.annotation,
                "-t", graph_title,
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]

            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'GSEA analysis was completed successfully for {wildcards.genome} and {wildcards.contrast}'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"GSEA analysis failed: {e}")
            raise RuntimeError(f"GSEA analysis failed: {e}")

def get_input_for_function_gsva(wildcards):
    """Resolve inputs for function_gsva: featureCounts counts matrix + GMT file."""
    logger.debug(f"[get_input_for_function_gsva] called with wildcards: {wildcards}")
    in_dict = {}
    # featureCounts merged counts (PE + SE combined)
    counts_file = counts_dir + f"/{wildcards.genome}/{wildcards.genome}_featureCounts.tsv"
    in_dict['counts'] = counts_file
    # GMT gene set file
    gmt = config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("function_gsva", {}).get("gmt", None)
    if not gmt or not os.path.exists(gmt):
        raise ValueError(
            f"GSVA requires a GMT file. {gmt} is not exists, please set Params.function.{wildcards.genome}.function_gsva.gmt in config."
        )
    in_dict['gmt'] = gmt
    # GTF for gene ID → gene name conversion
    gtf = genome.get("references", {}).get(wildcards.genome, {}).get("gtf")
    if not gtf:
        raise ValueError(f"GTF not found for genome {wildcards.genome} in config.genome.references.")
    in_dict['gtf'] = gtf
    return in_dict

rule function_gsva:
    """GSVA enrichment analysis on featureCounts normalized expression.

    Pipeline:
    1. featureCounts raw counts → normalization.py (TPM) → normalized matrix
    2. Normalized matrix + GMT → gsva.py → GSVA scores + heatmap

    Inputs:
        - counts: merged featureCounts output (PE + SE)
        - GMT gene set file
        - GTF annotation for gene ID conversion

    Outputs:
        - gsva_scores.tsv: enrichment score matrix (geneset × samples)
        - gsva_heatmap.png: clustered heatmap of enrichment scores
    """
    input:
        unpack(get_input_for_function_gsva)
    output:
        gsva_scores = outdir + "/{genome}/gsva/gsva_scores.tsv",
        # ② 按格式展开:每种格式一个输出文件(base 不带扩展名)
        gsva_heatmap = [outdir + "/{genome}/gsva/gsva_heatmap." + fmt
                        for fmt in _gsva_formats],
    log:
        logdir + "/function/{genome}/gsva.log"
    threads: rule_res("function_gsva").get("threads", 1)
    benchmark: logdir + "/function/{genome}/gsva.benchmark.tsv",
    resources:
        mem_mb=rule_res("function_gsva").get("mem_mb", 8000),
    params:
        norm_script = ROOT_DIR + "/modules/function/bin/normalization.py",
        gsva_script = ROOT_DIR + "/modules/function/bin/gsva.py",
        species = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("function_gsva", {}).get("species", "mouse") or "mouse",
        formats = _gsva_formats,   # ③ 转发(跨 genome 一致,见 _uniform_formats)
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml", "function_gsva")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("function_gsva", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start GSVA analysis at {current_time}")

            gsva_outdir = os.path.dirname(str(output.gsva_scores))
            os.makedirs(gsva_outdir, exist_ok=True)

            # Step 1: Normalize featureCounts → TPM matrix
            tpm_matrix = os.path.join(gsva_outdir, "expression_tpm.tsv")
            norm_cmd = [
                "python", params.norm_script,
                "-i", input.counts,
                "-o", tpm_matrix,
                "--method", "tpm",
                "--gtf", input.gtf,
                "--convert-gene-name",
            ]
            # Step 2: Run GSVA on normalized matrix
            species_name = params.species
            gsva_cmd = [
                "python", params.gsva_script,
                "-i", tpm_matrix,
                "-g", input.gmt,
                "-o", gsva_outdir,
                "--title", f"GSVA ({species_name}) - {wildcards.genome}",
            ]
            for fmt in params.formats:
                gsva_cmd += ["--format", str(fmt)]
            script = os.path.join(RECEIPT_DIR, f"gsva_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(gsva_outdir, f"gsva.sh")
            with open(script, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(norm_cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in gsva_cmd) + "\n")
                f.write(f"echo 'GSVA analysis was completed successfully for {wildcards.genome} !'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"GSVA analysis failed: {e}")
            raise RuntimeError(f"GSVA analysis failed: {e}")

rule function_result:
    input:
        # 消费端按 raster-first 约定解析同一变体(modules.md #9):
        # png 优先,否则取第一个光栅格式,否则取第一个格式
        go_kegg = outdir + "/{genome}/{contrast}/go_back_to_back." + _raster_pick(_result_go_fmts),
        gsea = outdir + "/{genome}/{contrast}/GSEA/TEcount_Gene_GSEA." + _raster_pick(_result_gsea_fmts)
