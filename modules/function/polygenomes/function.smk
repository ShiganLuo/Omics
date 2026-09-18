include: "../../common/common.smk"

from snakemake.logging import logger
indir = config.get("indir", "output/diff_expression")
counts_dir = config.get("counts_dir", "output/counts")
outdir = config.get("outdir", "output/function")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
group_pairs = config.get("group_pairs", {})
genome = config.get("genome", {})

func_params = config.get("Params", {}).get("function", {})

rule function_go_kegg:
    input:
        deseq2_result = indir + "/{genome}/{contrast}/{contrast}.TEcount_Gene.name.tsv",
    output:
        func_go_plot = outdir + "/{genome}/{contrast}/go_back_to_back.png",
        func_kegg_plot = outdir + "/{genome}/{contrast}/kegg_back_to_back.png",
        func_go_up = outdir + "/{genome}/{contrast}/go_up.csv",
        func_go_down = outdir + "/{genome}/{contrast}/go_down.csv",
        func_kegg_up = outdir + "/{genome}/{contrast}/kegg_up.csv",
        func_kegg_down = outdir + "/{genome}/{contrast}/kegg_down.csv",
        func_up_genes = outdir + "/{genome}/{contrast}/up_genes.txt",
        func_down_genes = outdir + "/{genome}/{contrast}/down_genes.txt"
    log:
        logdir + "/function/{genome}/{contrast}.go_kegg.log"
    threads: 1
    params:
        go_kegg_script = ROOT_DIR + "/modules/function/bin/go-kegg.r",
        species = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("go_kegg", {}).get("species") or "mouse",
        lfc_cut = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("go_kegg", {}).get("lfc_cut") or 1,
        p_cut = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("go_kegg", {}).get("p_cut") or 0.05,
        top = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("go_kegg", {}).get("top") or 10,
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("function_go_kegg", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start GO/KEGG analysis for {wildcards.contrast} at {current_time}")

            if wildcards.contrast not in group_pairs.get(wildcards.genome, {}):
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            ctrl = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name")
            exp = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name")
            sample_outdir = os.path.dirname(str(output.func_go_plot))

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

            script = os.path.join(sample_outdir, f"go_kegg_{current_time}.sh")
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                f.write("set -o pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'GO/KEGG analysis completed successfully at {time.strftime('%Y-%m-%d %H:%M:%S')}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as fh:
                fh.write(f"GO/KEGG analysis failed: {e}\n")
            logger.error(f"GO/KEGG analysis failed: {e}\n")
            raise e

def get_input_for_function_gsea(wildcards):
    logger.info(f"[get_input_for_function_gsea] called with wildcards: {wildcards}")
    in_dict = {}
    deseq2_result = indir + f"/{wildcards.genome}/{wildcards.contrast}/{wildcards.contrast}.TEcount_Gene.tsv"
    annotation = config.get('genome',{}).get('references', {}).get(wildcards.genome, {}).get('geneIDAnno')
    if not annotation or not os.path.exists(annotation):
        raise FileNotFoundError(f"Gene ID annotation not found for genome {wildcards.genome}.")
    gmt = config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("gsea", {}).get("gmt")
    if not gmt or not os.path.exists(gmt):
        raise ValueError(f"GSEA requires a GMT file. Please set Params.function.{wildcards.genome}.gsea.gmt in config.")
    in_dict['deseq2_result'] = deseq2_result
    in_dict['annotation'] = annotation
    in_dict['gmt'] = gmt
    return in_dict

rule function_gsea:
    input:
        unpack(get_input_for_function_gsea)
    output:
        func_gsea_plot = outdir + "/{genome}/{contrast}/GSEA/TEcount_Gene_GSEA.jpeg",
        func_gsea_csv = outdir + "/{genome}/{contrast}/GSEA/TEcount_Gene_GSEA.csv"
    log:
        logdir + "/function/{genome}/{contrast}.gsea.log"
    threads: 1
    params:
        gsea_script = ROOT_DIR + "/modules/function/bin/gsea.r",
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("function_gsea", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            if wildcards.contrast not in group_pairs.get(wildcards.genome, {}):
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            ctrl = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("control_group_name")
            exp = group_pairs.get(wildcards.genome, {}).get(wildcards.contrast, {}).get("experimental_group_name")
            rule_logger.info(f"Start GSEA analysis for {wildcards.contrast} at {current_time}")

            if wildcards.contrast not in group_pairs.get(wildcards.genome, {}):
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            sample_outdir = os.path.dirname(str(output.func_gsea_plot))
            script = os.path.join(sample_outdir, f"gsea_{current_time}.sh")
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

            
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                f.write("set -o pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'GSEA analysis completed successfully at {time.strftime('%Y-%m-%d %H:%M:%S')}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as fh:
                fh.write(f"GSEA analysis failed: {e}\n")
            logger.error(f"GSEA analysis failed: {e}")
            raise e


def get_input_for_function_gsva(wildcards):
    """Resolve inputs for function_gsva: featureCounts counts matrix + GMT file."""
    logger.info(f"[get_input_for_function_gsva] called with wildcards: {wildcards}")
    in_dict = {}
    # featureCounts merged counts (PE + SE combined)
    counts_file = counts_dir + f"/{wildcards.genome}/{wildcards.genome}_featureCounts.tsv"
    in_dict['counts'] = counts_file
    # GMT gene set file
    gmt = config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("gsva", {}).get("gmt")
    if not gmt or not os.path.exists(gmt):
        raise ValueError(
            f"GSVA requires a GMT file. {gmt} is not exists, please set Params.function.{wildcards.genome}.gsva.gmt in config."
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
        gsva_heatmap = outdir + "/{genome}/gsva/gsva_heatmap.png",
    log:
        logdir + "/function/{genome}/gsva.log"
    threads: 1
    params:
        norm_script = ROOT_DIR + "/modules/function/bin/normalization.py",
        gsva_script = ROOT_DIR + "/modules/function/bin/gsva.py",
        species = lambda wildcards: config.get("Params", {}).get("function", {}).get(wildcards.genome, {}).get("go_kegg", {}).get("species") or "mouse",
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("function_gsva", log_file=log_path)
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
            script = os.path.join(gsva_outdir, f"gsva_{current_time}.sh")
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -e\nset -o pipefail\n")
                f.write(" ".join(norm_cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in gsva_cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as fh:
                fh.write(f"GSVA analysis failed: {e}\n")
            logger.error(f"GSVA analysis failed: {e}")
            raise e


rule function_result:
    input:
        go_kegg = outdir + "/{genome}/{contrast}/go_back_to_back.png",
        gsea = outdir + "/{genome}/{contrast}/GSEA/TEcount_Gene_GSEA.jpeg"
