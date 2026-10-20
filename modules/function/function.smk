include: "../common/common.smk"
from snakemake.logging import logger

indir = config.get("indir", "output/diff_expression")
counts_dir = config.get("counts_dir", "output/counts")
outdir = config.get("outdir", "output/function")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
group_pairs = config.get("group_pairs", {})
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
species = config.get("Params", {}).get("function", {}).get("species", "mouse")



rule function_go_kegg:
    input:
        deseq2_result = indir + "/{contrast}/{contrast}.TEcount_Gene.name.tsv",
    output:
        func_go_plot = outdir + "/{contrast}/go_back_to_back.png",
        func_kegg_plot = outdir + "/{contrast}/kegg_back_to_back.png",
        func_go_up = outdir + "/{contrast}/go_up.csv",
        func_go_down = outdir + "/{contrast}/go_down.csv",
        func_kegg_up = outdir + "/{contrast}/kegg_up.csv",
        func_kegg_down = outdir + "/{contrast}/kegg_down.csv",
        func_up_genes = outdir + "/{contrast}/up_genes.txt",
        func_down_genes = outdir + "/{contrast}/down_genes.txt"
    log:
        logdir + "/../group/function/{contrast}.go_kegg.log"
    threads: 1
    params:
        go_kegg_script = ROOT_DIR + "/modules/function/bin/go-kegg.r",
        species = species,
        lfc_cut = config.get("Params", {}).get("function", {}).get("go_kegg", {}).get("lfc_cut", 1),
        p_cut = config.get("Params", {}).get("function", {}).get("go_kegg", {}).get("p_cut", 0.05),
        top = config.get("Params", {}).get("function", {}).get("go_kegg", {}).get("top", 10),
    conda:
        "function.yaml"
    container:
        sif("function.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("function_go_kegg", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start GO/KEGG analysis for {wildcards.contrast} at {current_time}")

            if wildcards.contrast not in group_pairs:
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            ctrl = group_pairs[wildcards.contrast]["control_group_name"]
            exp = group_pairs[wildcards.contrast]["experimental_group_name"]
            sample_outdir = os.path.join(outdir, wildcards.contrast)
            os.makedirs(sample_outdir, exist_ok=True)

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
                f.write(f"echo 'GO/KEGG analysis completed for {wildcards.contrast}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"GO/KEGG analysis failed: {e}")
            raise RuntimeError(f"GO/KEGG analysis failed: {e}\n")


def get_input_for_function_gsea(wildcards):
    """Resolve inputs for function_gsea: DESeq2 result + GMT + geneIDAnno."""
    logger.debug(f"[get_input_for_function_gsea] called with wildcards: {wildcards}")
    if wildcards.contrast not in group_pairs:
        raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
    in_dict = {}
    deseq2_result = os.path.join(indir, wildcards.contrast, f"{wildcards.contrast}.TEcount_Gene.tsv")
    if not os.path.exists(deseq2_result):
        raise FileNotFoundError(
            f"DESeq2 result not found: {deseq2_result}. "
            f"Ensure DESeq2_TEcount has completed for {wildcards.contrast}."
        )
    in_dict['deseq2_result'] = deseq2_result
    gmt = config.get("Params", {}).get("function", {}).get("gsea", {}).get("gmt")
    if not gmt or not os.path.exists(gmt):
        raise ValueError(
            f"GSEA requires a GMT file. {gmt} is not exists, please set Params.function.gsea.gmt in config."
        )
    in_dict['gmt'] = gmt
    geneIDAnno = config.get("genome",{}).get("references",{}).get("geneIDAnno")
    if not geneIDAnno:
        raise ValueError("GSEA requires geneIDAnno. Please set genome.geneIDAnno in config.")
    in_dict['annotation'] = geneIDAnno
    return in_dict

rule function_gsea:
    input:
        unpack(get_input_for_function_gsea)
    output:
        func_gsea_plot = outdir + "/{contrast}/GSEA/TEcount_Gene_GSEA.jpeg",
        func_gsea_csv = outdir + "/{contrast}/GSEA/TEcount_Gene_GSEA.csv"
    log:
        logdir + "/../group/function/{contrast}.gsea.log"
    threads: 1
    params:
        gsea_script = ROOT_DIR + "/modules/function/bin/gsea.r",
    conda:
        "function.yaml"
    container:
        sif("function.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("function_gsea", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            ctrl = group_pairs[wildcards.contrast]["control_group_name"]
            exp = group_pairs[wildcards.contrast]["experimental_group_name"]
            rule_logger.info(f"Start GSEA analysis for {wildcards.contrast} at {current_time}")

            if not os.path.exists(input.deseq2_result):
                raise FileNotFoundError(
                    f"DESeq2 result not found: {input.deseq2_result}. "
                    f"Ensure DESeq2_TEcount has completed for {wildcards.contrast}."
                )

            sample_outdir = os.path.join(outdir, wildcards.contrast)
            os.makedirs(sample_outdir, exist_ok=True)

            graph_title = wildcards.contrast

            cmd = [
                "Rscript",
                params.gsea_script,
                "-m", "Gene",
                "-g", input.gmt,
                "-i", input.deseq2_result,
                "-o", sample_outdir,
                "-a", input.annotation,
                "-t", graph_title,
            ]

            script = os.path.join(sample_outdir, f"gsea_{current_time}.sh")
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'GSEA analysis completed for {wildcards.contrast}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"GSEA analysis failed: {e}")
            raise RuntimeError(f"GSEA analysis failed: {e}")


def get_input_for_function_gsva(wildcards):
    """Resolve inputs for function_gsva: featureCounts counts matrix + GMT + GTF."""
    logger.debug(f"[get_input_for_function_gsva] called with wildcards: {wildcards}")
    in_dict = {}
    # featureCounts cross-sample counts, one matrix per existing layout
    counts = []
    if paired_samples:
        counts.append(counts_dir + "/all_paired_featureCounts.tsv")
    if single_samples:
        counts.append(counts_dir + "/all_single_featureCounts.tsv")
    if not counts:
        raise ValueError(
            "GSVA requires featureCounts counts. "
            "Please set paired_samples/single_samples in config."
        )
    in_dict['counts'] = counts
    # GMT gene set file
    gsva_gmt = config.get("Params", {}).get("function", {}).get("gsva", {}).get("gmt")
    if not gsva_gmt or not os.path.exists(gsva_gmt):
        raise ValueError(
            f"GSVA requires a GMT file. {gsva_gmt} is not exists, please set Params.function.gsva.gmt in config."
        )
    in_dict['gmt'] = gsva_gmt
    # GTF for gene ID → gene name conversion
    gtf = config.get("genome",{}).get("gtf")
    if not gtf:
        raise ValueError("GSVA requires a GTF. Please set genome.gtf in config.")
    in_dict['gtf'] = gtf
    return in_dict


rule function_gsva:
    """GSVA enrichment analysis on featureCounts normalized expression.

    Pipeline:
    1. featureCounts raw counts (PE + SE merged) → normalization.py (TPM) → normalized matrix
    2. Normalized matrix + GMT → gsva.py → GSVA scores + heatmap

    Inputs:
        - counts: featureCounts cross-sample matrices (PE and/or SE)
        - GMT gene set file
        - GTF annotation for gene ID conversion

    Outputs:
        - gsva_scores.tsv: enrichment score matrix (geneset × samples)
        - gsva_heatmap.png: clustered heatmap of enrichment scores
    """
    input:
        unpack(get_input_for_function_gsva)
    output:
        gsva_scores = outdir + "/gsva/gsva_scores.tsv",
        gsva_heatmap = outdir + "/gsva/gsva_heatmap.png",
    log:
        logdir + "/function/gsva.log"
    threads: 1
    params:
        norm_script = ROOT_DIR + "/modules/function/bin/normalization.py",
        gsva_script = ROOT_DIR + "/modules/function/bin/gsva.py",
        species = species,
    conda:
        "function.yaml"
    container:
        sif("function.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("function_gsva", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start GSVA analysis at {current_time}")

            gsva_outdir = os.path.dirname(str(output.gsva_scores))
            os.makedirs(gsva_outdir, exist_ok=True)

            # Step 0: merge PE/SE featureCounts matrices when both layouts exist
            counts_list = list(input.counts)
            if len(counts_list) > 1:
                import re
                import pandas as pd
                counts_file = os.path.join(gsva_outdir, "all_featureCounts.tsv")
                df_pe = pd.read_csv(counts_list[0], sep="\t", comment='#')
                df_pe.columns = [re.sub(r'.*\/([^\/]+)\.[^.]+$', r'\1', c) if '/' in c else c for c in df_pe.columns]
                df_se = pd.read_csv(counts_list[1], sep="\t", comment='#')
                df_se.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'], inplace=True, errors='ignore')
                df_se.columns = [re.sub(r'.*\/([^\/]+)\.[^.]+$', r'\1', c) if '/' in c else c for c in df_se.columns]
                df_merged = pd.merge(df_pe, df_se, on='Geneid')
                df_merged.to_csv(counts_file, sep="\t", index=False)
                rule_logger.info(f"PE+SE merged: {df_merged.shape[0]} genes × {df_merged.shape[1] - 6} samples")
            else:
                counts_file = counts_list[0]

            # Step 1: Normalize featureCounts → TPM matrix
            tpm_matrix = os.path.join(gsva_outdir, "expression_tpm.tsv")
            norm_cmd = [
                "python", params.norm_script,
                "-i", counts_file,
                "-o", tpm_matrix,
                "--method", "tpm",
                "--gtf", input.gtf,
                "--convert-gene-name",
            ]
            # Step 2: Run GSVA on normalized matrix
            gsva_cmd = [
                "python", params.gsva_script,
                "-i", tpm_matrix,
                "-g", input.gmt,
                "-o", gsva_outdir,
                "--title", f"GSVA ({params.species})",
            ]
            script = os.path.join(gsva_outdir, f"gsva_{current_time}.sh")
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -e\nset -o pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in norm_cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in gsva_cmd) + "\n")
                f.write(f"echo 'GSVA analysis completed for {params.species}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"GSVA analysis failed: {e}")
            raise RuntimeError(f"GSVA analysis failed: {e}")
