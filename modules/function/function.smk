include: "../common/common.smk"
from snakemake.logging import logger

indir = config.get("indir", "output/diff_expression")
outdir = config.get("outdir", "output/function")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
group_pairs = config.get("group_pairs", {})
geneIDAnno = config.get("genome", {}).get("geneIDAnno")

func_params = config.get("Params", {}).get("function", {})
species = func_params.get("species", "mouse")
gmt = func_params.get("gmt")
lfc_cut = func_params.get("lfc_cut", 1)
p_cut = func_params.get("p_cut", 0.05)
top = func_params.get("top", 10)


rule function_go_kegg:
    input:
        deseq2_done = indir + "/{control_group_name}_vs_{experimental_group_name}/DESeq2.done",
    output:
        go_png = outdir + "/{control_group_name}_vs_{experimental_group_name}/go_back_to_back.png",
        kegg_png = outdir + "/{control_group_name}_vs_{experimental_group_name}/kegg_back_to_back.png",
    log:
        logdir + "/function/{control_group_name}_vs_{experimental_group_name}.go_kegg.log"
    threads: 1
    params:
        go_kegg_script = ROOT_DIR + "/modules/function/bin/go-kegg.r",
        deseq2_result = indir + "/{control_group_name}_vs_{experimental_group_name}/TEcount_Gene.name.tsv",
        species = species,
        lfc_cut = lfc_cut,
        p_cut = p_cut,
        top = top,
    conda:
        "function.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("function_go_kegg", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            ctrl = wildcards.control_group_name
            exp = wildcards.experimental_group_name
            pair_key = f"{ctrl}_vs_{exp}"
            rule_logger.info(f"Start GO/KEGG analysis for {pair_key} at {current_time}")

            if pair_key not in group_pairs:
                raise ValueError(f"Group pair {pair_key} not found in group_pairs configuration.")

            sample_outdir = os.path.join(outdir, pair_key)
            os.makedirs(sample_outdir, exist_ok=True)

            if not os.path.exists(params.deseq2_result):
                raise FileNotFoundError(
                    f"DESeq2 result not found: {params.deseq2_result}. "
                    f"Ensure DESeq2_TEcount has completed for {pair_key}."
                )

            cmd = [
                "Rscript",
                params.go_kegg_script,
                "-i", params.deseq2_result,
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
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error in GO/KEGG analysis: {e}")
            raise e


rule function_gsea:
    input:
        deseq2_done = indir + "/{control_group_name}_vs_{experimental_group_name}/DESeq2.done",
    output:
        gsea_png = outdir + "/{control_group_name}_vs_{experimental_group_name}/GSEA/TEcount_Gene_GSEA.jpeg",
    log:
        logdir + "/function/{control_group_name}_vs_{experimental_group_name}.gsea.log"
    threads: 1
    params:
        gsea_script = ROOT_DIR + "/modules/function/bin/gsea.r",
        deseq2_result = indir + "/{control_group_name}_vs_{experimental_group_name}/TEcount_Gene.tsv",
        annotation = geneIDAnno,
        gmt = gmt,
    conda:
        "function.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("function_gsea", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            ctrl = wildcards.control_group_name
            exp = wildcards.experimental_group_name
            pair_key = f"{ctrl}_vs_{exp}"
            rule_logger.info(f"Start GSEA analysis for {pair_key} at {current_time}")

            if pair_key not in group_pairs:
                raise ValueError(f"Group pair {pair_key} not found in group_pairs configuration.")

            if not params.gmt:
                raise ValueError("GSEA requires a GMT file. Please set Params.function.gmt in config.")

            if not params.annotation:
                raise ValueError("GSEA requires geneIDAnno. Please set genome.geneIDAnno in config.")

            if not os.path.exists(params.deseq2_result):
                raise FileNotFoundError(
                    f"DESeq2 result not found: {params.deseq2_result}. "
                    f"Ensure DESeq2_TEcount has completed for {pair_key}."
                )

            sample_outdir = os.path.join(outdir, pair_key)
            os.makedirs(sample_outdir, exist_ok=True)

            graph_title = pair_key

            cmd = [
                "Rscript",
                params.gsea_script,
                "-m", "Gene",
                "-g", params.gmt,
                "-i", params.deseq2_result,
                "-o", sample_outdir,
                "-a", params.annotation,
                "-t", graph_title,
            ]

            script = os.path.join(sample_outdir, f"gsea_{current_time}.sh")
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                f.write("set -o pipefail\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error in GSEA analysis: {e}")
            raise e
