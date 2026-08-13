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
        lfc_cut = lfc_cut,
        p_cut = p_cut,
        top = top,
    conda:
        "function.yaml"
    container:
        sif("function.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("function_go_kegg", log_file=log_path)
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
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as fh:
                fh.write(f"GO/KEGG analysis failed: {e}\n")
            raise f"GO/KEGG analysis failed: {e}\n"


rule function_gsea:
    input:
        deseq2_result = indir + "/{contrast}/{contrast}.TEcount_Gene.tsv",
    output:
        func_gsea_plot = outdir + "/{contrast}/GSEA/TEcount_Gene_GSEA.jpeg",
        func_gsea_csv = outdir + "/{contrast}/GSEA/TEcount_Gene_GSEA.csv"
    log:
        logdir + "/../group/function/{contrast}.gsea.log"
    threads: 1
    params:
        gsea_script = ROOT_DIR + "/modules/function/bin/gsea.r",
        annotation = geneIDAnno,
        gmt = gmt,
    conda:
        "function.yaml"
    container:
        sif("function.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("function_gsea", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            ctrl = group_pairs[wildcards.contrast]["control_group_name"]
            exp = group_pairs[wildcards.contrast]["experimental_group_name"]
            rule_logger.info(f"Start GSEA analysis for {wildcards.contrast} at {current_time}")

            if wildcards.contrast not in group_pairs:
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")

            if not params.gmt:
                raise ValueError("GSEA requires a GMT file. Please set Params.function.gmt in config.")

            if not params.annotation:
                raise ValueError("GSEA requires geneIDAnno. Please set genome.geneIDAnno in config.")

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
                "-g", params.gmt,
                "-i", input.deseq2_result,
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
            with open(log_path, "a") as fh:
                fh.write(f"GSEA analysis failed: {e}\n")
            raise f"GSEA analysis failed: {e}\n"
