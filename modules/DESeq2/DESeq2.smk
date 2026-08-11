include: "../common/common.smk"
from snakemake.logging import logger
indir = config.get("indir", "data")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
group_pairs = config.get("group_pairs", [])
geneIDAnno = config.get('genome',{}).get('geneIDAnno')
gtf = config.get('genome',{}).get('gtf')
gene_map = config.get('genome',{}).get('gene_map')
rule DESeq2_TEcount:
    input:
        count_matrix = indir + "/TEcount/all_TEcount.tsv",
    output:
        deseq2_results = outdir + "/{contrast}/DESeq2.done",
        contrast_group = outdir + "/{contrast}/group.tsv",
        contrast_pca = outdir + "/{contrast}/PCA/{contrast}.cpmPCA.png",
        contrast_gene_volcano = outdir + "/{contrast}/volcano/{contrast}.TEcount_Gene_volcano.png",
        contrast_gene_heatmap = outdir + "/{contrast}/heatmap/{contrast}.TEcount_Gene_updown.png",
        contrast_gene_updown = outdir + "/{contrast}/upDown/{contrast}.TEcount_Gene_updown.tsv",
        contrast_te_updown = outdir + "/{contrast}/upDown/{contrast}.TEcount_TE_updown.tsv",
        contrast_gene_te_updown = outdir + "/{contrast}/upDown/{contrast}.TEcount_Gene_TE_updown.tsv",
        contrast_gene_name = outdir + "/{contrast}/{contrast}.TEcount_Gene.name.tsv",
        contrast_te_name = outdir + "/{contrast}/{contrast}.TEcount_TE.name.tsv"
    params:
        DESeq2_script = ROOT_DIR + "/modules/DESeq2/bin/DESeq2.r",
        write_group_script = ROOT_DIR + "/modules/DESeq2/bin/write_group_tsv.py",
        annotation_script = ROOT_DIR + "/modules/DESeq2/bin/gene_id2name.py",
        geneIDAnno = geneIDAnno,
        outdir = outdir,
        gtf = gtf,
        gene_map = gene_map

    conda:
        "DESeq2.yaml"
    container:
        sif("DESeq2.yaml")
    log:
        logdir + "/DESeq2/{contrast}.log"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("DESeq2_TEcount", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DESeq2 TEcount at {current_time}")
            sample_outdir = os.path.dirname(str(output.deseq2_results))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"DESeq2_TEcount_{current_time}.sh")
            control_group_name = wildcards.control_group_name
            experimental_group_name = wildcards.experimental_group_name
            if f"{control_group_name}_vs_{experimental_group_name}" not in group_pairs:
                raise ValueError(f"Group pair {control_group_name}_vs_{experimental_group_name} not found in group_pairs configuration.")
            control_samples = group_pairs[f"{control_group_name}_vs_{experimental_group_name}"]["control_samples"]
            experimental_samples = group_pairs[f"{control_group_name}_vs_{experimental_group_name}"]["experimental_samples"]
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
                "-f", "heatmap volcano pca",
                "-o", sample_outdir,
                "-a", params.geneIDAnno,
                "-Tcm", "all",
                "-px", f"{control_group_name}_vs_{experimental_group_name}"
            ]
            cmd3 = [
                "python", params.annotation_script,
                "-i", 
                f"{sample_outdir}/{control_group_name}_vs_{experimental_group_name}.TEcount_Gene.tsv",
                f"{sample_outdir}/{control_group_name}_vs_{experimental_group_name}.TEcount_TE.tsv",
                f"{sample_outdir}/{control_group_name}_vs_{experimental_group_name}.TEcount_Gene_TE.tsv",
                "-g", params.gtf,
                "-m", params.gene_map
            ]

            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                f.write(f"echo 'DESeq2 TEcount analysis completed successfully'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during DESeq2 TEcount: {e}\n")
            logger.error(f"Error occurred during DESeq2 TEcount: {e}\n")
            raise f"Error occurred during DESeq2 TEcount: {e}\n"
