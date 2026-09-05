include: "../../common/common.smk"
from snakemake.logging import logger
indir = config.get("indir", "data")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
group_pairs = config.get("group_pairs", [])

def get_input_for_DESeq2_TEcount(wildcards):
    """Dynamically determines the input count matrix for DESeq2 TEcount based on the contrast."""
    logger.info(f"[get_input_for_DESeq2_TEcount] called with wildcards: {wildcards}")
    in_dict = {}
    count_matrix = indir + "/{genome}/TEcount/all_TEcount.tsv"
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
        deseq2_results = outdir + "/{genome}/{contrast}/DESeq2.done",
        contrast_group = outdir + "/{genome}/{contrast}/group.tsv",
        contrast_pca = outdir + "/{genome}/{contrast}/PCA/{contrast}.cpmPCA.png",
        contrast_gene_volcano = outdir + "/{genome}/{contrast}/volcano/{contrast}.TEcount_Gene_volcano.png",
        contrast_gene_heatmap = outdir + "/{genome}/{contrast}/heatmap/{contrast}.TEcount_Gene_updown.png",
        contrast_gene_updown = outdir + "/{genome}/{contrast}/upDown/{contrast}.TEcount_Gene_updown.tsv",
        contrast_te_updown = outdir + "/{genome}/{contrast}/upDown/{contrast}.TEcount_TE_updown.tsv",
        contrast_gene_te_updown = outdir + "/{genome}/{contrast}/upDown/{contrast}.TEcount_Gene_TE_updown.tsv",
        contrast_gene_name = outdir + "/{genome}/{contrast}/{contrast}.TEcount_Gene.name.tsv",
        contrast_te_name = outdir + "/{genome}/{contrast}/{contrast}.TEcount_TE.name.tsv",
        contrast_gene_te_name = outdir + "/{genome}/{contrast}/{contrast}.TEcount_Gene_TE.name.tsv",
        contrast_gene = outdir + "/{genome}/{contrast}/{contrast}.TEcount_Gene.tsv",
        contrast_te = outdir + "/{genome}/{contrast}/{contrast}.TEcount_TE.tsv",
        contrast_gene_te = outdir + "/{genome}/{contrast}/{contrast}.TEcount_Gene_TE.tsv",
    params:
        DESeq2_script = ROOT_DIR + "/modules/DESeq2/bin/DESeq2.r",
        write_group_script = ROOT_DIR + "/modules/DESeq2/bin/write_group_tsv.py",
        annotation_script = ROOT_DIR + "/modules/DESeq2/bin/gene_id2name.py",
    conda:
        "../DESeq2.yaml"
    container:
        sif("../DESeq2.yaml")
    log:
        logdir + "/{genome}/DESeq2/{contrast}.log"
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
            control_group_name = group_pairs[wildcards.contrast]["control_group_name"]
            experimental_group_name = group_pairs[wildcards.contrast]["experimental_group_name"]
            if wildcards.contrast not in group_pairs:
                raise ValueError(f"Group pair {wildcards.contrast} not found in group_pairs configuration.")
            control_samples = group_pairs[wildcards.contrast]["control_samples"]
            experimental_samples = group_pairs[wildcards.contrast]["experimental_samples"]
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
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(shlex.join(cmd1) + "\n")
                f.write(shlex.join(cmd2) + "\n")
                f.write(shlex.join(cmd3) + "\n")
                f.write("echo 'DESeq2 TEcount analysis completed successfully'\n")

            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during DESeq2 TEcount: {e}\n")
            logger.error(f"Error occurred during DESeq2 TEcount: {e}\n")
            raise f"Error occurred during DESeq2 TEcount: {e}\n"
