from snakemake.logging import logger
indir = config.get("indir", "data")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
control_samples = config.get("control_samples", [])
control_group_name = config.get("control_group_name", "control")
treatment_samples = config.get("treatment_samples", [])
experimental_group_name = config.get("experimental_group_name", "treatment")
geneIDAnno = config.get('genome',{}).get('geneIDAnno')
rule DESeq2_TEcount:
    input:
        count_matrix = indir + "/TEcount/all_TEcount.tsv",
    output:
        deseq2_results = directory(outdir + "/TEcount")
    params:
        DESeq2_script = ROOT_DIR + "/DESeq2/bin/DESeq2.r",
        write_group_script = ROOT_DIR + "/DESeq2/bin/write_group_tsv.py",
    conda:
        "DESeq2.yaml"
    log:
        logdir + "/DESeq2/DESeq2.log"
    shell:
        """
        python {params.write_group_script} {outdir}/group.tsv
        Rscript {params.DESeq2_script} \
            -m TEcount \
            -i {input.count_matrix} \
            -g {outdir}/group.tsv \
            -p {control_group_name} {experimental_group_name} \
            -f heatmap volcano pca \
            -o {outdir}/TEcount \
            -a {geneIDAnno} \
            -Tcm all > {log} 2>&1
        """