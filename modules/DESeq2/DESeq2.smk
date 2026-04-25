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
        control_group_name = control_group_name,
        experimental_group_name = experimental_group_name,
        control_samples = ','.join(control_samples),
        treatment_samples = ','.join(treatment_samples),
        geneIDAnno = geneIDAnno,
        outdir = outdir
    conda:
        "DESeq2.yaml"
    log:
        logdir + "/DESeq2/DESeq2.log"
    shell:
        """
        python {params.write_group_script} \
            -o {params.outdir}/group.tsv \
            -c {params.control_samples} \
            -t {params.treatment_samples} \
            -p {params.control_group_name} \
            -e {params.experimental_group_name}
        Rscript {params.DESeq2_script} \
            -m TEcount \
            -i {input.count_matrix} \
            -g {params.outdir}/group.tsv \
            -p {params.control_group_name} {params.experimental_group_name} \
            -f heatmap volcano pca \
            -o {params.outdir}/TEcount \
            -a {params.geneIDAnno} \
            -Tcm all > {log} 2>&1
        """