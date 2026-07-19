include: "../common/common.smk"
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
        DESeq2_script = ROOT_DIR + "/modules/DESeq2/bin/DESeq2.r",
        write_group_script = ROOT_DIR + "/modules/DESeq2/bin/write_group_tsv.py",
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
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python {params.write_group_script} \\\n")
                f.write(f"    -o {params.outdir}/group.tsv \\\n")
                f.write(f"    -c {params.control_samples} \\\n")
                f.write(f"    -t {params.treatment_samples} \\\n")
                f.write(f"    -p {params.control_group_name} \\\n")
                f.write(f"    -e {params.experimental_group_name} > {log} 2>&1\n")
                f.write(f"Rscript {params.DESeq2_script} \\\n")
                f.write(f"    -m TEcount \\\n")
                f.write(f"    -i {input.count_matrix} \\\n")
                f.write(f"    -g {params.outdir}/group.tsv \\\n")
                f.write(f"    -p {params.control_group_name} {params.experimental_group_name} \\\n")
                f.write(f"    -f heatmap volcano pca \\\n")
                f.write(f"    -o {params.outdir}/TEcount \\\n")
                f.write(f"    -a {params.geneIDAnno} \\\n")
                f.write(f"    -Tcm all >> {log} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during DESeq2 TEcount: {e}\n")
            logger.error(f"Error occurred during DESeq2 TEcount: {e}")
            raise e