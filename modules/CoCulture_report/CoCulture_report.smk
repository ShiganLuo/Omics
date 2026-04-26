from snakemake.logging import logger
indir = config.get("indir", "data")
outdir = config.get("outdir", "results")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
genome1_samples = config.get("genome1_samples", [])
genome2_samples = config.get("genome2_samples", [])
genome1, genome2 = config.get("genome_pairs", ["genome1", "genome2"])
# need for test
rule HRT:
    input:
        genome1_count = indir + "/" + genome1 + "/all_TEcount_name.tsv",
        genome2_count = indir + "/" + genome2 + "/all_TEcount_name.tsv"
    output:
        report = outdir + "/CoCulture_report/CoCulture_report.html"
    log:
        logdir + "/CoCulture_report/CoCulture_report.log"
    params:
        genome1 = genome1,
        genome2 = genome2,
        genome1_samples = genome1_samples,
        genome2_samples = genome2_samples,
        HTR_script = ROOT_DIR + "/modules/CoCulture_report/HRT.py"
    shell:
        """
        python {params.HTR_script} \
            --genome1_count {input.genome1_count} \
            --genome2_count {input.genome2_count} \
            --genome1 {params.genome1} \
            --genome2 {params.genome2} \
            --genome1_samples {params.genome1_samples} \
            --genome2_samples {params.genome2_samples} \
            --output {output.report} > {log} 2>&1
        """