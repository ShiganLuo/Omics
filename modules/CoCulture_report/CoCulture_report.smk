include: "../common/common.smk"
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
    conda:
        "CoCulture_report.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="HRT_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start HRT run at {current_time}")
            script = os.path.join(outdir, f"HRT_{current_time}.sh")
            cmd = ["python", params.HTR_script,
                   "--genome1_count", input.genome1_count,
                   "--genome2_count", input.genome2_count,
                   "--genome1", params.genome1,
                   "--genome2", params.genome2,
                   "--genome1_samples", ' '.join(params.genome1_samples),
                   "--genome2_samples", ' '.join(params.genome2_samples),
                   "--output", output.report]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")
