include: "../common/common.smk"

from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")

# not test
rule fastx_quality_filter_single:
    input:
        fastq = indir + "/{sample_id}.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}.fq.gz"),
    params:
        fastx_toolkit = config.get('Procedure',{}).get('fastx_toolkit') or 'fastq_quality_filter',
        q = config.get('Params',{}).get("fastx_toolkit", {}).get('q') or 10,
        p = config.get('Params',{}).get("fastx_toolkit", {}).get('p') or 100,
        Q = config.get('Params',{}).get("fastx_toolkit", {}).get('Q') or 33,
        outdir = outdir
    threads: 2
    conda:
        "FASTX.yaml"
    log:
        logdir + "/{sample_id}/fastx.txt"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("fastx_quality_filter_single", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start fastx_quality_filter_single for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(params.outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fastx_quality_filter_single_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.fastx_toolkit} -Q {params.Q} -q {params.q} -p {params.p} \\\n")
                f.write(f"    -i {input.fastq} -o {output.fastq} > {log_path} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"fastx_quality_filter_single for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"fastx_quality_filter_single failed for sample {wildcards.sample_id}: {e}\n")
            raise e
