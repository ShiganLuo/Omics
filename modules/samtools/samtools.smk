include: "../common/common.smk"
from snakemake.logging import logger
indir = config.get("indir","input")
outdir = config.get("outdir", "output")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule bam_flagstat:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        flagstat = outdir + "/{sample_id}/{sample_id}_flagstat.txt"
    log:
        logdir + "/{sample_id}/flagstat.log"
    conda:
        "samtools.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("bam_flagstat", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start samtools flagstat for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.flagstat))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"bam_flagstat_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.samtools} flagstat {input.bam} > {output.flagstat}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during samtools flagstat for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during samtools flagstat for sample {wildcards.sample_id}: {e}")
            raise e

rule samtools_result:
    input:
        flagstat = outdir + "/{sample_id}/{sample_id}_flagstat.txt"

