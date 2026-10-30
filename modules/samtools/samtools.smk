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
    container:
        sif("samtools.yaml")
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
rule_logger = setup_logger("bam_flagstat", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start samtools flagstat for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.flagstat))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"bam_flagstat_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"bam_flagstat.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.samtools} flagstat {input.bam} > {output.flagstat}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during samtools flagstat for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule samtools_result:
    input:
        flagstat = outdir + "/{sample_id}/{sample_id}_flagstat.txt"

