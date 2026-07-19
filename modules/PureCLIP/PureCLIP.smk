include: "../common/common.smk"
from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")
rule pureclip:
    input:
        bam = indir + "/{sample_id}.dedup.bam",
        bai = indir + "/{sample_id}.dedup.bam.bai",
        fasta = config.get('genome',{}).get('fasta')
    output:
        sites = outdir + "/{sample_id}.pureclip.sites.bed",
        region = outdir + "/{sample_id}.pureclip.region.bed"
    params:
        pureclip = config.get('Procedure',{}).get('pureclip') or 'pureclip',
        ld = '--ld' if config.get('Params',{}).get('pureclip',{}).get('ld', True) else '',
        outdir = outdir
    threads: 8
    conda:
        "PureCLIP.yaml"
    log:
        log = logdir + "/{sample_id}/pureclip_run.txt"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("pureclip", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start pureclip for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.sites))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"pureclip_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.pureclip} {params.ld} -nt {threads} \\\n")
                f.write(f"    -i {input.bam} \\\n")
                f.write(f"    -bai {input.bai} \\\n")
                f.write(f"    -g {input.fasta} \\\n")
                f.write(f"    -o {output.sites} \\\n")
                f.write(f"    -or {output.region} \\\n")
                f.write(f"    > {log} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during pureclip for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during pureclip for sample {wildcards.sample_id}: {e}")
            raise e
