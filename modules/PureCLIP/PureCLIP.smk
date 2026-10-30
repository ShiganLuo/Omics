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
    container:
        sif("PureCLIP.yaml")
    log:
        log = logdir + "/{sample_id}/pureclip_run.txt"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
rule_logger = setup_logger("pureclip", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start pureclip for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.sites))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"pureclip_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"pureclip.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.pureclip} {params.ld} -nt {threads} \\\n")
                f.write(f"    -i {input.bam} \\\n")
                f.write(f"    -bai {input.bai} \\\n")
                f.write(f"    -g {input.fasta} \\\n")
                f.write(f"    -o {output.sites} \\\n")
                f.write(f"    -or {output.region} \\\n")
                f.write(f"    > {log} 2>&1\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during pureclip for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")
