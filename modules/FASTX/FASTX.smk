include: "../common/common.smk"

import shutil

from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")

# not test
rule fastx_quality_filter_single:
    input:
        fastq = indir + "/{sample_id}.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}/{sample_id}.fq.gz"),
    params:
        fastx_toolkit = config.get('Procedure',{}).get('fastx_toolkit') or 'fastq_quality_filter',
        q = config.get('Params',{}).get("fastx_toolkit", {}).get('q') or 10,
        p = config.get('Params',{}).get("fastx_toolkit", {}).get('p') or 100,
        Q = config.get('Params',{}).get("fastx_toolkit", {}).get('Q') or 33,
        outdir = outdir
    threads: 2
    conda:
        "FASTX.yaml"
    container:
        sif("FASTX.yaml")
    log:
        logdir + "/{sample_id}/fastx_quality_filter_single.log"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fastx_quality_filter_single", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start fastx_quality_filter_single for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(params.outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fastx_quality_filter_single_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "fastx_quality_filter_single.sh")
            cmd = [
                params.fastx_toolkit,
                "-Q", str(params.Q),
                "-q", str(params.q),
                "-p", str(params.p),
                "-i", str(input.fastq),
                "-o", str(output.fastq),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fastx_quality_filter_single for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            rule_logger.info("Executing: " + " ".join(cmd))
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"fastx_quality_filter_single for sample {wildcards.sample_id} completed")
        except Exception as e:
            rule_logger.error(f"Error occurred during fastx_quality_filter_single for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during fastx_quality_filter_single for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")
