include: "../common/common.smk"

import time

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
gtf = config.get("genome", {}).get("gtf")
read_num = config.get("Params", {}).get("tailer", {}).get("read", 2)
threshold = config.get("Params", {}).get("tailer", {}).get("threshold", 100)
rev_comp = config.get("Params", {}).get("tailer", {}).get("rev_comp", False)

rule tailer_global:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        bai = indir + "/{sample_id}/{sample_id}.bam.bai",
    output:
        tail = outdir + "/{sample_id}/{sample_id}_tail.csv",
    log:
        logdir + "/{sample_id}/tailer.log"
    threads: 1
    conda:
        "tailer.yaml"
    params:
        gtf = gtf,
        read_num = read_num,
        threshold = threshold,
        rev_comp = "--rev_comp" if rev_comp else "",
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("tailer_global", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start tailer_global for sample {wildcards.sample_id} at {current_time}")

            outdir_tail = os.path.dirname(str(output.tail))
            os.makedirs(outdir_tail, exist_ok=True)

            cmd = f"Tailer -a {params.gtf} -read {params.read_num} -t {params.threshold} {params.rev_comp} {input.bam} >> {log_path} 2>&1"
            rule_logger.info(f"Running: {cmd}")
            shell(cmd)

            rule_logger.info(f"tailer_global completed for sample {wildcards.sample_id}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during tailer_global for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during tailer_global for sample {wildcards.sample_id}: {e}")
            raise e
