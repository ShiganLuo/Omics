include: "../common/common.smk"

import time

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
abundant_rnas = config.get("Params", {}).get("subsample", {}).get("abundant_rnas", [])
subsample_n = config.get("Params", {}).get("subsample", {}).get("n", 100000)
seed = config.get("Params", {}).get("subsample", {}).get("seed", 42)

def is_abundant(sample_id):
    """Check if sample targets an abundant small RNA."""
    sample_upper = sample_id.upper()
    return any(rna.upper() in sample_upper for rna in abundant_rnas)

def get_subsample_n(wildcards):
    """Return subsample count: 100k for abundant, 0 (all) for others."""
    if is_abundant(wildcards.sample_id):
        return subsample_n
    return 0  # 0 means pass through all reads

rule subsample_fastq:
    input:
        r1 = indir + "/{sample_id}/{sample_id}_1.fq.gz",
        r2 = indir + "/{sample_id}/{sample_id}_2.fq.gz",
    output:
        r1 = outdir + "/{sample_id}/{sample_id}_1.fq.gz",
        r2 = outdir + "/{sample_id}/{sample_id}_2.fq.gz",
    log:
        logdir + "/{sample_id}/subsample.log"
    threads: 1
    conda:
        "subsample.yaml"
    params:
        n = get_subsample_n,
        seed = seed,
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("subsample", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start subsample for sample {wildcards.sample_id} at {current_time}")

            outdir_sample = os.path.dirname(str(output.r1))
            os.makedirs(outdir_sample, exist_ok=True)

            if params.n == 0:
                # Non-abundant: symlink to pass through
                rule_logger.info(f"Sample {wildcards.sample_id} is not abundant, symlinking")
                os.symlink(os.path.abspath(input.r1), str(output.r1))
                os.symlink(os.path.abspath(input.r2), str(output.r2))
            else:
                # Abundant: subsample using seqtk
                rule_logger.info(f"Sample {wildcards.sample_id} is abundant, subsampling to {params.n} reads (seed={params.seed})")

                # Decompress, subsample, recompress
                tmp_r1 = os.path.join(outdir_sample, f"{wildcards.sample_id}_1.tmp.fq")
                tmp_r2 = os.path.join(outdir_sample, f"{wildcards.sample_id}_2.tmp.fq")
                sub_r1 = os.path.join(outdir_sample, f"{wildcards.sample_id}_1.sub.fq")
                sub_r2 = os.path.join(outdir_sample, f"{wildcards.sample_id}_2.sub.fq")

                # Decompress
                shell(f"zcat {input.r1} > {tmp_r1}")
                shell(f"zcat {input.r2} > {tmp_r2}")

                # Subsample with same seed for R1 and R2
                shell(f"seqtk sample -s {params.seed} {tmp_r1} {params.n} > {sub_r1} 2>> {log_path}")
                shell(f"seqtk sample -s {params.seed} {tmp_r2} {params.n} > {sub_r2} 2>> {log_path}")

                # Compress output
                shell(f"gzip -c {sub_r1} > {output.r1}")
                shell(f"gzip -c {sub_r2} > {output.r2}")

                # Cleanup
                for f in [tmp_r1, tmp_r2, sub_r1, sub_r2]:
                    if os.path.exists(f):
                        os.remove(f)

            rule_logger.info(f"subsample completed for sample {wildcards.sample_id}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during subsample for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during subsample for sample {wildcards.sample_id}: {e}")
            raise e
