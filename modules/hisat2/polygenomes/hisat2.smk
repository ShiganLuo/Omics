include: "../../common/common.smk"

from snakemake.logging import logger

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
paired_samples = config.get('paired_samples', [])
single_samples = config.get('single_samples', [])

rule hisat2_index:
    input:
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output:
        ix1 = outdir + "/index/{genome}/{genome}.1.ht2",
        ix2 = outdir + "/index/{genome}/{genome}.2.ht2",
        ix3 = outdir + "/index/{genome}/{genome}.3.ht2",
        ix4 = outdir + "/index/{genome}/{genome}.4.ht2",
        ix5 = outdir + "/index/{genome}/{genome}.5.ht2",
        ix6 = outdir + "/index/{genome}/{genome}.6.ht2",
        ix7 = outdir + "/index/{genome}/{genome}.7.ht2",
        ix8 = outdir + "/index/{genome}/{genome}.8.ht2"
    threads: 8
    conda:
        "../hisat2.yaml"
    params:
        prefix = lambda wildcards: outdir + f"/index/{wildcards.genome}/{wildcards.genome}",
        HISAT2_BUILD = config.get('Procedure', {}).get('hisat2-build') or 'hisat2-build'
    log:
        logdir + "/index/{genome}/hisat2_build.log"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("hisat2_index", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start hisat2_index for genome {wildcards.genome} at {current_time}")

            sample_outdir = os.path.dirname(params.prefix)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"hisat2_index_{current_time}.sh")

            cmd = [
                params.HISAT2_BUILD, "-p", str(threads),
                input.fasta, params.prefix,
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"hisat2_index for genome {wildcards.genome} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"hisat2_index failed for genome {wildcards.genome}: {e}\n")
            raise e

def get_hisat2_index(wildcards):
    logger.info(f"[get_hisat2_index] called with wildcards: {wildcards}")
    config_index_prefix = config.get('genome', {}).get(wildcards.genome, {}).get('index_prefix') or None
    if config_index_prefix:
        first_file = f"{config_index_prefix}.1.ht2"
        if os.path.exists(first_file):
            logger.info(f"genome {wildcards.genome}'s hisat index exists, use it")
            return [f"{config_index_prefix}.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]
        else:
            logger.info(f"genome {wildcards.genome}'s hisat index doesn't exist, generate it")
    return [outdir + f"/index/{wildcards.genome}/{wildcards.genome}.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]


def get_alignment_input(wildcards):
    """Dynamically determines the input file type: paired-end or single-end sequencing."""
    logger.info(f"[get_alignment_input] called with wildcards: {wildcards}")
    paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
    single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"

    if wildcards.sample_id in paired_samples:
        logger.info(f"Paired-end: {[paired_r1, paired_r2]}")
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        logger.info(f"Single-end: {[single]}")
        return [single]
    else:
        logger.error(f"Sample {wildcards.sample_id} not in paired_samples or single_samples")
        raise ValueError(f"Sample {wildcards.sample_id} not defined in paired_samples or single_samples")

rule hisat2_align:
    input:
        fastq = get_alignment_input,
        index = get_hisat2_index
    output:
        outfile = outdir + "/{genome}/{sample_id}/{sample_id}.bam"
    log:
        logdir + "/{sample_id}/{genome}/hisat2_align.log"
    threads: 12
    conda:
        "../hisat2.yaml"
    params:
        HISAT2 = config.get('Procedure', {}).get('hisat2') or 'hisat2',
        SAMTOOLS = config.get('Procedure', {}).get('samtools') or 'samtools',
        index_prefix = lambda wildcards, input: input.index[0].rsplit('.', 2)[0],
        input_params = lambda wildcards, input:
            f"-1 {input.fastq[0]} -2 {input.fastq[1]}" if len(input.fastq) == 2 else f"-U {input.fastq[0]}"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("hisat2_align", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start hisat2_align for sample {wildcards.sample_id} at {current_time}")

            sample_outdir = os.path.dirname(str(output.outfile))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"hisat2_align_{current_time}.sh")

            # Build hisat2 command
            cmd_hisat2 = [
                params.HISAT2, "-x", params.index_prefix,
                params.input_params,
                "-p", str(threads),
            ]
            # Build samtools sort command
            cmd_samtools = [
                params.SAMTOOLS, "sort", "-@", str(threads), "-o", str(output.outfile),
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd_hisat2) + " 2>> " + log_path + " | \\\n")
                f.write(" ".join(cmd_samtools) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"hisat2_align for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"hisat2_align failed for sample {wildcards.sample_id}: {e}\n")
            raise e

rule hisat2_result:
    input:
        bam = outdir + "/{genome}/{sample_id}.bam"
