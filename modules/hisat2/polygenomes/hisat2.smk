include: "../../common/common.smk"
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
logdir_combine = config.get("logdir_combine", "log/combine")
genome_paired_samples = config.get('genome_paired_samples', {})
genome_single_samples = config.get('genome_single_samples', {})
def get_input_for_hisat2_index(wildcards):
    """Dynamically determines the input fasta file for hisat2 index based on the genome."""
    logger.info(f"[get_input_for_hisat2_index] called with wildcards: {wildcards}")
    fasta = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('fasta')
    if not fasta:
        logger.error(f"Fasta file for genome {wildcards.genome} not found in config")
        raise ValueError(f"Fasta file for genome {wildcards.genome} not found in config")
    return fasta
rule hisat2_index:
    input:
        fasta = get_input_for_hisat2_index
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
    container:
        sif("../hisat2.yaml")
    params:
        prefix = lambda wildcards: outdir + f"/index/{wildcards.genome}/{wildcards.genome}",
        HISAT2_BUILD = config.get('Procedure', {}).get('hisat2-build') or 'hisat2-build'
    log:
        logdir_combine + "/index/{genome}/hisat2_build.log"
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
                f.write(f"echo 'hisat2_index for genome {wildcards.genome} completed'\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"hisat2_index for genome {wildcards.genome} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"hisat2_index failed for genome {wildcards.genome}: {e}\n")
            raise e

def get_hisat2_index(wildcards):
    logger.debug(f"[get_hisat2_index] called with wildcards: {wildcards}")
    config_index_prefix = config.get('genome', {}).get(wildcards.genome, {}).get('hisat2_index_prefix') or None
    if config_index_prefix:
        first_file = f"{config_index_prefix}.1.ht2"
        if os.path.exists(first_file):
            logger.debug(f"genome {wildcards.genome}'s hisat index exists, use it")
            return [f"{config_index_prefix}.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]
        else:
            logger.debug(f"genome {wildcards.genome}'s hisat index doesn't exist, generate it")
    return [outdir + f"/index/{wildcards.genome}/{wildcards.genome}.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]


def get_alignment_input(wildcards):
    """Dynamically determines the input file type: paired-end or single-end sequencing."""
    logger.debug(f"[get_alignment_input] called with wildcards: {wildcards}")
    paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
    single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"

    if wildcards.sample_id in genome_paired_samples.get(wildcards.genome, []):
        logger.info(f"Paired-end: {[paired_r1, paired_r2]}")
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in genome_single_samples.get(wildcards.genome, []):
        logger.info(f"Single-end: {[single]}")
        return [single]
    else:
        logger.error(f"Sample {wildcards.sample_id} not in paired_samples: {genome_paired_samples.get(wildcards.genome, [])} or single_samples: {genome_single_samples.get(wildcards.genome, [])}")
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
        hisat2 = config.get('Procedure',{}).get('hisat2') or 'hisat2',
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools',
        score_min = config.get('Params',{}).get('hisat2', {}).get('score_min') or "L,0,-0.2",
        no_spliced_alignment = config.get('Params',{}).get('hisat2', {}).get('no-spliced-alignment') or False,
        flag_params = config.get('Params',{}).get('hisat2', {}).get('flag_params') or "",
        k = config.get('Params',{}).get('hisat2', {}).get('k') or 5,
        unmapped_prefix = lambda wildcards: f"{outdir}/{wildcards.genome}/{wildcards.sample_id}/unmapped",
        index_prefix = lambda wildcards, input: input.index[0].rsplit('.', 2)[0],
        input_params = lambda wildcards, input: \
            f"-1 {input.fastq[0]} -2 {input.fastq[1]}" if len(input.fastq) == 2 else f"-U {input.fastq[0]}"
    conda:
        "hisat2.yaml"
    container:
        sif("hisat2.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close
            rule_logger = seup_logger("hisat2_align", log_file=log_path)
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            sample_outdir = os.path.dirname(outfile)
            script = f"{sample_outdir}/hisat2_align.{current_time}.sh"
            cmd1 = [
                f"{params.hisat2}",
                "-x", params.index_prefix,
                "--score-min", params.score_min,
                "-k", str(params.k),
                "--novel-splicesite-outfile", output.splice,
                "--un-conc-gz", params.unmapped_prefix,
                params.flag_params,
                params.input_params,
                "-p", str(threads),
            ]
            if params.no_spliced_alignment:
                cmd1.append("--no-spliced-alignment")
            cmd2 = [
                "|", f"{params.samtools}", "sort", "-@", str(threads), "-o", output.outfile
            ]
            cmd = cmd1 + cmd2
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) +"\n")
                f.write(f"echo 'hisat_align for {wildcards.sample_id} was completed successfully at {current_time}'")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"histat_align failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"histat_align failed for sample {wildcards.sample_id}: {e}")


rule hisat2_result:
    input:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.bam"