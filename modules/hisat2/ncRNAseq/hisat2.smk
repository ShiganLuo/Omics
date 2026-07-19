include: "../../common/common.smk"

rule hisat2_index_ncRNAseq:
    input:
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output:
        index = expand(
            outdir + "/genome/{{genome}}/index/hista2/{{genome}}.{idx}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
    threads: 8
    conda:
        "../../hisat2.yaml"
    params:
        prefix = lambda wildcards: outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}",
        HISAT2_BUILD = config.get('Procedure', {}).get('hisat2-build') or 'hisat2-build'
    log:
        outdir + "/log/genome/{genome}/hisat2_build.log"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("hisat2_index_ncRNAseq", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start hisat2_index_ncRNAseq for genome {wildcards.genome} at {current_time}")

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

            rule_logger.info(f"hisat2_index_ncRNAseq for genome {wildcards.genome} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"hisat2_index_ncRNAseq failed for genome {wildcards.genome}: {e}\n")
            raise e

def get_hisat2_index_ncRNAseq(wildcards):
    logger.info(f"[get_hisat2_index] called with wildcards: {wildcards}")
    config_index_prefix = config.get('genome', {}).get(wildcards.genome, {}).get('hisat2_index_prefx') or None
    if config_index_prefix:
        first_file = f"{config_index_prefix}.1.ht2"
        if os.path.exists(first_file):
            return [f"{config_index_prefix}.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]
    return expand(
        outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}.{{idx}}.ht2",
        idx = [1, 2, 3, 4, 5, 6, 7, 8]
    )

rule hisat2_align_ncRNAseq_single:
    """
    HISAT2 alignment for single-end ncRNA-seq reads (STAR-equivalent).
    --no-spliced-alignment : disable spliced alignment
    -k 99999 : allow multiple mapping for small RNA
    --score-min L,0,-0.6 : at least 2/3 of read length should match
    """
    input:
        fastq = outdir + "/ncRNAseq/cutadapt/{sample_id}_cutadapt2_trimmed.fq.gz",
        genome_index = get_hisat2_index_ncRNAseq
    output:
        bam = outdir + "/ncRNAseq/bam/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam",
        unmapped = outdir + "/ncRNAseq/bam/{genome}/{sample_id}.unmapped.fq.gz"
    log:
        outdir + "/log/ncRNAseq/HISAT2/{genome}/{sample_id}.log"
    threads: 12
    params:
        HISAT2 = config.get('Procedure', {}).get('hisat2') or 'hisat2',
        SAMTOOLS = config.get('Procedure', {}).get('samtools') or 'samtools',
    conda:
        config['conda']['smallRNAseq']
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("hisat2_align_ncRNAseq_single", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start hisat2_align_ncRNAseq_single for sample {wildcards.sample_id} at {current_time}")

            sample_outdir = os.path.dirname(str(output.bam))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"hisat2_align_{current_time}.sh")

            cmd = [
                params.HISAT2,
                "-x", str(input.genome_index[0]).rsplit('.', 2)[0],
                "-U", str(input.fastq),
                "-p", str(threads),
                "--no-spliced-alignment",
                "-k", "99999",
                "--score-min", "L,0,-0.6",
                "--un-gz", str(output.unmapped),
            ]
            cmd_samtools = [
                params.SAMTOOLS, "sort", "-@", str(threads), "-o", str(output.bam),
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd) + " 2>> " + log_path + " | \\\n")
                f.write(" ".join(cmd_samtools) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"hisat2_align_ncRNAseq_single for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"hisat2_align_ncRNAseq_single failed for sample {wildcards.sample_id}: {e}\n")
            raise e
