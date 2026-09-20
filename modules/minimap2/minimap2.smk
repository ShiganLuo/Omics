include: "../common/common.smk"
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")
if not fasta:
    raise ValueError(
        "minimap2 module requires 'genome.fasta' in config. "
        "Please provide a valid reference genome FASTA path."
    )

def get_minimap2_input(wildcards):
    """Dynamically determine input FASTQ for minimap2.

    Checks for single-end (.single.fq.gz) and paired-end (_1.fq.gz) patterns.
    Falls back to a configurable suffix or default .single.fq.gz.
    """
    suffix = config.get("Params", {}).get("minimap2", {}).get("fastq_suffix") or ".single.fq.gz"
    fq_single = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}{suffix}"
    fq_r1 = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
    if os.path.exists(fq_single):
        return {"fastq": fq_single, "fasta": fasta}
    elif os.path.exists(fq_r1):
        return {"fastq": fq_r1, "fasta": fasta}
    else:
        # Default to single-end naming for long-read data
        return {"fastq": fq_single, "fasta": fasta}

rule minimap2_align:
    input:
        unpack(get_minimap2_input)
    output:
        bam = outdir + "/{sample_id}/{sample_id}.sorted.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted.bam.bai"
    log:
        logdir + "/{sample_id}/minimap2_align.log"
    threads: 16
    conda:
        "minimap2.yaml"
    container:
        sif("minimap2.yaml")
    params:
        minimap2 = config.get("Procedure", {}).get("minimap2") or "minimap2",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools",
        preset = config.get("Params", {}).get("minimap2", {}).get("preset") or "map-ont",
        extra = config.get("Params", {}).get("minimap2", {}).get("extra") or ""
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("minimap2_align", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start minimap2 alignment for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.bam))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"minimap2_align_{wildcards.sample_id}_{current_time}.sh")
            raw_bam = outdir + f"/{wildcards.sample_id}/{wildcards.sample_id}.raw.bam"
            # minimap2 align → samtools sort → samtools index
            cmd_align = [
                params.minimap2, "-ax", params.preset,
                "--eqx", "-L",
                "-t", str(threads),
                input.fasta, input.fastq,
                "|", params.samtools, "view", "-@ 4", "-bS", "-",
                "|", params.samtools, "sort", "-@ 4", "-o", raw_bam,
            ]
            cmd_index = [
                params.samtools, "index", "-@ 4", raw_bam,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd_align) + "\n")
                f.write(" ".join(cmd_index) + "\n")
                f.write(f"mv {raw_bam} {output.bam}\n")
                f.write(f"mv {raw_bam}.bai {output.bai}\n")
                f.write(f'echo "minimap2 alignment for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during minimap2 alignment for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during minimap2 alignment for sample {wildcards.sample_id}: {e}")
            raise e
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"minimap2_align finished for sample {wildcards.sample_id} at {current_time}")

rule minimap2_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.sorted.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted.bam.bai"