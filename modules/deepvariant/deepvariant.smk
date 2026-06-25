from snakemake.logging import logger

include: "../common/common.smk"
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")
fai = config.get("genome", {}).get("fai")
bam_substring = config.get("bam_substring") or ""

def get_input_for_deepvariant_run(wildcards):
    if bam_substring != "":
        bam = os.path.join(indir, wildcards.sample_id, f"{wildcards.sample_id}.{bam_substring}.bam")
        bai = os.path.join(indir, wildcards.sample_id, f"{wildcards.sample_id}.{bam_substring}.bai")
    else:
        bam = os.path.join(indir, wildcards.sample_id, f"{wildcards.sample_id}.bam")
        bai = os.path.join(indir, wildcards.sample_id, f"{wildcards.sample_id}.bai")
    return {"bam": bam, "bai": bai, "fasta": fasta, "fai": fai}

rule deepvariant_run:
    input:
        unpack(get_input_for_deepvariant_run)
    output:
        vcf = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi",
        gvcf = outdir + "/{sample_id}/{sample_id}.g.vcf.gz"
    log:
        logdir + "/{sample_id}/deepvariant.log"
    threads: 8
    container:
        config.get("container", {}).get("deepvariant") or "docker://google/deepvariant:1.10.0"
    params:
        deepvariant = config.get("Procedure", {}).get("deepvariant") or "run_deepvariant",
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools",
        model_type = config.get("Params", {}).get("deepvariant", {}).get("model_type") or "PACBIO",
        outdir_sample = outdir + "/{sample_id}"
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="deepvariant_run", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start DeepVariant for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir,f"{wildcards.sample_id}/deepvariant_{current_time}.sh")
            cmd1 = [
                params.deepvariant,
                "--num_shards", str(threads),
                "--model_type", params.model_type,
                "--ref", input.fasta,
                "--reads", input.bam,
                "--output_vcf", output.vcf,
                "--output_gvcf", output.gvcf
            ]
            cmd2 = [
                params.bcftools, "index", output.vcf
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
            shell("bash {script} > {log} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during DeepVariant execution: {str(e)}\n")
            raise f"Error occurred while running DeepVariant for sample {wildcards.sample_id}: {e}, you can check the log file {log[0]} for more details."

rule deepvariant_result:
    input:
        vcf = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
