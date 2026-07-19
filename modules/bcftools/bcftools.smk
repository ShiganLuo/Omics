include: "../common/common.smk"
from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule bcftools_sort:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf"
    output:
        vcf = outdir + "/sort/{sample_id}/{sample_id}.sorted.vcf.gz"
    log:
        logdir + "/{sample_id}/bcftools_sort.log"
    conda:
        "bcftools.yaml"
    params:
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("bcftools_sort", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start bcftools sort for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.vcf))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"bcftools_sort_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"mkdir -p $(dirname {output.vcf})\n")
                f.write(f"{params.bcftools} sort \\\n")
                f.write(f"    -O z \\\n")
                f.write(f"    -o {output.vcf} \\\n")
                f.write(f"    {input.vcf} \\\n")
                f.write(f"    > {log} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during bcftools sort for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during bcftools sort for sample {wildcards.sample_id}: {e}")
            raise e

rule bcftools_index_csi:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf.gz"
    output:
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
    log:
        logdir + "/{sample_id}/bcftools_index.log"
    conda:
        "bcftools.yaml"
    params:
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("bcftools_index", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start bcftools index for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.csi))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"bcftools_index_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.bcftools} index {input.vcf} > {log} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during bcftools index for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during bcftools index for sample {wildcards.sample_id}: {e}")
            raise e

rule bcftools_result:
    input:
        vcf = outdir + "/sort/{sample_id}/{sample_id}.sorted.vcf.gz",
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
