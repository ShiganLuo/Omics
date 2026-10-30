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
    container:
        sif("bcftools.yaml")
    params:
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
rule_logger = setup_logger("bcftools_sort", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start bcftools sort for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.vcf))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"bcftools_sort_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"bcftools_sort.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"mkdir -p $(dirname {output.vcf})\n")
                f.write(f"{params.bcftools} sort \\\n")
                f.write(f"    -O z \\\n")
                f.write(f"    -o {output.vcf} \\\n")
                f.write(f"    {input.vcf} \\\n")
                f.write(f"    > {log} 2>&1\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during bcftools sort for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule bcftools_index_csi:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf.gz"
    output:
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
    log:
        logdir + "/{sample_id}/bcftools_index.log"
    conda:
        "bcftools.yaml"
    container:
        sif("bcftools.yaml")
    params:
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
rule_logger = setup_logger("bcftools_index", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start bcftools index for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.csi))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"bcftools_index_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"bcftools_index.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.bcftools} index {input.vcf} > {log} 2>&1\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during bcftools index for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule bcftools_result:
    input:
        vcf = outdir + "/sort/{sample_id}/{sample_id}.sorted.vcf.gz",
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
