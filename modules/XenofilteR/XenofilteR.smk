include: "../common/common.smk"
from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
# first col: target(human) genome,second col: contaminating genome. human sample may contaminated by mouse genome

def get_input_for_XenofilteR(wildcards):
    logger.info(f"[get_input_for_XenofilteR] called with wildcards: {wildcards}")
    host = config.get('Params', {}).get('XenofilteR', {}).get('sample_contamination', {}).get(wildcards.sample_id, {}).get('host')
    contaminant = config.get('Params', {}).get('XenofilteR', {}).get('sample_contamination', {}).get(wildcards.sample_id, {}).get('contaminant')
    host_bam = indir + f"/{host}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    if host == contaminant:
        contaminant_bam = host_bam
    else:
        contaminant_bam = indir + f"/{contaminant}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    in_dict = {
        "contaminant_bam": contaminant_bam,
        "host_bam": host_bam
    }
    logger.info(f"[get_input_for_XenofilteR] returning input dictionary: {in_dict}")
    return in_dict

rule XenofilteR:
    input:
        unpack(get_input_for_XenofilteR)
    output:
        csvIn = outdir + "/{genome}/{sample_id}/{sample_id}.csv",
        outBam = temp(outdir + "/{genome}/{sample_id}/{sample_id}.bam"),
        outBai = temp(outdir + "/{genome}/{sample_id}/{sample_id}.bam.bai"),
    log:
        logdir + "/{sample_id}/{genome}/XenofilteR.log"
    threads: 8
    params:
        outdir = lambda wildcards: f"{outdir}/{wildcards.sample_id}",
        outSampleName = lambda wildcards: wildcards.sample_id,
        tempBam = lambda wildcards: f"{outdir}/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam",
        tempBai = lambda wildcards: f"{outdir}/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam.bai",
        MM = config.get('Parameters', {}).get('XenofilteR', {}).get('MM', 8),
        script = ROOT_DIR + "/modules/XenofilteR/bin/XenofilteR.r",
        Rscript = config.get('Procedure',{}).get('Rscript') or 'Rscript'
    conda:
        "XenofilteR.yaml"
    container:
        sif("XenofilteR.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("XenofilteR", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start XenofilteR for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.csvIn))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"XenofilteR_{current_time}.sh")
            cmd1 = [
                "echo", f"{input.contaminant_bam},{input.host_bam}", ">", output.csvIn
            ]
            cmd2 = [
                params.Rscript, params.script,
                "--inputFile", output.csvIn,
                "--outputDir", params.outdir,
                "--renameSamples", params.outSampleName,
                "--MM", str(params.MM),
                "--workers", "1"
            ]
            cmd3 = [
                "mv", params.tempBam, output.outBam
            ]
            cmd4 = [
                "mv", params.tempBai, output.outBai
            ]
            cmd5 = [
                "ln", "-s", input.host_bam, output.outBam
            ]
            host_genome = config.get('Params', {}).get('XenofilteR', {}).get('sample_contamination', {}).get(wildcards.sample_id, {}).get('host')
            contaminant_genome = config.get('Params', {}).get('XenofilteR', {}).get('sample_contamination', {}).get(wildcards.sample_id, {}).get('contaminant')
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                if input.contaminant_bam == input.host_bam:
                    rule_logger.info(f"Host genome and pollution source genome are the same for sample {wildcards.sample_id}. Skipping XenofilteR filtering.")
                    f.write(" ".join(cmd5) + "\n")
                    f.write(f"touch {output.csvIn}\n")
                    f.write(f"ln -s {input.host_bam}.bai {output.outBai}\n")
                else:
                    rule_logger.info(f"Running XenofilteR for sample {wildcards.sample_id} with host genome {host_genome} and pollution source genome {contaminant_genome}.")
                    f.write(" ".join(cmd1) + "\n")
                    f.write(" ".join(cmd2) + "\n")
                    f.write(" ".join(cmd3) + "\n")
                    f.write(" ".join(cmd4) + "\n")
                f.write(f"echo 'XenofilteR completed successfully for sample {wildcards.sample_id}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during XenofilteR for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during XenofilteR for sample {wildcards.sample_id}: {e}")
            raise e

