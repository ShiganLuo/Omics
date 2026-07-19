include: "../common/common.smk"
from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
# first col: target(human) genome,second col: contaminating genome. human sample may contaminated by mouse genome

def get_inputFile_for_XenofilterR(wildcards):
    logger.info(f"[get_inputFile_for_XenofilterR] called with wildcards: {wildcards}")
    row = [
        f"{indir}/{wildcards.pollution_source_genome}/{wildcards.sample_id}.bam",
        f"{indir}/{wildcards.host_genome}/{wildcards.sample_id}.bam"    
        ]
    return row

rule XenofilterR:
    input:
        bams = get_inputFile_for_XenofilterR
    output:
        csvIn = outdir + "/xenofilterR/{sample_id}/{sample_id}.csv",
        outBam = temp(outdir + "/xenofilterR/{sample_id}/{sample_id}_Filtered.bam"),
        outBai = temp(outdir + "/xenofilterR/{sample_id}/{sample_id}_Filtered.bam.bai")
    log:
        outdir + "/log/XenofilterR/{sample_id}/XenofilterR.log"
    threads: 8
    params:
        csv_content = lambda wildcards, input: ",".join(input.bams),
        outdir = lambda wildcards: f"{outdir}/xenofilterR/{wildcards.sample_id}",
        outSampleName = lambda wildcards: wildcards.sample_id,
        tempBam = lambda wildcards: f"{outdir}/xenofilterR/{wildcards.host_genome}/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam",
        tempBai = lambda wildcards: f"{outdir}/xenofilterR/{wildcards.host_genome}/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam.bai",
        MM = 8,
        script = os.path.join(ROOT_DIR, "modules", "XenofilteR", "utils", "XenofilteR.r"),xenofilterR/{wildcards.sample_id}/Filtered_bams
        Rscript = config.get('Procedure',{}).get('Rscript') or 'Rscript'
    conda:
        "XenofilterR.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("XenofilterR", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start XenofilterR for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.csvIn))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"XenofilterR_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"echo \"{params.csv_content}\" > {output.csvIn}\n")
                f.write(f"# rename ignorme .bam\n")
                f.write(f"{params.Rscript} {params.script} \\\n")
                f.write(f"    --inputFile {output.csvIn} \\\n")
                f.write(f"    --outputDir {params.outdir} \\\n")
                f.write(f"    --renameSamples {params.outSampleName} \\\n")
                f.write(f"    --MM {params.MM} \\\n")
                f.write(f"    --workers 1 > {log} 2>&1\n")
                f.write(f"mv {params.tempBam} {output.outBam} # XenofilteR would run failed if it find Filtered_bams dir exist\n")
                f.write(f"mv {params.tempBai} {output.outBai}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during XenofilterR for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during XenofilterR for sample {wildcards.sample_id}: {e}")
            raise e

