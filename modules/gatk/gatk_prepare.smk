from snakemake.logging import logger
include: "../common/common.smk"
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
fasta = config.get("genome",{}).get("fasta")
input_bam_substring = config.get("input_bam_substring") or ""
rule gatk_index:
    input:
        fasta = fasta
    output:
        fai_index = outdir + "/index/genome.fa.fai",
        dict_index = outdir + "/index/genome.dict",
        fasta_link = outdir + "/index/genome.fa"
    log:
        logdir + "/index/gatk_index.log"
    threads: 4
    params:
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="gatk_index", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start gatk_index at {current_time}")
            script = os.path.join(outdir,f"index/gatk_index_{current_time}.sh")
            cmd1 = ["ln", "-s", input.fasta, output.fasta_link]
            cmd2 = [
                params.gatk, "CreateSequenceDictionary",
                "--java-options", params.javaOptions,
                "-R", input.fasta,
                "-O", output.dict_index
            ]
            if params.tmp_dir:
                cmd2.extend(["--tmp-dir", params.tmp_dir])
            cmd3 = [
                "samtools", "faidx",
                input.fasta,
                "-o", output.fai_index
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during gatk_index execution: {str(e)}\n")
            raise f"Error occurred while running gatk_index: {e}, you can check the log file {log[0]} for more details."


def get_input_for_addReadsGroup(wildcards):
    logger.info(f"addReadsGroup called with {wildcards}")
    in_dict = {}
    if input_bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + input_bam_substring + ".bam")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
    return in_dict

rule addReadsGroup:
    input:
        unpack(get_input_for_addReadsGroup)
    output:
        bam = temp(outdir + "/{sample_id}/{sample_id}.addReadsGroup.bam"),
        bai = temp(outdir + "/{sample_id}/{sample_id}.addReadsGroup.bai")
    log:
        logdir + "/{sample_id}/addReadsGroup.log"
    threads: 8
    params:
        id = "{sample_id}",
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        RGLB = config.get("addReadsGroup", {}).get("RGLB") or "lib1",
        RGPL = config.get("addReadsGroup", {}).get("RGPL") or "illumina",
        RGPU = config.get("addReadsGroup", {}).get("RGPU") or "unit1",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="addReadsGroup", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start addReadsGroup for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir,f"{wildcards.sample_id}/addReadsGroup_{current_time}.sh")
            cmd1 = [
                params.gatk, "AddOrReplaceReadGroups", params.javaOptions,
                "--INPUT", input.bam,
                "--OUTPUT", output.bam,
                "-SO", "coordinate",
                "--RGLB", params.RGLB,
                "--RGPL", params.RGPL,
                "--RGPU", params.RGPU,
                "--RGSM", params.id
            ]
            if params.tmp_dir:
                cmd1.extend(["--tmp-dir", params.tmp_dir])
            cmd2 = [
                params.samtools, "index",
                "-@", str(threads),
                output.bam,
                "-o", output.bai
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
            shell("bash {script} > {log} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during addReadsGroup execution: {str(e)}\n")
            raise f"Error occurred while running addReadsGroup for sample {wildcards.sample_id}: {e}, you can check the log file {log[0]} for more details."


rule MarkDuplicates:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.addReadsGroup.bam",
        bai = outdir + "/{sample_id}/{sample_id}.addReadsGroup.bai"
    output:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bai",
        metrics = outdir + "/{sample_id}/{sample_id}.Markdup-metrics.txt"
    log:
        logdir + "/{sample_id}/MarkDuplicates.log"
    threads: 8
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="MarkDuplicates", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start MarkDuplicates for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir,f"{wildcards.sample_id}/MarkDuplicates_{current_time}.sh")
            cmd = [
                params.gatk, "MarkDuplicates", "--java-options", f'"{params.javaOptions}"',
                "--INPUT", input.bam,
                "--OUTPUT", output.bam,
                "--CREATE_INDEX", "true",
                "--VALIDATION_STRINGENCY", "SILENT",
                "--METRICS_FILE", output.metrics
            ]
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell("bash {script} > {log} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during MarkDuplicates execution: {str(e)}\n")
            raise f"Error occurred while running MarkDuplicates for sample {wildcards.sample_id}: {e}, you can check the log file {log[0]} for more details."


rule gatk_prepare_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bai",
        metrics = outdir + "/{sample_id}/{sample_id}.Markdup-metrics.txt"




