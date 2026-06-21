include: "../../common/common.smk"
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
interval = config.get("reference", {}).get("interval")
known_sites = config.get("reference", {}).get("known_sites", [])


rule BaseRecalibrator:
    input:
        bam = indir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = indir + "/{sample_id}/{sample_id}.sorted_markdup.bai",
        ref = fasta,
        interval = interval
    output:
        table = outdir + "/bqsr/{sample_id}.recal_data.table"
    log:
        logdir + "/{sample_id}/gatk_bqsr.log"
    threads: 8
    conda:
        "../gatk.yaml"
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        known_sites = " ".join([f"--known-sites {ks}" for ks in known_sites])
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="BaseRecalibrator", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start BaseRecalibrator for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir,f"{wildcards.sample_id}/BaseRecalibrator_{current_time}.sh")
            cmd = [
                params.gatk, "BaseRecalibrator",
                "--java-options", params.javaOptions,
                "-R", input.ref,
                "-I", input.bam,
                params.known_sites
            ]
            if input.interval:
                cmd.extend(["-L", input.interval])
            cmd.extend([
                "-O", output.table
            ])
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during BaseRecalibrator execution: {str(e)}\n")
            raise f"Error occurred while running BaseRecalibrator for sample {wildcards.sample_id}: {e}, you can check the log file {log[0]} for more details."



rule ApplyBQSR:
    input:
        bam = indir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        table = outdir + "/{sample_id}/{sample_id}.recal_data.table",
        ref = fasta
    output:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam"
    log:
        logdir + "/{sample_id}/ApplyBQSR.log"
    conda:
        "../gatk.yaml"
    threads:
        8
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        interval = interval
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="ApplyBQSR", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start ApplyBQSR for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir,f"{wildcards.sample_id}/ApplyBQSR_{current_time}.sh")
            cmd = [
                params.gatk, "ApplyBQSR",
                "--java-options", params.javaOptions,
                "-R", input.ref,
                "-I", input.bam,
                "-bqsr", input.table,
                "-O", output.bam
            ]
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            if params.interval:
                cmd.extend(["-L", params.interval])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during ApplyBQSR execution: {str(e)}\n")
            raise f"Error occurred while running ApplyBQSR for sample {wildcards.sample_id}: {e}, you can check the log file {log[0]} for more details."

rule gatk_bqsr_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam.bai"
