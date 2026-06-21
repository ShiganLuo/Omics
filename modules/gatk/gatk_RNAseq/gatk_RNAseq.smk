include: "../../common/common.smk"
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"


rule SplitNCigarReads:
    input:
        bam = indir + "/bam-sorted-Markdup/{sample_id}.bam",
    output:
        bam = outdir + "/Split/{sample_id}/{sample_id}.bam"
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        gatk = config["Procedure"]["gatk"],
        genome = lambda wildcards: config['genome']['fasta'],
        indict = lambda wildcards: config['genome']['dict_index'],
        fai = lambda wildcards: config['genome']['fai_index']
    threads: 8 
    log:
        logdir + "/{sample_id}/SplitNCigarReads.log"
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="SplitNCigarReads", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start SplitNCigarReads for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir,f"{wildcards.sample_id}/SplitNCigarReads_{current_time}.sh")
            cmd = [
                params.gatk, "SplitNCigarReads",
                "--java-options", params.javaOptions,
                "-R", params.genome,
                "-I", input.bam,
                "-O", output.bam
            ]
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during SplitNCigarReads execution: {str(e)}\n")
            raise f"Error occurred while running SplitNCigarReads for sample {wildcards.sample_id}: {e}, you can check the log file {log[0]} for more details."


rule VarientCalling:
    input:
        bam = outdir + "/Split/{sample_id}/{sample_id}.bam"
    output:
        vcf = outdir + "/vcf/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/VarientCalling.log"
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        gatk = config["Procedure"]["gatk"],
        genome = lambda wildcards: config['genome']['fasta'],
        fai = lambda wildcards: config['genome']['fai_index']
    threads: 8
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="VarientCalling", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start VarientCalling for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir,f"{wildcards.sample_id}/VarientCalling_{current_time}.sh")
            cmd = [
                params.gatk, "HaplotypeCaller",
                "--java-options", params.javaOptions,
                "-R", params.genome,
                "-I", input.bam,
                "-O", output.vcf,
                "--dont-use-soft-clipped-bases",
                "--standard-min-confidence-threshold-for-calling", "20"
            ]
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during VarientCalling execution: {str(e)}\n")
            raise f"Error occurred while running VarientCalling for sample {wildcards.sample_id}: {e}, you can check the log file {log[0]} for more details."

rule vcf_filter:
    input:
        vcf = outdir + "/vcf/{sample_id}/{sample_id}.vcf.gz"
    output:
        vcf = outdir + "/vcf-filtered/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/vcf_filter.log"
    threads: 8
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        vcf = outdir + "/SNP/vcf/filter/{genome}/{sample_id}.vcf",
        gatk = config["Procedure"]["gatk"],
        bgzip = config["Procedure"]["bgzip"],
        genome = lambda wildcards: config['genome']['fasta'],
        fai = lambda wildcards: config['genome']['fai_index']
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="vcf_filter", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start vcf_filter for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir,f"{wildcards.sample_id}/vcf_filter_{current_time}.sh")
            cmd1 = [
                params.gatk, "VariantFiltration",
                "--java-options", params.javaOptions,
                "-R", params.genome,
                "-V", input.vcf,
                "--window", "35",
                "--cluster", "3",
                "--filter-name", "FS",
                "--filter", "FS > 30.0",
                "--filter-name", "QD",
                "--filter", "QD < 2.0",
                "-O", output.vcf
            ]
            cmd2 = [
                params.bgzip, output.vcf
            ]
            if params.tmp_dir:
                cmd1.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during vcf_filter execution: {str(e)}\n")
            raise f"Error occurred while running vcf_filter for sample {wildcards.sample_id}: {e}, you can check the log file {log[0]} for more details."


