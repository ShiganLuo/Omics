include: "../common/common.smk"

outdir = config.get("outdir", "output")
indir = config.get("indir", "input")
logdir = config.get("logdir", "logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
mode = config.get("mode") or None


def get_cutadapt_input_paired(wildcards):
    """Return paired-end FASTQ inputs, with optional UMI suffix."""
    if mode == "UMI":
        return [
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.umi.fq.gz",
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.umi.fq.gz",
        ]
    return [
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz",
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz",
    ]


def get_cutadapt_input_single(wildcards):
    """Return single-end FASTQ input, with optional UMI suffix."""
    if mode == "UMI":
        return f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.umi.single.fq.gz"
    return f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"


rule trimming_Paired:
    input:
        get_cutadapt_input_paired
    output:
        fastq1 = outdir + "/{sample_id}/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/{sample_id}/{sample_id}_2.fq.gz",
        report  = outdir + "/{sample_id}/{sample_id}.cutadapt_report.txt",
    log:
        logdir + "/{sample_id}/cutadapt_paired.log"
    threads: 4
    conda:
        "cutadapt.yaml"
    params:
        cutadapt     = config.get("Procedure", {}).get("cutadapt") or "cutadapt",
        quality      = config.get("Params", {}).get("cutadapt", {}).get("quality") or 20,
        adapter_r1   = config.get("Params", {}).get("cutadapt", {}).get("adapter_r1") or None,
        adapter_r2   = config.get("Params", {}).get("cutadapt", {}).get("adapter_r2") or None,
        minimum_length = config.get("Params", {}).get("cutadapt", {}).get("minimum_length") or None,
        maximum_length = config.get("Params", {}).get("cutadapt", {}).get("maximum_length") or None,
        match_read_wildcards = config.get("Params", {}).get("cutadapt", {}).get("match_read_wildcards") or False,
        cut = config.get("Params", {}).get("cutadapt", {}).get("cut") or None,
        trimmed_only = config.get("Params", {}).get("cutadapt", {}).get("trimmed_only") or False
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="cutadapt_paired", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start cutadapt paired trimming for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/cutadapt_paired_{current_time}.sh")
            cmd = [
                params.cutadapt,
                "-j", str(threads),
                "-q", str(params.quality),
                "--report=minimal",
                "-o", output.fastq1,
                "-p", output.fastq2,
            ]
            if params.adapter_r1:
                for adapter in params.adapter_r1.split(" "):
                    cmd += ["-a", adapter]
            if params.adapter_r2:
                for adapter in params.adapter_r2.split(" "):
                    cmd += ["-A", adapter]
            if params.minimum_length:
                cmd += ["-m", str(params.minimum_length)]
            if params.maximum_length:
                cmd += ["-M", str(params.maximum_length)]
            if params.match_read_wildcards:
                cmd += ["--match-read-wildcards"]
            if params.cut:
                cmd += ["--cut", str(params.cut)]
            if params.trimmed_only:
                cmd += ["--trimmed-only"]
            cmd += [input[0], input[1]]
            cmd_str = " ".join(cmd)
            with open(script, "w") as f:
                f.write(cmd_str + "\n")
                f.write(f"echo 'cutadapt paired trimming finished for {wildcards.sample_id}'\n")
            shell(f"bash {script} > {output.report} 2>> {log_path}")
            logger.info(f"Finished cutadapt paired trimming for sample {wildcards.sample_id}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"cutadapt paired trimming failed for sample {wildcards.sample_id} with error: {e}\n")
            raise


rule trimming_Single:
    input:
        get_cutadapt_input_single
    output:
        fastq  = outdir + "/{sample_id}/{sample_id}.single.fq.gz",
        report = outdir + "/{sample_id}/{sample_id}.single.cutadapt_report.txt",
    log:
        logdir + "/{sample_id}/cutadapt_single.log"
    threads: 4
    conda:
        "cutadapt.yaml"
    params:
        cutadapt       = config.get("Procedure", {}).get("cutadapt") or "cutadapt",
        quality        = config.get("Params", {}).get("cutadapt", {}).get("quality") or 20,
        adapter_r1     = config.get("Params", {}).get("cutadapt", {}).get("adapter_r1") or None,
        minimum_length = config.get("Params", {}).get("cutadapt", {}).get("minimum_length") or None,
        maximum_length = config.get("Params", {}).get("cutadapt", {}).get("maximum_length") or None,
        match_read_wildcards = config.get("Params", {}).get("cutadapt", {}).get("match_read_wildcards") or False,
        cut = config.get("Params", {}).get("cutadapt", {}).get("cut") or None,
        trimmed_only = config.get("Params", {}).get("cutadapt", {}).get("trimmed_only") or False,
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="cutadapt_single", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start cutadapt single trimming for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/cutadapt_single_{current_time}.sh")
            cmd = [
                params.cutadapt,
                "-j", str(threads),
                "-q", str(params.quality),
                "--report=minimal",
                "-o", output.fastq,
            ]
            if params.adapter_r1:
                for adapter in params.adapter_r1.split(" "):
                    cmd += ["-a", adapter]
            if params.minimum_length:
                cmd += ["-m", str(params.minimum_length)]
            if params.maximum_length:
                cmd += ["-M", str(params.maximum_length)]
            if params.match_read_wildcards:
                cmd += ["--match-read-wildcards"]
            if params.cut:
                cmd += ["--cut", str(params.cut)]
            if params.trimmed_only:
                cmd += ["--trimmed-only"]
            cmd += [input]
            cmd_str = " ".join(cmd)
            with open(script, "w") as f:
                f.write(cmd_str + "\n")
                f.write(f"echo 'cutadapt single trimming finished for {wildcards.sample_id}'\n")
            shell(f"bash {script} > {output.report} 2>> {log_path}")
            logger.info(f"Finished cutadapt single trimming for sample {wildcards.sample_id}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"cutadapt single trimming failed for sample {wildcards.sample_id} with error: {e}\n")
            raise
