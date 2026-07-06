from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")
mode = config.get("mode") or None
def get_input_for_trimming_Paired(wildcards):
    logger.info(f"Getting input for trimming_Paired with mode: {mode}")
    if mode == "UMI":
        return [
            f"{indir}/{wildcards.sample_id}_1.umi.fq.gz",
            f"{indir}/{wildcards.sample_id}_2.umi.fq.gz"
        ]
    else:
        return [
            f"{indir}/{wildcards.sample_id}_1.fq.gz",
            f"{indir}/{wildcards.sample_id}_2.fq.gz"
        ]

rule trimming_Paired:
    input:
        get_input_for_trimming_Paired
    output:
        fastq1 = outdir + "/{sample_id}/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/{sample_id}/{sample_id}_2.fq.gz",
        report1 = outdir + "/{sample_id}/trimming_statistics_1.txt",
        report2 = outdir + "/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir = outdir + "/{sample_id}",
        quality = config.get('Params',{}).get("trim_galore", {}).get('quality') or 25,
        trim_galore = config.get('Procedure',{}).get('trim_galore') or 'trim_galore',
        adapters = config.get('Params',{}).get("trim_galore", {}).get('adapters') or None
    threads: 6
    log:
        log = logdir + "/{sample_id}/trimming.txt"
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"{wildcards.sample_id}/trimming_Paired_{current_time}.sh")

            cmd1 = [
                params.trim_galore, "--paired",
                "--cores", str(threads),
                "--quality", str(params.quality),
                "-o", params.outdir,
                "--basename", wildcards.sample_id,
                input[0], input[1]
            ]
            if params.adapters:
                cmd1 += ["--adapter", params.adapters]
            suffix1 = os.path.basename(input[0])
            suffix2 = os.path.basename(input[1])

            cmd2 = [
                mv, f"{params.outdir}/{wildcards.sample_id}_val_1.fq.gz", output.fastq1
            ]
            cmd3 = [
                mv, f"{params.outdir}/{wildcards.sample_id}_val_2.fq.gz", output.fastq2
            ]
            cmd4 = [
                mv, f"{params.outdir}/{suffix1}_trimming_report.txt", output.report1
            ]
            cmd5 = [
                mv, f"{params.outdir}/{suffix2}_trimming_report.txt", output.report2
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                f.write(" ".join(cmd4) + "\n")
                f.write(" ".join(cmd5) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"trimming_Paired failed for {wildcards.sample_id}: {e}\n")
            raise

def get_input_for_trimming_Single(wildcards):
    logger.info(f"Getting input for trimming_Single with mode: {mode}")
    if mode == "UMI":
        return f"{indir}/{wildcards.sample_id}.umi.single.fq.gz",
    else:
        return f"{indir}/{wildcards.sample_id}.single.fq.gz",

rule trimming_Single:
    input:
        fastq = get_input_for_trimming_Single
    output:
        fastq = outdir + "/{sample_id}/{sample_id}.single.fq.gz",
        report = outdir + "/{sample_id}/trimming_statistics.txt"
    params:
        outdir = outdir + "/{sample_id}",
        quality = config.get('Params',{}).get("trim_galore", {}).get('quality') or 25,
        trim_galore = config.get('Procedure',{}).get('trim_galore') or 'trim_galore',
        adapters = config.get('Params',{}).get("trim_galore", {}).get('adapters') or None
    threads: 6
    log:
        log = logdir + "/{sample_id}/trimming.txt"
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"{wildcards.sample_id}/trimming_Single_{current_time}.sh")
            cmd1 = [
                params.trim_galore, "--phred33",
                "--cores", str(threads),
                "--quality", str(params.quality),
                "-o", params.outdir,
                "--basename", wildcards.sample_id,
                input.fastq
            ]
            if params.adapters:
                cmd1 += ["--adapter", params.adapters]
            cmd2 = [
                mv, f"{params.outdir}/{wildcards.sample_id}_trimmed.fq.gz", output.fastq
            ]
            suffix = os.path.basename(input[0])
            cmd3 = [
                mv, f"{params.outdir}/{suffix}_trimming_report.txt", output.report
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"trimming_Single failed for {wildcards.sample_id}: {e}\n")
            raise




