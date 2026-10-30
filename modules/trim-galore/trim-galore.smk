include: "../common/common.smk"

import shutil
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")
mode = config.get("mode") or None
def get_input_for_trimming_Paired(wildcards):
    logger.debug(f"Getting input for trimming_Paired with mode: {mode}")
    if mode == "UMI":
        return [
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.umi.fq.gz",
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.umi.fq.gz"
        ]
    else:
        return [
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz",
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
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
    threads: rule_res("trimming_Paired").get("threads", 6)
    resources:
        mem_mb=rule_res("trimming_Paired").get("mem_mb", 8000),
    log:
        log = logdir + "/{sample_id}/trimming.txt"
    conda:
        "trim-galore.yaml"
    container:
        sif("trim-galore.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("trimming_Paired", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start trimming_Paired for sample {wildcards.sample_id} at {current_time}")
            os.makedirs(os.path.join(logdir, "cmd", str(wildcards.sample_id)), exist_ok=True)
            script_path = os.path.join(logdir, "cmd", str(wildcards.sample_id), f"{wildcards.sample_id}/trimming_Paired_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/trimming_Paired.sh")
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
                "mv", f"{params.outdir}/{wildcards.sample_id}_val_1.fq.gz", output.fastq1
            ]
            cmd3 = [
                "mv", f"{params.outdir}/{wildcards.sample_id}_val_2.fq.gz", output.fastq2
            ]
            cmd4 = [
                "mv", f"{params.outdir}/{suffix1}_trimming_report.txt", output.report1
            ]
            cmd5 = [
                "mv", f"{params.outdir}/{suffix2}_trimming_report.txt", output.report2
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                f.write(" ".join(cmd4) + "\n")
                f.write(" ".join(cmd5) + "\n")
                f.write(f'echo "trimming_Paired for sample {wildcards.sample_id} successfully completed !"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"trimming_Paired failed for {wildcards.sample_id}: {e}")
            raise RuntimeError(f"trimming_Paired failed for {wildcards.sample_id}: {e}\n")

def get_input_for_trimming_Single(wildcards):
    logger.debug(f"Getting input for trimming_Single with mode: {mode}")
    if mode == "UMI":
        return f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.umi.single.fq.gz"
    else:
        return f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"

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
    threads: rule_res("trimming_Single").get("threads", 6)
    resources:
        mem_mb=rule_res("trimming_Single").get("mem_mb", 4000),
    log:
        log = logdir + "/{sample_id}/trimming.txt"
    conda:
        "trim-galore.yaml"
    container:
        sif("trim-galore.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("trimming_Single", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            os.makedirs(os.path.join(logdir, "cmd", str(wildcards.sample_id)), exist_ok=True)
            script_path = os.path.join(logdir, "cmd", str(wildcards.sample_id), f"{wildcards.sample_id}/trimming_Single_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/trimming_Single.sh")
            cmd1 = [
                params.trim_galore, "--phred33",
                "--cores", str(threads),
                "--quality", str(params.quality),
                "-o", params.outdir,
                "--basename", wildcards.sample_id,
                str(input.fastq)
            ]
            if params.adapters:
                cmd1 += ["--adapter", params.adapters]
            cmd2 = [
                "mv", f"{params.outdir}/{wildcards.sample_id}_trimmed.fq.gz", output.fastq
            ]
            suffix = os.path.basename(input[0])
            cmd3 = [
                "mv", f"{params.outdir}/{suffix}_trimming_report.txt", output.report
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                f.write(f'echo "trimming_Single for sample {wildcards.sample_id} successfully completed !"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"trimming_Single failed for {wildcards.sample_id}: {e}")
            raise RuntimeError(f"trimming_Single failed for {wildcards.sample_id}: {e}\n")




