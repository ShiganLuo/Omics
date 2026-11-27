include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
# temp 规范(modules.md):本模块的中间文件用 maybe_temp(path, "<key>") 包裹,
# 是否标 temp 由 Params.<key>.keep_intermediates 决定(缺省 True=兜底保留)。


# first col: target(human) genome,second col: contaminating genome. human sample may contaminated by mouse genome

def get_input_for_XenofilteR(wildcards):
    logger.debug(f"[get_input_for_XenofilteR] called with wildcards: {wildcards}")
    host = config.get('Params', {}).get('XenofilteR', {}).get('sample_contamination', {}).get(wildcards.sample_id, {}).get('host')
    contaminant = config.get('Params', {}).get('XenofilteR', {}).get('sample_contamination', {}).get(wildcards.sample_id, {}).get('contaminant')
    host_bam = indir + f"/{host}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    if host == contaminant:
        contaminant_bam = host_bam
    else:
        contaminant_bam = indir + f"/{contaminant}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    in_dict = {
        "contaminant_bam": contaminant_bam,
        "host_bam": host_bam,
        # .bai 必须声明为 input:star_align 的 .bai 是 temp 输出,无人消费时
        # 该 job 一结束即被删除,skip 分支会 ln -s 出空链接,产物校验
        # MissingOutputException(实测 hTBLC3,3 次尝试全挂)。
        "host_bai": host_bam + ".bai",
        "contaminant_bai": contaminant_bam + ".bai",
    }
    logger.debug(f"[get_input_for_XenofilteR] returning input dictionary: {in_dict}")
    return in_dict

rule XenofilteR:
    input:
        unpack(get_input_for_XenofilteR)
    output:
        csvIn = maybe_temp(outdir + "/{genome}/{sample_id}/{sample_id}.csv", "XenofilteR"),
        outBam = maybe_temp(outdir + "/{genome}/{sample_id}/{sample_id}.bam", "XenofilteR"),
        outBai = maybe_temp(outdir + "/{genome}/{sample_id}/{sample_id}.bam.bai", "XenofilteR"),
    log:
        logdir + "/{sample_id}/{genome}/XenofilteR.log"
    threads: rule_res("XenofilteR").get("threads", 8)
    benchmark: logdir + "/{sample_id}/{genome}/XenofilteR.benchmark.tsv",
    resources:
        mem_mb=rule_res("XenofilteR").get("mem_mb", 16000),
    params:
        outdir = lambda wildcards: f"{outdir}/{wildcards.sample_id}",
        outSampleName = lambda wildcards: wildcards.sample_id,
        tempBam = lambda wildcards: f"{outdir}/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam",
        tempBai = lambda wildcards: f"{outdir}/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam.bai",
        MM = config.get('Parameters', {}).get('XenofilteR', {}).get('MM', 8),
        script = ROOT_DIR + "/modules/XenofilteR/bin/XenofilteR.r",
        Rscript = config.get('Procedure',{}).get('Rscript') or 'Rscript',
        host_genome = lambda wildcards: config.get("Params", {}).get("XenofilteR", {}).get("sample_contamination", {}).get(wildcards.sample_id, {}).get("host"),
        contaminant_genome = lambda wildcards: config.get("Params", {}).get("XenofilteR", {}).get("sample_contamination", {}).get(wildcards.sample_id, {}).get("contaminant"),
    conda:
        "XenofilteR.yaml"
    container:
        sif("XenofilteR.yaml", "XenofilteR")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("XenofilteR", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start XenofilteR for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.csvIn))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"XenofilteR_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"XenofilteR.sh")
            cmd1 = [
                "echo", f"{input.host_bam},{input.contaminant_bam}", ">", output.csvIn
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
                "touch", output.csvIn
            ]
            cmd6 = [
                "ln", "-s", input.host_bam, output.outBam
            ]
            # 复制而非软链:host_bai 是 temp,本 job 结束后即删,软链会让
            # 下游 addReadsGroup 拿到空链接;host_bam 非 temp,软链安全。
            cmd7 = [
                "cp", input.host_bai, output.outBai
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                if input.contaminant_bam == input.host_bam:
                    rule_logger.info(f"Host genome and pollution source genome are the same for sample {wildcards.sample_id}. Skipping XenofilteR filtering.")
                    f.write(" ".join(cmd5) + "\n")
                    f.write(" ".join(cmd6) + "\n")
                    f.write(" ".join(cmd7) + "\n")
                else:
                    rule_logger.info(f"Running XenofilteR for sample {wildcards.sample_id} with host genome {params.host_genome} and pollution source genome {params.contaminant_genome}.")
                    f.write(" ".join(cmd1) + "\n")
                    f.write(" ".join(cmd2) + "\n")
                    f.write(" ".join(cmd3) + "\n")
                    f.write(" ".join(cmd4) + "\n")
                f.write(f"echo 'XenofilteR completed successfully for sample {wildcards.sample_id}'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during XenofilteR for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during XenofilteR for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

