include: "../../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
_bam_suffix = config.get("bam_suffix", ".bam")

rule pbmm2_align:
    input:
        bam = indir + "/{sample_id}/{sample_id}" + _bam_suffix,
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.sorted.bam",
        bai = outdir + "/{genome}/{sample_id}/{sample_id}.sorted.bai"
    log:
        logdir + "/{sample_id}/{genome}/pbmm2_align.log"
    threads: rule_res("pbmm2_align").get("threads", 16)
    benchmark: logdir + "/{sample_id}/{genome}/pbmm2_align.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "pbmm2_align", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 32000),
    conda: "../pbmm2.yaml"
    container:
        sif("../pbmm2.yaml", "pbmm2_align")
    params:
        pbmm2 = config.get("Procedure", {}).get("pbmm2") or "pbmm2",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("pbmm2_align", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start pbmm2 alignment for sample {wildcards.sample_id} on genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.dirname(str(output.bam))
            os.makedirs(sample_outdir, exist_ok=True)
            # 带时间戳的备份(log/cmd 历史);outdir 只放无时间戳的当前命令(见"命令回执双写")
            script = receipt_script(f"pbmm2_align_{wildcards.genome}_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"pbmm2_align_{wildcards.sample_id}.sh")
            raw_align_bam = outdir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.raw.bam"
            cmd1 = [
                params.pbmm2, "align",
                "--sort",
                "--num-threads", str(threads),
                "--log-level", "INFO",
                input.fasta,
                input.bam,
                raw_align_bam
            ]
            cmd2 = [
                params.samtools, "sort",
                "-@", str(threads),
                raw_align_bam,
                "-o", output.bam
            ]
            cmd3 = [
                params.samtools, "index",
                "-@", str(threads),
                output.bam,
                "-o", output.bai
            ]
            cmd4 = [
                # pbmm2 --sort also writes the .bai companion — remove both
                "rm", "-f", raw_align_bam, raw_align_bam + ".bai",
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                f.write(" ".join(cmd4) + "\n")
                f.write(f'echo "pbmm2 alignment for sample {wildcards.sample_id} on genome {wildcards.genome} was completed successfully"\n')
            rule_logger.info("Executing: " + " && ".join(" ".join(str(x) for x in c) for c in (cmd1, cmd2, cmd3, cmd4)))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during pbmm2 alignment for sample {wildcards.sample_id} on genome {wildcards.genome}: {e}\n")
            raise RuntimeError(f"Error occurred during pbmm2 alignment for sample {wildcards.sample_id} on genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.\n")

rule pbmm2_result:
    input:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.sorted.bam",
        bai = outdir + "/{genome}/{sample_id}/{sample_id}.sorted.bai"
