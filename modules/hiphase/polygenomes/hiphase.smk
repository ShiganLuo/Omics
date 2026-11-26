include: "../../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
bam_dir = config.get("bam_dir", indir)
vcf_dir = config.get("vcf_dir", indir)
input_vcf_substring = config.get("input_vcf_substring") or ""
input_bam_substring = config.get("input_bam_substring") or ""
output_substring = config.get("output_substring") or ""
def get_input_for_hiphase_phase(wildcards):
    logger.info(f"hiphase_phase called with {wildcards}")
    in_dict = {}
    if input_bam_substring != "":
        in_dict["bam"] = os.path.join(bam_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.{input_bam_substring}.bam")
        in_dict["bai"] = os.path.join(bam_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.{input_bam_substring}.bai")
    else:
        in_dict["bam"] = os.path.join(bam_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(bam_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam.bai")
    if input_vcf_substring != "":
        in_dict["vcf"] = os.path.join(vcf_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.{input_vcf_substring}.vcf.gz")
        in_dict["csi"] = os.path.join(vcf_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.{input_vcf_substring}.vcf.gz.csi")
    else:
        in_dict["vcf"] = os.path.join(vcf_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.vcf.gz")
        in_dict["csi"] = os.path.join(vcf_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.vcf.gz.csi")
    in_dict["fasta"] = lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    logger.info(f"hiphase_phase input for sample {wildcards.sample_id} genome {wildcards.genome}: {in_dict}")
    return in_dict

rule hiphase_phase:
    input:
        unpack(get_input_for_hiphase_phase)
    output:
        vcf = outdir + "/{genome}/{sample_id}/{sample_id}.phased.vcf.gz" if output_substring == "" else outdir + "/{genome}/{sample_id}/{sample_id}." + output_substring + ".phased.vcf.gz",
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.phased.bam" if output_substring == "" else outdir + "/{genome}/{sample_id}/{sample_id}." + output_substring + ".phased.bam"
    log:
        logdir + "/{genome}/{sample_id}/hiphase.log"
    threads: rule_res("hiphase_phase").get("threads", 8)
    benchmark: logdir + "/{genome}/{sample_id}/hiphase.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "hiphase_phase", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 16000),
    params:
        hiphase = config.get("Procedure", {}).get("hiphase") or "hiphase",
    conda:
        "../hiphase.yaml"
    container:
        sif("../hiphase.yaml", "hiphase_phase")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start hiphase for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
        script = receipt_script(f"hiphase_{wildcards.sample_id}_{current_time}.sh")
        os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
        plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.sample_id}/hiphase.sh")
        cmd1 = [
            params.hiphase,
            "--threads", str(threads),
            "--reference", input.fasta,
            "--bam", input.bam,
            "--output-bam", output.bam,
            "--vcf", input.vcf,
            "--output-vcf", output.vcf
        ]
        with open(script, "w") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(f'echo "script completed successfully"\n')
        shutil.copyfile(script, plain_script)
        if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
            script_only_touch(output)
        else:
            shell("bash {plain_script} >> {log} 2>&1")
rule hiphase_result:
    input:
        vcf = outdir + "/{genome}/{sample_id}/{sample_id}.phased.vcf.gz" if output_substring == "" else outdir + "/{genome}/{sample_id}/{sample_id}." + output_substring + ".phased.vcf.gz",
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.phased.bam" if output_substring == "" else outdir + "/{genome}/{sample_id}/{sample_id}." + output_substring + ".phased.bam"
