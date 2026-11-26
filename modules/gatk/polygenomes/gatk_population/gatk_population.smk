
include: "../../../common/common.smk"


outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
sample_bams = config.get("sample_bams", {})
gatk = config.get("Procedure", {}).get("gatk") or "gatk"
# 参数精确到 rule:javaOptions/tmp-dir 按 rule 独立配置(Params.gatk.<rule_name>
# .javaOptions/.tmp-dir),无 module 级共享。

if not samples or not sample_bams:
    raise ValueError("gatk_population requires samples and sample_bams")

def get_genome_fasta(wildcards):
    return config["genome"][wildcards.genome]["fasta"]

def get_genome_fai(wildcards):
    return config["genome"][wildcards.genome].get("fai_index")

def get_genome_dict(wildcards):
    return config["genome"][wildcards.genome].get("dict_index") or config["genome"][wildcards.genome].get("dict")

def get_genome_interval(wildcards):
    return config["genome"][wildcards.genome].get("interval")

rule gatk_population_haplotype_caller:
    input:
        bam=lambda wc: sample_bams[wc.sample],
        bai=lambda wc: sample_bams[wc.sample] + ".bai",
        fasta=get_genome_fasta,
        fai=get_genome_fai,
        dict=get_genome_dict
    output:
        vcf=outdir + "/{genome}/variants/gvcf/{sample}.g.vcf.gz",
        tbi=outdir + "/{genome}/variants/gvcf/{sample}.g.vcf.gz.tbi"
    log:
        logdir + "/{genome}/variants/{sample}/gatk_population_haplotype_caller.log"
    threads: rule_res("gatk_population_haplotype_caller").get("threads", 8)
    benchmark: logdir + "/{genome}/variants/{sample}/gatk_population_haplotype_caller.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "gatk_population_haplotype_caller", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../../gatk.yaml"
    container:
        sif("../../gatk.yaml", "gatk_population_haplotype_caller")
    params:
        gatk=gatk,
        java_options=config.get("Params", {}).get("gatk", {}).get("gatk_population_haplotype_caller", {}).get("javaOptions", "-Xmx30g"),
        tmp_dir=config.get("Params", {}).get("gatk", {}).get("gatk_population_haplotype_caller", {}).get("tmp-dir", None),
        interval=get_genome_interval
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger(logger_name="gatk_population_haplotype_caller", log_file=log_path)
        try:
            os.makedirs(os.path.dirname(output.vcf), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start gatk_population_haplotype_caller for sample {wildcards.sample} genome {wildcards.genome} at {current_time}")
            script = receipt_script(f"gatk_population_haplotype_caller_{wildcards.sample}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(os.path.dirname(output.vcf), f"gatk_population_haplotype_caller_{wildcards.sample}.sh")
            cmd = [params.gatk, "--java-options", params.java_options, "HaplotypeCaller", "-R", input.fasta, "-I", input.bam, "-O", output.vcf, "-ERC", "GVCF", "--native-pair-hmm-threads", str(threads)]
            if params.interval:
                cmd += ["-L", params.interval(wildcards)]
            if params.tmp_dir:
                cmd += ["--tmp-dir", params.tmp_dir]
            with open(script, "w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                handle.write(" ".join(str(item) for item in cmd) + "\n")
                handle.write(f"{params.gatk} IndexFeatureFile -I {output.vcf}\n")
                handle.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as exc:
            rule_logger.error(f"HaplotypeCaller failed for {wildcards.sample} genome {wildcards.genome}: {exc}")
            raise RuntimeError(f"HaplotypeCaller failed for {wildcards.sample} genome {wildcards.genome}: {exc}, you can check the log file {log_path} for more details.")

rule gatk_population_joint_genotyping:
    input:
        gvcfs=expand(outdir + "/{{genome}}/variants/gvcf/{sample}.g.vcf.gz", sample=samples),
        indexes=expand(outdir + "/{{genome}}/variants/gvcf/{sample}.g.vcf.gz.tbi", sample=samples),
        fasta=get_genome_fasta,
        fai=get_genome_fai,
        dict=get_genome_dict
    output:
        vcf=outdir + "/{genome}/variants/joint/joint.vcf.gz",
        tbi=outdir + "/{genome}/variants/joint/joint.vcf.gz.tbi"
    log:
        logdir + "/{genome}/variants/gatk_population_joint_genotyping.log"
    threads: rule_res("gatk_population_joint_genotyping").get("threads", 8)
    benchmark: logdir + "/{genome}/variants/gatk_population_joint_genotyping.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "gatk_population_joint_genotyping", input_size=input_size_named(input, 'vcf'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../../gatk.yaml"
    container:
        sif("../../gatk.yaml", "gatk_population_joint_genotyping")
    params:
        gatk=gatk,
        java_options=config.get("Params", {}).get("gatk", {}).get("gatk_population_joint_genotyping", {}).get("javaOptions", "-Xmx30g"),
        tmp_dir=config.get("Params", {}).get("gatk", {}).get("gatk_population_joint_genotyping", {}).get("tmp-dir", None),
        interval=get_genome_interval,
        workspace=lambda wildcards: outdir + f"/{wildcards.genome}/variants/joint/genomicsdb"
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger(logger_name="gatk_population_joint_genotyping", log_file=log_path)
        try:
            os.makedirs(os.path.dirname(output.vcf), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start gatk_population_joint_genotyping for genome {wildcards.genome} at {current_time}")
            sample_map = os.path.join(os.path.dirname(output.vcf), "sample_map.tsv")
            with open(sample_map, "w") as handle:
                for sample, gvcf in zip(samples, input.gvcfs):
                    handle.write(f"{sample}\t{gvcf}\n")
            script = receipt_script(f"gatk_population_joint_genotyping_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(os.path.dirname(output.vcf), f"gatk_population_joint_genotyping.sh")
            import_cmd = [params.gatk, "--java-options", params.java_options, "GenomicsDBImport", "--sample-name-map", sample_map, "--genomicsdb-workspace-path", params.workspace, "--reader-threads", str(threads)]
            genotype_cmd = [params.gatk, "--java-options", params.java_options, "GenotypeGVCFs", "-R", input.fasta, "-V", "gendb://" + params.workspace, "-O", output.vcf]
            interval_val = params.interval(wildcards)
            if interval_val:
                import_cmd += ["-L", interval_val]
                genotype_cmd += ["-L", interval_val]
            if params.tmp_dir:
                import_cmd += ["--tmp-dir", params.tmp_dir]
                genotype_cmd += ["--tmp-dir", params.tmp_dir]
            with open(script, "w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                handle.write("rm -rf " + params.workspace + "\n")
                handle.write(" ".join(str(item) for item in import_cmd) + "\n")
                handle.write(" ".join(str(item) for item in genotype_cmd) + "\n")
                handle.write(f"{params.gatk} IndexFeatureFile -I {output.vcf}\n")
                handle.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as exc:
            rule_logger.error(f"Joint genotyping failed for genome {wildcards.genome}: {exc}")
            raise RuntimeError(f"Joint genotyping failed for genome {wildcards.genome}: {exc}, you can check the log file {log_path} for more details.")
