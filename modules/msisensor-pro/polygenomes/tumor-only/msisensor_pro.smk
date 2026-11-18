include: "../../../common/common.smk"

# ======================================================================================================================
# Project: Project_ABC
# Script : build_baseline.smk
# Author : Peng Jia
# Date   :  2024/5/27
# Email  : pengjia@xjtu.edu.cn
# Description: Pipeline of baseline building of MSIsensor-pro
# ======================================================================================================================
import pandas as pd

print(config)
outdir = config["outdir"]
outdir = outdir if outdir.endswith("/") else outdir + "/"
baseline_bams = config["samples4baseline"]
baseline_samples = pd.read_table(f"{baseline_bams}",index_col=0)
samples = pd.read_table(f"{config['samples4test']}",index_col=0)
msisensor_pro = config["msisensor-pro"]
genome_version = config["genome_version"]
print(samples)
print(baseline_samples)
genome_reference = {config["genome_version"]: config["reference"]}
# print(outdir + f"baselines/{genome_version}.baseline.file")

rule all:
    input:
        outdir+ f"tumor_only_output.{genome_version}.merge.tsv"

        # outdir + f"baselines/{genome_version}.baseline.file"


rule scan:
    input:
        lambda wildcards: genome_reference[wildcards.genome_version]
    output:
        outdir + "reference/{genome_version}.msisensor.scan.list"
    log:
        outdir + "reference/{genome_version}.msisensor.scan.log"
    benchmark: outdir + "reference/{genome_version}.msisensor.scan.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "scan", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../../msisensor-pro.yaml"
    container:
        sif("../../msisensor-pro.yaml", "scan")
    run:
        if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
            script_only_touch(output)
        else:
            shell("{msisensor_pro} scan -d {input} -o {output} 2>{log} 1>{log}")
rule pro_pre:
    input:
        bam=lambda wildcards: baseline_samples.loc[wildcards.sample, "bam_path"],
        ms=outdir + "reference/{genome_version}.msisensor.scan.list",
        ref=lambda wildcards: genome_reference[wildcards.genome_version],
    output:
        outdir + "baselines/details/{sample}.{genome_version}.baseline.out"
    log:
        outdir + "baselines/details/{sample}.{genome_version}.baseline.out.log"
    threads: rule_res("pro_pre").get("threads", 2)
    benchmark: outdir + "baselines/details/{sample}.{genome_version}.baseline.out.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "pro_pre", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../../msisensor-pro.yaml"
    container:
        sif("../../msisensor-pro.yaml", "pro_pre")
    run:
        if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
            script_only_touch(output)
        else:
            shell("{msisensor_pro} pro -d {input.ms} -t {input.bam} -g {input.ref} -o {output} >> {log} 2>&1")
rule merge:
    input:
        expand(outdir + "baselines/details/{sample}.{{genome_version}}.baseline.out",sample=baseline_samples.index)
    output:
        outdir + "baselines/{genome_version}.baseline.samples.list"
    log:
        outdir + "baselines/{genome_version}.baseline.samples.list.log"
    benchmark: outdir + "baselines/{genome_version}.baseline.samples.list.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "merge", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../../msisensor-pro.yaml"
    container:
        sif("../../msisensor-pro.yaml", "merge")
    run:
        with open(f"{output}","w") as f:
            for i in input:
                f.write(f"{i.split('/')[-1]} \t{i}_all\n")

rule baseline:
    input:
        conf=outdir + "baselines/{genome_version}.baseline.samples.list",
        ms=outdir + "reference/{genome_version}.msisensor.scan.list",
    output:
        outdir + "baselines/{genome_version}.baseline.tsv"
    log:
        outdir + "baselines/{genome_version}.baseline.tsv.log"
    benchmark: outdir + "baselines/{genome_version}.baseline.tsv.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "baseline", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../../msisensor-pro.yaml"
    container:
        sif("../../msisensor-pro.yaml", "baseline")
    run:
        if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
            script_only_touch(output)
        else:
            shell(f"{msisensor_pro} baseline  -d {input.ms} -i {input.conf} -o {output} -s 1 >> {log} 2>&1")
rule run_pro:
    input:
        ms = outdir + "baselines/{genome_version}.baseline.tsv",
        t=lambda wildcards: samples.loc[wildcards.case, "tumor_path"],
        ref=lambda wildcards: genome_reference[wildcards.genome_version]
    output:
        outdir + "tumor_only_output/{case}/{case}.{genome_version}.msisensor-pro"
    log:
        outdir + "tumor_only_output/{case}/{case}.{genome_version}.msisensor-pro.log"
    benchmark: outdir + "tumor_only_output/{case}/{case}.{genome_version}.msisensor-pro.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "run_pro", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../../msisensor-pro.yaml"
    container:
        sif("../../msisensor-pro.yaml", "run_pro")
    run:
        if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
            script_only_touch(output)
        else:
            shell("{msisensor_pro} pro -d {input.ms} -t {input.t} -g {input.ref}  -o {output} 2>{log}")
rule merge_msi_result:
    input:
        expand(outdir + "tumor_only_output/{case}/{case}.{{genome_version}}.msisensor-pro",case=samples.index)
    output:
        outdir + "tumor_only_output.{genome_version}.merge.tsv"
    log:
        outdir + "tumor_only_output.{genome_version}.merge.tsv.log"
    benchmark: outdir + "tumor_only_output.{genome_version}.merge.tsv.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "merge_msi_result", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    run:
        output_info = pd.DataFrame(columns=["Total_number_of_sites", "Number_of_unstable_sites", "MSI_score"])
        output_info.index.name="case_name"
        for case in samples.index:
            value_info = [i.split() for i in open(outdir + f"tumor_only_output/{case}/{case}.{wildcards.genome_version}.msisensor-pro")]
            output_info.loc[case] = value_info[1]
        output_info.to_csv(f"{output}",sep="\t")
