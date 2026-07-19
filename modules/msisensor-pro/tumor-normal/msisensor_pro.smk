include: "../../common/common.smk"

# ======================================================================================================================
# Project: Project_ABC
# Script : build_baseline.smk
# Author : Peng Jia
# Date   :  2024/5/27
# Email  : pengjia@xjtu.edu.cn
# Description: Pipeline of baseline building of MSIsensor-pro
# ======================================================================================================================
import pandas as pd


outdir = config["outdir"]
outdir = outdir if outdir.endswith("/") else outdir + "/"
samples = pd.read_table(f"{config['sample_path']}",index_col=0)

msisensor_pro = config["msisensor-pro"]
genome_version = config["genome_version"]

genome_reference = {config["genome_version"]: config["reference"]}


rule all:
    input:
        outdir + f"tumor_normal_output.{genome_version}.merge.tsv"


rule scan:
    input:
        lambda wildcards: genome_reference[wildcards.genome_version]
    output:
        outdir + "reference/{genome_version}.msisensor.scan.list"
    log:
        outdir + "reference/{genome_version}.msisensor.scan.log"
    conda:
        "../msisensor_pro.schema.yaml"
    run:
        shell("{msisensor_pro} scan -d {input} -o {output} 2>{log} 1>{log}")

rule msisensor_msi:
    input:
        t=lambda wildcards: samples.loc[wildcards.case, "tumor_path"],
        n=lambda wildcards: samples.loc[wildcards.case, "normal_path"],
        ms=outdir + "reference/{genome_version}.msisensor.scan.list",
        ref=lambda wildcards: genome_reference[wildcards.genome_version]
    output:
        outdir + "tumor_normal_output/{case}/{case}.{genome_version}.msisensor"
    log:
        outdir + "tumor_normal_output/{case}/{case}.{genome_version}.msisensor.log"
    conda:
        "../msisensor_pro.schema.yaml"
    run:
        shell("{msisensor_pro} msi -d {input.ms} -n {input.n} -t {input.t} -g {input.ref}  -o {output} 2>{log}")

rule merge_msi_result:
    input:
        expand(outdir + "tumor_normal_output/{case}/{case}.{{genome_version}}.msisensor",case=samples.index)
    output:
        outdir + "tumor_normal_output.{genome_version}.merge.tsv"
    run:
        output_info = pd.DataFrame(columns=["Total_number_of_sites", "Number_of_unstable_sites", "MSI_score"])
        output_info.index.name="case_name"
        for case in samples.index:
            value_info = [i.split() for i in open(outdir + f"tumor_normal_output/{case}/{case}.{wildcards.genome_version}.msisensor")]
            output_info.loc[case] = value_info[1]
        output_info.to_csv(f"{output}",sep="\t")
# print(value_info)