
include: "../../../common/common.smk"
import os
import time
from snakemake.logging import logger

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
sample_bams = config.get("sample_bams", {})
gatk = config.get("Procedure", {}).get("gatk") or "gatk"
gatk_params = config.get("Params", {}).get("gatk", {})
java_options = gatk_params.get("javaOptions") or "-Xmx30g"
tmp_dir = gatk_params.get("tmp-dir")

if not samples or not sample_bams:
    raise ValueError("gatk_population requires samples and sample_bams")

def get_genome_fasta(wildcards):
    return config["genome"][wildcards.genome]["fasta"]

def get_genome_fai(wildcards):
    return config["genome"][wildcards.genome].get("fai_index") or config["genome"][wildcards.genome].get("fai")

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
    threads: 8
    conda:
        "../../gatk.yaml"
    container:
        sif("../../gatk.yaml")
    params:
        gatk=gatk,
        java_options=java_options,
        tmp_dir=tmp_dir,
        interval=get_genome_interval
    run:
        log_path = str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            os.makedirs(os.path.dirname(output.vcf), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script = os.path.join(os.path.dirname(output.vcf), f"gatk_population_haplotype_caller_{wildcards.sample}_{current_time}.sh")
            cmd = [params.gatk, "--java-options", params.java_options, "HaplotypeCaller", "-R", input.fasta, "-I", input.bam, "-O", output.vcf, "-ERC", "GVCF", "--native-pair-hmm-threads", str(threads)]
            if params.interval:
                cmd += ["-L", params.interval(wildcards)]
            if params.tmp_dir:
                cmd += ["--tmp-dir", params.tmp_dir]
            with open(script, "w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                handle.write(" ".join(str(item) for item in cmd) + "\n")
                handle.write(f"{params.gatk} IndexFeatureFile -I {output.vcf}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as exc:
            with open(log_path, "a") as handle:
                handle.write(f"HaplotypeCaller failed for {wildcards.sample} genome {wildcards.genome}: {exc}\n")
            raise

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
    threads: 8
    conda:
        "../../gatk.yaml"
    container:
        sif("../../gatk.yaml")
    params:
        gatk=gatk,
        java_options=java_options,
        tmp_dir=tmp_dir,
        interval=get_genome_interval,
        workspace=lambda wildcards: outdir + f"/{wildcards.genome}/variants/joint/genomicsdb"
    run:
        log_path = str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            os.makedirs(os.path.dirname(output.vcf), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_map = os.path.join(os.path.dirname(output.vcf), "sample_map.tsv")
            with open(sample_map, "w") as handle:
                for sample, gvcf in zip(samples, input.gvcfs):
                    handle.write(f"{sample}\t{gvcf}\n")
            script = os.path.join(os.path.dirname(output.vcf), f"gatk_population_joint_genotyping_{current_time}.sh")
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
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as exc:
            with open(log_path, "a") as handle:
                handle.write(f"Joint genotyping failed for genome {wildcards.genome}: {exc}\n")
            raise
