
include: "../../common/common.smk"
import os
import time
from snakemake.logging import logger

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
sample_bams = config.get("sample_bams", {})
fasta = config.get("genome", {}).get("fasta")
fai = config.get("genome", {}).get("fai_index") or config.get("genome", {}).get("fai")
dict_index = config.get("genome", {}).get("dict_index") or config.get("genome", {}).get("dict")
interval = config.get("genome", {}).get("interval")
gatk = config.get("Procedure", {}).get("gatk") or "gatk"
gatk_params = config.get("Params", {}).get("gatk", {})
java_options = gatk_params.get("javaOptions") or "-Xmx30g"
tmp_dir = gatk_params.get("tmp-dir")

if not samples or not sample_bams or not fasta or not fai or not dict_index:
    raise ValueError("gatk_population requires samples, sample_bams, genome.fasta, genome.fai_index, and genome.dict_index")

rule gatk_population_haplotype_caller:
    input:
        bam=lambda wc: sample_bams[wc.sample],
        bai=lambda wc: sample_bams[wc.sample] + ".bai",
        fasta=fasta,
        fai=fai,
        dict=dict_index
    output:
        vcf=outdir + "/variants/gvcf/{sample}.g.vcf.gz",
        tbi=outdir + "/variants/gvcf/{sample}.g.vcf.gz.tbi"
    log:
        logdir + "/variants/{sample}/gatk_population_haplotype_caller.log"
    threads: 8
    conda:
        "../gatk.yaml"
    params:
        gatk=gatk,
        java_options=java_options,
        tmp_dir=tmp_dir,
        interval=interval
    run:
        log_path = str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            os.makedirs(os.path.dirname(output.vcf), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script = os.path.join(os.path.dirname(output.vcf), f"gatk_population_haplotype_caller_{wildcards.sample}_{current_time}.sh")
            cmd = [params.gatk, "--java-options", params.java_options, "HaplotypeCaller", "-R", input.fasta, "-I", input.bam, "-O", output.vcf, "-ERC", "GVCF", "--native-pair-hmm-threads", str(threads)]
            if params.interval:
                cmd += ["-L", params.interval]
            if params.tmp_dir:
                cmd += ["--tmp-dir", params.tmp_dir]
            with open(script, "w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                handle.write(" ".join(str(item) for item in cmd) + "\n")
                handle.write(f"{params.gatk} IndexFeatureFile -I {output.vcf}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as exc:
            with open(log_path, "a") as handle:
                handle.write(f"HaplotypeCaller failed for {wildcards.sample}: {exc}\n")
            raise

rule gatk_population_joint_genotyping:
    input:
        gvcfs=expand(outdir + "/variants/gvcf/{sample}.g.vcf.gz", sample=samples),
        indexes=expand(outdir + "/variants/gvcf/{sample}.g.vcf.gz.tbi", sample=samples),
        fasta=fasta,
        fai=fai,
        dict=dict_index
    output:
        vcf=outdir + "/variants/joint/joint.vcf.gz",
        tbi=outdir + "/variants/joint/joint.vcf.gz.tbi"
    log:
        logdir + "/variants/gatk_population_joint_genotyping.log"
    threads: 8
    conda:
        "../gatk.yaml"
    params:
        gatk=gatk,
        java_options=java_options,
        tmp_dir=tmp_dir,
        interval=interval,
        workspace=outdir + "/variants/joint/genomicsdb"
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
            if params.interval:
                import_cmd += ["-L", params.interval]
                genotype_cmd += ["-L", params.interval]
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
                handle.write(f"Joint genotyping failed: {exc}\n")
            raise
