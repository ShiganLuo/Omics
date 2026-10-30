
include: "../../common/common.smk"

import shutil
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
    container:
        sif("../gatk.yaml")
    params:
        gatk=gatk,
        java_options=java_options,
        tmp_dir=tmp_dir,
        interval=interval
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="gatk_population_haplotype_caller", log_file=log_path)
        try:
            os.makedirs(os.path.dirname(output.vcf), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start gatk_population_haplotype_caller for sample {wildcards.sample} at {current_time}")
            script = os.path.join(cmd_dir(logdir, str(wildcards.sample)), f"gatk_population_haplotype_caller_{wildcards.sample}_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(output.vcf), f"gatk_population_haplotype_caller_{wildcards.sample}.sh")
            cmd = [params.gatk, "--java-options", params.java_options, "HaplotypeCaller", "-R", input.fasta, "-I", input.bam, "-O", output.vcf, "-ERC", "GVCF", "--native-pair-hmm-threads", str(threads)]
            if params.interval:
                cmd += ["-L", params.interval]
            if params.tmp_dir:
                cmd += ["--tmp-dir", params.tmp_dir]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(str(item) for item in cmd) + "\n")
                f.write(f"{params.gatk} IndexFeatureFile -I {output.vcf}\n")
                f.write(f"echo 'gatk_population_haplotype_caller completed successfully for sample {wildcards.sample}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as exc:
            rule_logger.error(f"HaplotypeCaller failed for {wildcards.sample}: {exc}")
            raise RuntimeError(f"Error occurred while running HaplotypeCaller for sample {wildcards.sample}: {exc}, you can check the log file {log_path} for more details.")

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
    container:
        sif("../gatk.yaml")
    params:
        gatk=gatk,
        java_options=java_options,
        tmp_dir=tmp_dir,
        interval=interval,
        workspace=outdir + "/variants/joint/genomicsdb"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="gatk_population_joint_genotyping", log_file=log_path)
        try:
            os.makedirs(os.path.dirname(output.vcf), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_map = os.path.join(os.path.dirname(output.vcf), "sample_map.tsv")
            with open(sample_map, "w") as handle:
                for sample, gvcf in zip(samples, input.gvcfs):
                    handle.write(f"{sample}\t{gvcf}\n")
            script = os.path.join(cmd_dir(logdir, "group"), f"gatk_population_joint_genotyping_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(output.vcf), f"gatk_population_joint_genotyping.sh")
            import_cmd = [params.gatk, "--java-options", params.java_options, "GenomicsDBImport", "--sample-name-map", sample_map, "--genomicsdb-workspace-path", params.workspace, "--reader-threads", str(threads)]
            genotype_cmd = [params.gatk, "--java-options", params.java_options, "GenotypeGVCFs", "-R", input.fasta, "-V", "gendb://" + params.workspace, "-O", output.vcf]
            if params.interval:
                import_cmd += ["-L", params.interval]
                genotype_cmd += ["-L", params.interval]
            if params.tmp_dir:
                import_cmd += ["--tmp-dir", params.tmp_dir]
                genotype_cmd += ["--tmp-dir", params.tmp_dir]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write("rm -rf " + params.workspace + "\n")
                f.write(" ".join(str(item) for item in import_cmd) + "\n")
                f.write(" ".join(str(item) for item in genotype_cmd) + "\n")
                f.write(f"{params.gatk} IndexFeatureFile -I {output.vcf}\n")
                f.write(f"echo 'gatk_population_joint_genotyping completed successfully for genome {wildcards.genome}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as exc:
            rule_logger.error(f"Joint genotyping failed for {wildcards.sample}: {exc}")
            raise RuntimeError(f"Error occurred while running joint genotyping for sample {wildcards.sample}: {exc}, you can check the log file {log_path} for more details.")
