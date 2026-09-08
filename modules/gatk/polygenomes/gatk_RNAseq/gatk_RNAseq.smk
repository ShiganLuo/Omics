include: "../../../common/common.smk"
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"

def get_input_for_SplitNCigarReads(wildcards):
    logger.info(f"[get_input_for_SplitNCigarReads] called with wildcards: {wildcards}")
    bam = indir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.sorted_markdup.bam"
    fasta = config.get('genome', {}).get('references',{}).get(wildcards.genome,{}).get('fasta')
    if not fasta or not os.path.exists(fasta):
        raise ValueError(f"Fasta file for genome is not specified or does not exist: {fasta}")
    dict_index = config.get('genome', {}).get('references',{}).get(wildcards.genome,{}).get('dict_index')
    fai_index = config.get('genome', {}).get('references',{}).get(wildcards.genome,{}).get('fai_index')
    if not dict_index or not os.path.exists(dict_index) or not fai_index or not os.path.exists(fai_index):
        fasta_link = indir + f"/index/{wildcards.genome}/{wildcards.genome}.fa"
        if not os.path.exists(fasta):
            os.symblink(fasta, fasta_link)
        fasta = fasta_link
        dict_index = indir + f"/index/{wildcards.genome}/{wildcards.genome}.dict"
        fai_index = indir + f"/index/{wildcards.genome}/{wildcards.genome}.fa.fai"        
    in_dict = {
        "bam": bam,
        "fasta": fasta,
        "dict_index": dict_index,
        "fai_index": fai_index,
    }
    return in_dict

rule SplitNCigarReads:
    input:
        unpack(get_input_for_SplitNCigarReads)
    output:
        bam = temp(outdir + "/{genome}/{sample_id}/{sample_id}.split.bam")
    params:
        javaOptions = config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        gatk = config.get("Procedure", {}).get("gatk") or "gatk"
    conda: "../../gatk.yaml"
    container:
        sif("../../gatk.yaml")
    threads: 8 
    log:
        logdir + "/{genome}/{sample_id}/SplitNCigarReads.log"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger(logger_name="SplitNCigarReads", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start SplitNCigarReads for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.dirname(str(output.bam))
            script = os.path.join(sample_outdir, f"SplitNCigarReads_{current_time}.sh")
            fasta = input.fasta

            cmd = [
                params.gatk, "SplitNCigarReads",
                "--java-options", params.javaOptions,
                "-R", fasta,
                "-I", input.bam,
                "-O", output.bam
            ]
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'SplitNCigarReads completed successfully for sample {wildcards.sample_id} genome {wildcards.genome}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during SplitNCigarReads execution: {str(e)}\n")
            logger.error(f"Error occurred during SplitNCigarReads for sample {wildcards.sample_id} genome {wildcards.genome}: {e}")
            raise e

def get_input_for_VarientCalling(wildcards):
    logger.info(f"[get_input_for_VarientCalling] called with wildcards: {wildcards}")
    bam = outdir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.split.bam"
    fasta = config.get('genome', {}).get('references',{}).get(wildcards.genome,{}).get('fasta')
    if not fasta or not os.path.exists(fasta):
        raise ValueError(f"Fasta file for genome is not specified or does not exist: {fasta}")
    fai_index = config.get('genome', {}).get('references',{}).get(wildcards.genome,{}).get('fai_index')
    if not fai_index or not os.path.exists(fai_index):
        fai_index = indir + f"/index/{wildcards.genome}/{wildcards.genome}.fa.fai"
        fasta_link = indir + f"/index/{wildcards.genome}/{wildcards.genome}.fa"
        if not os.path.exists(fasta):
            os.symblink(fasta, fasta_link)
        fasta = fasta_link
    in_dict = {
        "bam": bam,
        "fasta": fasta,
        "fai_index": fai_index
    }
    return in_dict

rule VarientCalling:
    input:
        unpack(get_input_for_VarientCalling)
    output:
        vcf = outdir + "/{genome}/{sample_id}/{sample_id}.raw.vcf.gz"
    log:
        logdir + "/{genome}/{sample_id}/VarientCalling.log"
    conda: "../../gatk.yaml"
    container:
        sif("../../gatk.yaml")
    params:
        javaOptions = config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
    threads: 8
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger(logger_name="VarientCalling", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start VarientCalling for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.dirname(str(output.vcf))
            script = os.path.join(sample_outdir, f"VarientCalling_{current_time}.sh")
            cmd = [
                params.gatk, "HaplotypeCaller",
                "--java-options", params.javaOptions,
                "-R", input.fasta,
                "-I", input.bam,
                "-O", output.vcf,
                "--dont-use-soft-clipped-bases",
                "--standard-min-confidence-threshold-for-calling", "20"
            ]
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'VarientCalling completed successfully for sample {wildcards.sample_id} genome {wildcards.genome}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during VarientCalling execution: {str(e)}\n")
            logger.error(f"Error occurred during VarientCalling for sample {wildcards.sample_id} genome {wildcards.genome}: {e}")
            raise e
def get_input_for_vcf_filter(wildcards):
    logger.info(f"[get_input_for_vcf_filter] called with wildcards: {wildcards}")
    vcf = outdir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.raw.vcf.gz"
    fasta = config.get('genome', {}).get('references',{}).get(wildcards.genome,{}).get('fasta')
    if not fasta or not os.path.exists(fasta):
        raise ValueError(f"Fasta file for genome is not specified or does not exist: {fasta}")
    fai_index = config.get('genome', {}).get('references',{}).get(wildcards.genome,{}).get('fai_index')
    if not fai_index or not os.path.exists(fai_index):
        fai_index = indir + f"/index/{wildcards.genome}/{wildcards.genome}.fa.fai"
        fasta_link = indir + f"/index/{wildcards.genome}/{wildcards.genome}.fa"
        if not os.path.exists(fasta):
            os.symblink(fasta, fasta_link)
        fasta = fasta_link
    in_dict = {
        "vcf": vcf,
        "fasta": fasta,
        "fai_index": fai_index
    }
    return in_dict

rule vcf_filter:
    input:
        unpack(get_input_for_vcf_filter)
    output:
        vcf = outdir + "/{genome}/{sample_id}/{sample_id}.filtered.vcf.gz"
    log:
        logdir + "/{genome}/{sample_id}/vcf_filter.log"
    conda: "../../gatk.yaml"
    container:
        sif("../../gatk.yaml")
    threads: 8
    params:
        javaOptions = config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        vcf = lambda wildcards: outdir + f"/{wildcards.genome}/SNP/vcf/filter/{wildcards.sample_id}.vcf",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        bgzip = config.get("Procedure", {}).get("bgzip") or "bgzip",
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger(logger_name="vcf_filter", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start vcf_filter for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.dirname(str(output.vcf))
            script = os.path.join(sample_outdir, f"vcf_filter_{current_time}.sh")
            cmd1 = [
                params.gatk, "VariantFiltration",
                "--java-options", params.javaOptions,
                "-R", input.fasta,
                "-V", input.vcf,
                "--window", "35",
                "--cluster", "3",
                "--filter-name", "FS",
                "--filter", "FS > 30.0",
                "--filter-name", "QD",
                "--filter", "QD < 2.0",
                "-O", output.vcf
            ]
            cmd2 = [
                params.bgzip, output.vcf
            ]
            if params.tmp_dir:
                cmd1.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:                
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(f"echo 'vcf_filter completed successfully for sample {wildcards.sample_id} genome {wildcards.genome}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during vcf_filter execution: {str(e)}\n")
            logger.error(f"Error occurred during vcf_filter for sample {wildcards.sample_id} genome {wildcards.genome}: {e}")
            raise e
