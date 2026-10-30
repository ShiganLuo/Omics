include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir= config.get("indir", "output/raw_fastq")
logdir_combine = config.get("logdir_combine", "log/combine")
ROOT_DIR = config.get("ROOT_DIR", "./")
genome_samples = config.get("genome_samples", {})

def get_input_for_TEcount(wildcards):
    logger.debug(f"[get_input_for_TEcount] called with wildcards: {wildcards}")
    bam = indir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    TE_gtf = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('TE_gtf')
    if not TE_gtf or not os.path.exists(TE_gtf):
        raise ValueError(f"TE_gtf file not found for genome {wildcards.genome}: {TE_gtf}")
    gtf = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"gtf file not found for genome {wildcards.genome}: {gtf}")
    in_dict = {
        "bam": bam,
        "TE_gtf": TE_gtf,
        "gtf": gtf
    }
    return in_dict

rule TEcount:
    input:
        unpack(get_input_for_TEcount)
    output:
        project = outdir + "/{genome}/{sample_id}/{sample_id}.TEcount.cntTable"
    threads: rule_res("TEcount").get("threads", 8)
    resources:
        mem_mb=rule_res("TEcount").get("mem_mb", 16000),
    params:
        project = "{sample_id}.TEcount",
        outdir = outdir + "/{genome}/{sample_id}",
        TEcount = config.get('Procedure',{}).get('TEcount') or 'TEcount'
    log:
        logdir + "/{sample_id}/{genome}/TEcount.log"
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("TEcount", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TEcount for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.project))
            script = os.path.join(RECEIPT_DIR, f"TEcount_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"TEcount_{wildcards.sample_id}.sh")
            cmd = [
                params.TEcount,
                "--sortByPos",
                "--format", "BAM",
                "--mode", "multi",
                "-b", input.bam,
                "--GTF", input.gtf,
                "--TE", input.TE_gtf,
                "--project", params.project,
                "--outdir", params.outdir
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'TEcount for sample {wildcards.sample_id} was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during TEcount for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"TEcount failed for sample {wildcards.sample_id}: {e}")

def get_input_for_combine_TEcount(wildcards):
    logger.debug(f"[get_input_for_combine_TEcount] called with wildcards: {wildcards}")
    cntTables = []
    for sample_id in genome_samples.get(wildcards.genome, []):
        cntTables.append(f"{outdir}/{wildcards.genome}/{sample_id}/{sample_id}.TEcount.cntTable")
    if len(cntTables) == 0:
        raise ValueError(f"rule combine_TEcount didn't get any input files,samples:{genome_samples.get(wildcards.genome, [])}")
    gtf = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"gtf file not found for genome {wildcards.genome}: {gtf}")
    in_dict = {
        "cntTable": cntTables,
        "gtf": gtf
    }
    return in_dict

rule combine_TEcount:
    input:
        unpack(get_input_for_combine_TEcount)
    output:
        outfile_raw = outdir + "/{genome}/{genome}_TEcount.tsv",
        outfile_annotate = outdir + "/{genome}/{genome}_TEcount_name.tsv",
    threads: rule_res("combine_TEcount").get("threads", 2)
    resources:
        mem_mb=rule_res("combine_TEcount").get("mem_mb", 8000),
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml")
    params:
        combineTE = ROOT_DIR + "/modules/TEtranscripts/bin/combineTE.py",
        geneId2Name = ROOT_DIR +"/modules/TEtranscripts/bin/geneId2Name.py",
        indir = outdir + "/{genome}"
    log:
        logdir_combine + "/TEtranscripts/{genome}/combine_TEcount.log"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("combine_TEcount", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TEcount at {current_time}")
            sample_outdir = os.path.dirname(str(output.outfile_raw))
            script = os.path.join(RECEIPT_DIR, f"combine_TEcount_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"combine_TEcount.sh")
            cmd1 = [
                "python", params.combineTE,
                "-p", "TEcount",
                "-i", params.indir,
                "-o", output.outfile_raw
            ]
            cmd2 = [
                "python", params.geneId2Name,
                "-c", output.outfile_raw,
                "-g", input.gtf,
                "-o", output.outfile_annotate
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(f"echo 'combine_TEcount for {wildcards.genome} completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during combine_TEcount: {e}")
            raise RuntimeError(f"Error occurred during combine_TEcount: {e}")

def get_input_for_TElocal(wildcards):
    logger.debug(f"[get_input_for_TElocal] called with wildcards: {wildcards}")
    bam = indir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    TEind = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('TEind')
    if not TEind or not os.path.exists(TEind):
        raise ValueError(f"TEind file not found for genome {wildcards.genome}: {TEind}")
    gtf = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"gtf file not found for genome {wildcards.genome}: {gtf}")
    in_dict = {
        "bam": bam,
        "TEind": TEind,
        "gtf": gtf
    }
    return in_dict

rule TElocal:
    input:
        unpack(get_input_for_TElocal)
    output:
        project = outdir + "/{genome}/{sample_id}/{sample_id}.TElocal.cntTable"
    log:
        logdir + "/{sample_id}/{genome}/TElocal.log"
    params:
        project = "{sample_id}.TElocal",
        TElocal = config.get('Procedure',{}).get('TElocal') or 'TElocal'
    threads: rule_res("TElocal").get("threads", 8)
    resources:
        mem_mb=rule_res("TElocal").get("mem_mb", 16000),
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("TElocal", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TElocal for sample {wildcards.sample_id} at {current_time}")
            outdir_local = os.path.dirname(str(output.project))
            script = os.path.join(RECEIPT_DIR, f"TElocal_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(outdir_local, f"TElocal_{wildcards.sample_id}.sh")
            cmd1 = [
                params.TElocal,
                "--sortByPos",
                "-b", input.bam,
                "--GTF", input.gtf,
                "--TE", input.TEind,
                "--project", params.project
            ]
            cmd2 = [
                "mv", f"{params.project}.cntTable", output.project
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(f"echo 'TElocal for sample {wildcards.sample_id} was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during TElocal for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during TElocal for sample {wildcards.sample_id}: {e}")

def get_input_for_combine_TElocal(wildcards):
    logger.debug(f"[get_input_for_combine_TElocal] called with wildcards: {wildcards}")
    cntTables = []
    for sample_id in genome_samples.get(wildcards.genome, []):
        cntTables.append(f"{outdir}/{wildcards.genome}/{sample_id}/{sample_id}.TElocal.cntTable")
    if len(cntTables) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,samples:{genome_samples.get(wildcards.genome, [])}")
    gtf = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"gtf file not found for genome {wildcards.genome}: {gtf}")
    in_dict = {
        "cntTables": cntTables,
        "gtf": gtf
    }
    return in_dict

rule combine_TElocal:
    input:
        unpack(get_input_for_combine_TElocal)
    output:
        outfile_raw = outdir + "/{genome}/{genome}_TElocal.tsv",
        outfile_annotate = outdir + "/{genome}/{genome}_TElocal_name.tsv",
    threads: rule_res("combine_TElocal").get("threads", 2)
    resources:
        mem_mb=rule_res("combine_TElocal").get("mem_mb", 8000),
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml")
    params:
        combineTE = ROOT_DIR + "/modules/TEtranscripts/bin/combineTE.py",
        indir = outdir + "/{genome}",
        geneId2Name = ROOT_DIR +"/modules/TEtranscripts/bin/geneId2Name.py"
    log:
        logdir_combine + "/TEtranscripts/{genome}/combine_TElocal.log"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("combine_TElocal", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TElocal at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"combine_TElocal_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(params.indir, f"combine_TElocal.sh")
            cmd1 = [
                "python", str(params.combineTE),
                "-p", "TElocal",
                "-i", str(params.indir),
                "-o", str(output.outfile_raw)
            ]
            cmd2 = [
                "python", params.geneId2Name,
                "-c", output.outfile_raw,
                "-g", input.gtf,
                "-o", output.outfile_annotate
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(f"echo 'combine_TElocal for {wildcards.genome} was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} > {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during combine_TElocal: {e}")
            raise RuntimeError(f"Error occurred during combine_TElocal: {e}")

rule TEtranscripts_result:
    input:
        TEcount = outdir + "/{genome}/{genome}_TEcount.tsv",
        TElocal = outdir + "/{genome}/{genome}_TElocal.tsv"
