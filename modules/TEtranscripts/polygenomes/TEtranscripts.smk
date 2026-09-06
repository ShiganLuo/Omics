include: "../../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir= config.get("indir", "output/raw_fastq")
logdir_combine = config.get("logdir_combine", "log/combine")
ROOT_DIR = config.get("ROOT_DIR", "./")
genome_samples = config.get("genome_samples", [])

def get_input_for_TEcount(wildcards):
    logger.info(f"[get_input_for_TEcount] called with wildcards: {wildcards}")
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
        project = outdir + "/{genome}/TEcount/{sample_id}.TEcount.cntTable"
    params:
        project = "{sample_id}.TEcount",
        outdir = outdir + "/{genome}/TEcount",
        TEcount = config.get('Procedure',{}).get('TEcount') or 'TEcount'
    log:
        logdir + "/{sample_id}/{genome}/TEcount.log"
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("TEcount", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TEcount for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.project))
            script = os.path.join(sample_outdir, f"TEcount_{wildcards.sample_id}_{current_time}.sh")
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
                f.write(f"echo 'TEcount completed for sample {wildcards.sample_id} at {current_time}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during TEcount for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during TEcount for sample {wildcards.sample_id}: {e}")
            raise e

def get_input_for_combine_TEcount(wildcards):
    logger.info(f"[get_input_for_combine_TEcount] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in genome_samples.get(wildcards.genome, []):
        cntTable.append(f"{outdir}/{wildcards.genome}/TEcount/{sample_id}.TEcount.cntTable")
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TEcount didn't get any input files,samples:{genome_samples.get(wildcards.genome, [])}")
    return cntTable

rule combine_TEcount:
    input:
        fileList = get_input_for_combine_TEcount
    output:
        outfile = outdir + "/{genome}/TEcount/all_TEcount.tsv"
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml")
    params:
        combineTE = ROOT_DIR + "/modules/TEtranscripts/bin/combineTE.py",
        indir = outdir + "/{genome}/TEcount"
    log:
        logdir_combine + "/TEtranscripts/{genome}/combine_TEcount.log"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("combine_TEcount", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TEcount at {current_time}")
            sample_outdir = os.path.dirname(str(output.outfile))
            script = os.path.join(sample_outdir, f"combine_TEcount_{current_time}.sh")
            cmd = [
                "python", params.combineTE,
                "-p", "TEcount",
                "-i", params.indir,
                "-o", output.outfile
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'combine_TEcount completed at {current_time}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during combine_TEcount: {e}\n")
            logger.error(f"Error occurred during combine_TEcount: {e}")
            raise e

def get_input_for_TElocal(wildcards):
    logger.info(f"[get_input_for_TElocal] called with wildcards: {wildcards}")
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
        project = outdir + "/{genome}/TElocal/{sample_id}.TElocal.cntTable"
    log:
        logdir + "/{sample_id}/{genome}/TElocal.log"
    params:
        project = "{sample_id}.TElocal",
        TElocal = config.get('Procedure',{}).get('TElocal') or 'TElocal'
    threads: 2
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("TElocal", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TElocal for sample {wildcards.sample_id} at {current_time}")
            outdir_local = os.path.dirname(str(output.project))
            script = os.path.join(outdir_local, f"TElocal_{wildcards.sample_id}_{current_time}.sh")
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
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(f"echo 'TElocal completed for sample {wildcards.sample_id} at {current_time}'\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during TElocal for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during TElocal for sample {wildcards.sample_id}: {e}")
            raise e

def get_input_for_combine_TElocal(wildcards):
    logger.info(f"[get_input_for_combine_TElocal] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in genome_samples.get(wildcards.genome, []):
        cntTable.append(f"{outdir}/{wildcards.genome}/TElocal/{sample_id}.TElocal.cntTable")
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,samples:{genome_samples.get(wildcards.genome, [])}")
    return cntTable

rule combine_TElocal:
    input:
        fileList = get_input_for_combine_TElocal
    output:
        outfile = outdir + "/{genome}/TElocal/all_TElocal.tsv"
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml")
    params:
        combineTE = ROOT_DIR + "/modules/TEtranscripts/bin/combineTE.py",
        indir = outdir + "/{genome}/TElocal"
    log:
        logdir_combine + "/TEtranscripts/{genome}/combine_TElocal.log"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("combine_TElocal", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TElocal at {current_time}")
            script = os.path.join(params.indir, f"combine_TElocal_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during combine_TElocal: {e}\n")
            logger.error(f"Error occurred during combine_TElocal: {e}")
            raise e

rule TEtranscripts_result:
    input:
        TEcount = outdir + "/{genome}/TEcount/all_TEcount.tsv",
        TElocal = outdir + "/{genome}/TElocal/all_TElocal.tsv"
