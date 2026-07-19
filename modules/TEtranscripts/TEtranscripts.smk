include: "../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir= config.get("indir", "output/raw_fastq")
ROOT_DIR = config.get("ROOT_DIR", "./")
samples = config.get("samples", [])


rule TEcount:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        project = outdir + "/TEcount/{sample_id}.TEcount.cntTable"
    params:
        project = "{sample_id}.TEcount",
        outdir = outdir + "/TEcount",
        TE_gtf = lambda wildcards: config['genome']['TE_gtf'],
        gtf = lambda wildcards: config['genome']['gtf'],
        TEcount = config.get('Procedure',{}).get('TEcount') or 'TEcount'
    log:
        logdir + "/{sample_id}/TEcount.log"
    conda:
        "TEtranscripts.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("TEcount", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TEcount for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(params.outdir, f"TEcount_{wildcards.sample_id}_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.TEcount} --sortByPos --format BAM --mode multi \\\n")
                f.write(f"    -b {input.bam} --GTF {params.gtf} --TE {params.TE_gtf} \\\n")
                f.write(f"    --project {params.project} --outdir {params.outdir}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during TEcount for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during TEcount for sample {wildcards.sample_id}: {e}")
            raise e

def get_cntTable_for_TEcount(wildcards):
    logger.info(f"[get_cntTable_for_TEcount] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in samples:
        cntTable.append(f"{outdir}/TEcount/{sample_id}.TEcount.cntTable")

    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TEcount didn't get any input files,samples:{samples}")
    return cntTable

rule combine_TEcount:
    input:
        fileList = get_cntTable_for_TEcount
    output:
        outfile = outdir + "/TEcount/all_TEcount.tsv"
    conda:
        "TEtranscripts.yaml"
    params:
        combineTE = ROOT_DIR + "/modules/TEtranscripts/bin/combineTE.py",
        indir = outdir + "/TEcount"
    log:
        logdir + "/all/TEtranscripts/combine_TEcount.log"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("combine_TEcount", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TEcount at {current_time}")
            script = os.path.join(params.indir, f"combine_TEcount_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during combine_TEcount: {e}\n")
            logger.error(f"Error occurred during combine_TEcount: {e}")
            raise e

rule TElocal:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        project = outdir + "/TElocal/{sample_id}.TElocal.cntTable"
    log:
        logdir + "/{sample_id}/TElocal.log"
    params:
        project = "{sample_id}.TElocal",
        TE = lambda wildcards: config['genome']['TEind'],
        GTF = lambda wildcards: config['genome']['gtf'],
        TElocal = config.get('Procedure',{}).get('TElocal') or 'TElocal'
    threads: 2
    conda:
        "TEtranscripts.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("TElocal", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TElocal for sample {wildcards.sample_id} at {current_time}")
            outdir_local = os.path.dirname(str(output.project))
            script = os.path.join(outdir_local, f"TElocal_{wildcards.sample_id}_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.TElocal} --sortByPos -b {input.bam} \\\n")
                f.write(f"    --GTF {params.GTF} --TE {params.TE} \\\n")
                f.write(f"    --project {params.project}\n")
                f.write(f"mv {params.project}.cntTable {output.project}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during TElocal for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during TElocal for sample {wildcards.sample_id}: {e}")
            raise e

def get_cntTable_for_TElocal(wildcards):
    logger.info(f"[get_cntTable_for_TElocal] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in samples:
        cntTable.append(f"{outdir}/TElocal/{sample_id}.TElocal.cntTable")

    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,samples:{samples}")
    return cntTable

rule combine_TElocal:
    input:
        fileList = get_cntTable_for_TElocal
    output:
        outfile = outdir + "/TElocal/all_TElocal.tsv"
    conda:
        "TEtranscripts.yaml"
    params:
        combineTE = ROOT_DIR + "/modules/TEtranscripts/bin/combineTE.py",
        indir = outdir + "/TElocal"
    log:
        logdir + "/all/TEtranscripts/combine_TElocal.log"
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
        TEcount = outdir + "/TEcount/all_TEcount.tsv",
        TElocal = outdir + "/TElocal/all_TElocal.tsv"
