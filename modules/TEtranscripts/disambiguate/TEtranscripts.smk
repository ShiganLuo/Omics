include: "../../common/common.smk"
from typing import Tuple
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "output/disambiguate")
ROOT_DIR = config.get("ROOT_DIR", "./")
genome_pairs: Tuple[str, str] = config.get("genome_pairs", ())
genomeA, genomeB = genome_pairs
samples = config.get("samples", [])

def get_input_for_TEcount(wildcards):
    logger.info(f"[get_input_for_TEcount] called with wildcards: {wildcards}")
    if wildcards.genome == genomeA:
        bam = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.disambiguatedSpecies_{wildcards.genome}.bam"
    elif wildcards.genome == genomeB:
        bam = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.disambiguatedSpecies_{wildcards.genome}.bam"
    else:
        raise ValueError(f"wildcards.genome {wildcards.genome} is not in genome_pairs {genome_pairs}")
    return bam

rule TEcount:
    input:
        bamA = get_input_for_TEcount
    output:
        project = outdir + "/TEcount/{genome}/{sample_id}.TEcount.cntTable"
    params:
        project = lambda wildcards: f"{wildcards.sample_id}.TEcount",
        outdir = lambda wildcards: outdir + f"/TEcount/{wildcards.genome}",
        TE_gtf = lambda wildcards: config['genome'][wildcards.genome]['TE_gtf'],
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf'],
        TEcount = config.get('Procedure',{}).get('TEcount') or 'TEcount'
    log:
        logdir + "/{sample_id}/{genome}/TEcount.log"
    threads: 2
    conda:
        "../TEtranscripts.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("TEcount", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TEcount for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            script = os.path.join(params.outdir, f"TEcount_{wildcards.sample_id}_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.TEcount} --sortByPos --format BAM --mode multi \\\n")
                f.write(f"    -b {input.bamA} --GTF {params.gtf} --TE {params.TE_gtf} \\\n")
                f.write(f"    --project {params.project} --outdir {params.outdir}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during TEcount for sample {wildcards.sample_id} genome {wildcards.genome}: {e}\n")
            logger.error(f"Error occurred during TEcount for sample {wildcards.sample_id} genome {wildcards.genome}: {e}")
            raise e

def get_cntTable_for_TEcount(wildcards):
    logger.info(f"[get_cntTable_for_TEcount] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in samples:
        cntTable.append(f"{outdir}/TEcount/{wildcards.genome}/{sample_id}.TEcount.cntTable")
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TEcount didn't get any input files,genome: {wildcards.genome}")
    return cntTable

rule combine_TEcount:
    input:
        get_cntTable_for_TEcount
    output:
        outfile_id = outdir + "/TEcount/{genome}/all_TEcount.tsv",
        outfile_name = outdir + "/TEcount/{genome}/all_TEcount_name.tsv"
    conda:
        "../TEtranscripts.yaml"
    params:
        combineTE = ROOT_DIR +"/modules/TEtranscripts/bin/combineTE.py",
        geneId2Name = ROOT_DIR +"/modules/TEtranscripts/bin/geneId2Name.py",
        indir = outdir + "/TEcount/{genome}",
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    threads: 1
    log:
        logdir + "/all/TEtranscripts/{genome}_combine_TEcount.log"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("combine_TEcount", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TEcount for genome {wildcards.genome} at {current_time}")
            script = os.path.join(params.indir, f"combine_TEcount_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile_id}\n")
                f.write(f"python {params.geneId2Name} -c {output.outfile_id} -g {params.gtf} -o {output.outfile_name}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during combine_TEcount for genome {wildcards.genome}: {e}\n")
            logger.error(f"Error occurred during combine_TEcount for genome {wildcards.genome}: {e}")
            raise e

def get_input_for_TElocal(wildcards):
    logger.info(f"[get_input_for_TElocal] called with wildcards: {wildcards}")
    if wildcards.genome == genomeA:
        bam = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.disambiguatedSpecies_{wildcards.genome}.bam"
    elif wildcards.genome == genomeB:
        bam = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.disambiguatedSpecies_{wildcards.genome}.bam"
    else:
        raise ValueError(f"wildcards.genome {wildcards.genome} is not in genome_pairs {genome_pairs}")
    return bam

rule TElocal:
    input:
        bam = get_input_for_TElocal
    output:
        project = outdir + "/TElocal/{genome}/{sample_id}.TElocal.cntTable"
    log:
        logdir + "/{sample_id}/{genome}/TElocal.log"
    params:
        project = lambda wildcards: f"{wildcards.sample_id}.TElocal",
        TE = lambda wildcards: config['genome'][wildcards.genome]['TEind'],
        GTF = lambda wildcards: config['genome'][wildcards.genome]['gtf'],
        TElocal = config.get('Procedure',{}).get('TElocal') or 'TElocal'
    threads: 2
    conda:
        "../TEtranscripts.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("TElocal", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TElocal for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
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
                f.write(f"Error occurred during TElocal for sample {wildcards.sample_id} genome {wildcards.genome}: {e}\n")
            logger.error(f"Error occurred during TElocal for sample {wildcards.sample_id} genome {wildcards.genome}: {e}")
            raise e

def get_cntTable_for_TElocal(wildcards):
    logger.info(f"[get_cntTable_for_TElocal] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in samples:
        cntTable.append(f"{outdir}/TElocal/{wildcards.genome}/{sample_id}.TElocal.cntTable")

    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,genome: {wildcards.genome}")
    return cntTable

rule combine_TElocal:
    input:
        get_cntTable_for_TElocal
    output:
        outfile_id = outdir + "/TElocal/{genome}/all_TElocal.tsv",
        outfile_name = outdir + "/TElocal/{genome}/all_TElocal_name.tsv"
    conda:
        "../TEtranscripts.yaml"
    params:
        combineTE = ROOT_DIR +"/modules/TEtranscripts/bin/combineTE.py",
        geneId2Name = ROOT_DIR +"/modules/TEtranscripts/bin/geneId2Name.py",
        indir = outdir + "/TElocal/{genome}",
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    threads: 1
    log:
        logdir + "/all/TEtranscripts/{genome}_combine_TElocal.log"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("combine_TElocal", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TElocal for genome {wildcards.genome} at {current_time}")
            script = os.path.join(params.indir, f"combine_TElocal_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile_id}\n")
                f.write(f"python {params.geneId2Name} -c {output.outfile_id} -g {params.gtf} -o {output.outfile_name}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during combine_TElocal for genome {wildcards.genome}: {e}\n")
            logger.error(f"Error occurred during combine_TElocal for genome {wildcards.genome}: {e}")
            raise e

rule TEtranscripts_result:
    input:
        TEcount = outdir + "/TEcount/{genome}/all_TEcount.tsv",
        TElocal = outdir + "/TElocal/{genome}/all_TElocal.tsv"
