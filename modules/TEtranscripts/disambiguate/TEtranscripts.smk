include: "../../common/common.smk"

from typing import Tuple
# temp 规范(modules.md):本模块的中间文件用 maybe_temp(path, "TEtranscripts") 包裹,
# 是否标 temp 由 Params.TEtranscripts.keep_intermediates 决定(缺省 True=兜底保留)。
# per-sample cntTable 是 combine 的中间输入;all_TEcount/all_TElocal 合并表是交付物。
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
        project = maybe_temp(outdir + "/TEcount/{genome}/{sample_id}.TEcount.cntTable", "TEtranscripts")
    params:
        project = lambda wildcards: f"{wildcards.sample_id}.TEcount",
        outdir = lambda wildcards: outdir + f"/TEcount/{wildcards.genome}",
        TE_gtf = lambda wildcards: config['genome'][wildcards.genome]['TE_gtf'],
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf'],
        TEcount = config.get('Procedure',{}).get('TEcount') or 'TEcount'
    log:
        logdir + "/{sample_id}/{genome}/TEcount.log"
    threads: rule_res("TEcount").get("threads", 2)
    benchmark: logdir + "/{sample_id}/{genome}/TEcount.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "TEcount", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml", "TEcount")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("TEcount", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TEcount for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"TEcount_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(params.outdir, f"TEcount_{wildcards.sample_id}.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.TEcount} --sortByPos --format BAM --mode multi \\\n")
                f.write(f"    -b {input.bamA} --GTF {params.gtf} --TE {params.TE_gtf} \\\n")
                f.write(f"    --project {params.project} --outdir {params.outdir}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during TEcount for sample {wildcards.sample_id} genome {wildcards.genome}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

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
    container:
        sif("../TEtranscripts.yaml", "combine_TEcount")
    params:
        combineTE = ROOT_DIR +"/modules/TEtranscripts/bin/combineTE.py",
        geneId2Name = ROOT_DIR +"/modules/TEtranscripts/bin/geneId2Name.py",
        indir = outdir + "/TEcount/{genome}",
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    threads: rule_res("combine_TEcount").get("threads", 1)
    benchmark: logdir + "/all/TEtranscripts/{genome}_combine_TEcount.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "combine_TEcount", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    log:
        logdir + "/all/TEtranscripts/{genome}_combine_TEcount.log"
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("combine_TEcount", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TEcount for genome {wildcards.genome} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"combine_TEcount_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(params.indir, f"combine_TEcount.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile_id}\n")
                f.write(f"python {params.geneId2Name} -c {output.outfile_id} -g {params.gtf} -o {output.outfile_name}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during combine_TEcount for genome {wildcards.genome}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

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
        project = maybe_temp(outdir + "/TElocal/{genome}/{sample_id}.TElocal.cntTable", "TEtranscripts")
    log:
        logdir + "/{sample_id}/{genome}/TElocal.log"
    params:
        project = lambda wildcards: f"{wildcards.sample_id}.TElocal",
        TE = lambda wildcards: config['genome'][wildcards.genome]['TEind'],
        GTF = lambda wildcards: config['genome'][wildcards.genome]['gtf'],
        TElocal = config.get('Procedure',{}).get('TElocal') or 'TElocal'
    threads: rule_res("TElocal").get("threads", 2)
    benchmark: logdir + "/{sample_id}/{genome}/TElocal.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "TElocal", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../TEtranscripts.yaml"
    container:
        sif("../TEtranscripts.yaml", "TElocal")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("TElocal", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TElocal for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            outdir_local = os.path.dirname(str(output.project))
            script = os.path.join(RECEIPT_DIR, f"TElocal_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir_local, f"TElocal_{wildcards.sample_id}.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.TElocal} --sortByPos -b {input.bam} \\\n")
                f.write(f"    --GTF {params.GTF} --TE {params.TE} \\\n")
                f.write(f"    --project {params.project}\n")
                f.write(f"mv {params.project}.cntTable {output.project}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during TElocal for sample {wildcards.sample_id} genome {wildcards.genome}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

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
    container:
        sif("../TEtranscripts.yaml", "combine_TElocal")
    params:
        combineTE = ROOT_DIR +"/modules/TEtranscripts/bin/combineTE.py",
        geneId2Name = ROOT_DIR +"/modules/TEtranscripts/bin/geneId2Name.py",
        indir = outdir + "/TElocal/{genome}",
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    threads: rule_res("combine_TElocal").get("threads", 1)
    benchmark: logdir + "/all/TEtranscripts/{genome}_combine_TElocal.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "combine_TElocal", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    log:
        logdir + "/all/TEtranscripts/{genome}_combine_TElocal.log"
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("combine_TElocal", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start combine_TElocal for genome {wildcards.genome} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"combine_TElocal_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(params.indir, f"combine_TElocal.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile_id}\n")
                f.write(f"python {params.geneId2Name} -c {output.outfile_id} -g {params.gtf} -o {output.outfile_name}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during combine_TElocal for genome {wildcards.genome}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule TEtranscripts_result:
    input:
        TEcount = outdir + "/TEcount/{genome}/all_TEcount.tsv",
        TElocal = outdir + "/TElocal/{genome}/all_TElocal.tsv"
