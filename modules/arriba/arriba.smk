include: "../common/common.smk"
import os
import time
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
samples = config.get("samples",[])
bam_substring = config.get("bam_substring") or ""
def get_input_for_arriba(wildcards):
    logger.info("called rule arriba by {wildcards}")
    in_dict = {}
    if bam_substring != "" :
        in_dict["bam"] = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.{bam_substring}.bam"
    else:
        in_dict["bam"] = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    in_dict["fasta"] = config.get('genome',{}).get('fasta')
    in_dict["gtf"] = config.get('genome',{}).get('gtf')
    in_dict["blacklist"] = config.get('Params',{}).get('arriba',{}).get('blacklist')
    in_dict["known_fusions"] = config.get('Params',{}).get('arriba',{}).get('known_fusions')
    return in_dict
rule arriba:
    input:
        unpack(get_input_for_arriba)
    output:
        passed_fusion_tsv = outdir + "/{sample_id}/{sample_id}_passed_fusions.tsv",
        discarded_fusion_tsv = outdir + "/{sample_id}/{sample_id}_discarded_fusions.tsv",
    log:
        logdir + "/{sample_id}/arriba.log"
    threads: 4
    params:
        arriba = config.get('Procedure',{}).get('arriba') or 'arriba',
        t = config.get('Params',{}).get('arriba',{}).get('t') or None,
        d = config.get('Params',{}).get('arriba',{}).get('d') or None,
        E = config.get('Params',{}).get('arriba',{}).get('E') or 0.3,
        p = config.get('Params',{}).get('arriba',{}).get('p') or None
    conda:
        "arriba.yaml"
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/{wildcards.sample_id}/{wildcards.sample_id}_arriba.{current_time}.sh"
        if not os.path.exists(input.blacklist):
            raise ValueError(f"Blacklist file not found: {input.blacklist}")
        if not os.path.exists(input.known_fusions):
            raise ValueError(f"Known fusions file not found: {input.known_fusions}")
        cmd = [
            params.arriba,
            "-x", input.bam,
            "-o", output.passed_fusion_tsv,
            "-O", output.discarded_fusion_tsv,
            "-a", input.fasta,
            "-g", input.gtf,
            "-b", input.blacklist,
            "-k", input.known_fusions,
            "-E", str(params.E)
        ]
        if params.t:
            cmd += ["-t", params.t]
        if params.d:
            cmd += ["-d", params.d]
        if params.p:
            cmd += ["-p", params.p]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule arriba_report:
    input:
        passed_fusions = expand(outdir + "/{sample_id}/{sample_id}_passed_fusions.tsv", sample_id=samples),
        discarded_fusions = expand(outdir + "/{sample_id}/{sample_id}_discarded_fusions.tsv", sample_id=samples)
    output:
        report = outdir + "/../arriba_report/arriba_fusion_report.html"
    log:
        logdir + "/all/arriba_report.log"
    params:
        summary_script = os.path.join(ROOT_DIR, "modules/arriba/bin/summarize_arriba_fusions.py")
    conda:
        "arriba.yaml"
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/arriba_report.{current_time}.sh"
        cmd = [
            "python", params.summary_script,
            "-p", ",".join(input.passed_fusions),
            "-d", ",".join(input.discarded_fusions),
            "-o", outdir + "/../arriba_report"
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")
rule arriba_result:
    input:
        passed_fusion_tsv = outdir + "/{sample_id}/{sample_id}_passed_fusions.tsv",
        discarded_fusion_tsv = outdir + "/{sample_id}/{sample_id}_discarded_fusions.tsv"
