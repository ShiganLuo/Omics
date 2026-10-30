include: "../common/common.smk"

import shutil

from snakemake.logging import logger
indir = config.get("indir", "output")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", "./")
single_samples = config.get("single_samples", [])
paired_samples = config.get("paired_samples", [])

# need test
def get_bams_for_featureCounts_single(wildcards):
    logger.debug(f"[get_bams_for_featureCounts_single] called with wildcards: {wildcards}")
    bams = []
    suffix = config.get("bam_substring", "")
    suffix_str = f".{suffix}" if suffix else ""
    for sample_id in single_samples:
        bams.append(f"{indir}/{sample_id}/{sample_id}{suffix_str}.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_single_noMultiple didn't get any input bams")
    return bams

rule featureCounts_single_noMultiple:
    input:
        bams = get_bams_for_featureCounts_single
    output:
        outfile = outdir + "/all_single_featureCounts.tsv"
    log:
        logdir + "/all/featureCounts/featureCounts_single_noMultiple.log"
    conda: 
        "featureCounts.yaml"
    container:
        sif("featureCounts.yaml")
    threads:
        10
    params:
        featureCounts = config.get('Procedure',{}).get('featureCounts') or 'featureCounts',
        gtf = config["genome"]["gtf"]
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("featureCounts_single_noMultiple", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start featureCounts_single_noMultiple at {current_time}")
            os.makedirs(os.path.join(logdir, "cmd", "group"), exist_ok=True)
            script = os.path.join(logdir, "cmd", "group", f"featureCounts_single_noMultiple_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(output.outfile), f"featureCounts_single_noMultiple.sh")
            os.makedirs(os.path.dirname(output.outfile), exist_ok=True)
            bams_str = " ".join(input.bams)
            with open(script, "w") as f:
                f.write(f"{params.featureCounts} -T {threads} -t exon -g gene_id -a {params.gtf} -o {output.outfile} {bams_str} > {log_path} 2>&1\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"featureCounts_single_noMultiple completed")
        except Exception as e:
            rule_logger.error(f"featureCounts_single_noMultiple failed: {e}\n")
            raise RuntimeError(f"featureCounts_single_noMultiple failed: {e}\n")

def get_bams_for_featureCounts_paired(wildcards):
    logger.debug(f"[get_bams_for_featureCounts_paired] called with wildcards: {wildcards}")
    bams = []
    suffix = config.get("bam_substring", "")
    suffix_str = f".{suffix}" if suffix else ""
    for sample_id in paired_samples:
        bams.append(f"{indir}/{sample_id}/{sample_id}{suffix_str}.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_paired_noMultiple didn't get any input bams")
    return bams

rule featureCounts_paired_noMultiple:
    input:
        bams = get_bams_for_featureCounts_paired
    output:
        outfile = outdir + "/all_paired_featureCounts.tsv",
    log:
        logdir + "/all/featureCounts/featureCounts_paired_noMultiple.log"
    conda:
        "featureCounts.yaml"
    container:
        sif("featureCounts.yaml")
    threads:
        10
    params:
        featureCounts = config.get('Procedure',{}).get('featureCounts') or 'featureCounts',
        gtf = config["genome"]["gtf"]
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("featureCounts_paired_noMultiple", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start featureCounts_paired_noMultiple at {current_time}")
            os.makedirs(os.path.join(logdir, "cmd", "group"), exist_ok=True)
            script = os.path.join(logdir, "cmd", "group", f"featureCounts_paired_noMultiple_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(output.outfile), f"featureCounts_paired_noMultiple.sh")
            os.makedirs(os.path.dirname(output.outfile), exist_ok=True)
            bams_str = " ".join(input.bams)
            with open(script, "w") as f:
                # for multiple -M -O
                f.write("/bin/bash\nset -euo pipefail\n")
                f.write(f"{params.featureCounts} -T {threads} -B -p --countReadPairs -t exon -g gene_id -a {params.gtf} -o {output.outfile} {bams_str}\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"featureCounts_paired_noMultiple completed")
        except Exception as e:
            rule_logger.error(f"featureCounts_paired_noMultiple failed: {e}\n")
            raise RuntimeError(f"featureCounts_paired_noMultiple failed: {e}\n")


rule featureCounts_result:
    input:
        paired = outdir + "/all_paired_featureCounts.tsv",
        single = outdir + "/all_single_featureCounts.tsv"
    output:
        touch(outdir + "/featureCounts.done")
