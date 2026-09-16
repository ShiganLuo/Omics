include: "../../common/common.smk"

import re
from snakemake.logging import logger
indir = config.get("indir", "output")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", "./")
genome_paired_samples = config.get("genome_paired_samples", {})
genome_single_samples = config.get("genome_single_samples", {})

def get_bams_for_featureCounts_single(wildcards):
    logger.info(f"[get_bams_for_featureCounts_single] called with wildcards: {wildcards}")
    bams = []
    for sample_id in genome_single_samples.get(wildcards.genome, []):
        bams.append(f"{indir}/{wildcards.genome}/{sample_id}/{sample_id}.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_single_noMultiple: no single-end BAMs for genome {wildcards.genome}")
    gtf = config.get("genome",{}).get("references",{}).get(wildcards.genome,{}).get("gtf")
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"rule featureCounts_single_noMultiple: GTF file for genome {wildcards.genome} not found in config or does not exist: {gtf}")
    in_dict = {
        "bams": bams,
        "gtf": gtf
    }
    return in_dict

rule featureCounts_single_noMultiple:
    input:
        unpack(get_bams_for_featureCounts_single)
    output:
        outfile = outdir + "/{genome}/{genome}_single_featureCounts.tsv"
    log:
        logdir + "/{genome}/featureCounts_single_noMultiple.log"
    conda: 
        "../featureCounts.yaml"
    container:
        sif("../featureCounts.yaml")
    threads:
        10
    params:
        featureCounts = config.get('Procedure',{}).get('featureCounts') or 'featureCounts',
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("featureCounts_single_noMultiple", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start featureCounts_single_noMultiple at {current_time}")
            script = os.path.join(os.path.dirname(output.outfile), f"featureCounts_single_noMultiple_{current_time}.sh")
            os.makedirs(os.path.dirname(output.outfile), exist_ok=True)
            bams_str = " ".join(input.bams)
            cmd = [
                params.featureCounts,
                "-T", str(threads),
                "-t", "exon",
                "-g", "gene_id",
                "-a", str(input.gtf),
                "-o", str(output.outfile),
                bams_str
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'featureCounts_single_noMultiple for {wildcards.genome} completed at {current_time}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"featureCounts_single_noMultiple failed: {e}\n")
            logger.error(f"featureCounts_single_noMultiple failed: {e}\n")
            raise e

def get_bams_for_featureCounts_paired(wildcards):
    logger.info(f"[get_bams_for_featureCounts_paired] called with wildcards: {wildcards}")
    bams = []
    for sample_id in genome_paired_samples.get(wildcards.genome, []):
        bams.append(f"{indir}/{wildcards.genome}/{sample_id}/{sample_id}.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_paired_noMultiple: no paired-end BAMs for genome {wildcards.genome}")
    gtf =  config.get("genome", {}).get("references", {}).get(wildcards.genome, {}).get("gtf")
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"rule featureCounts_paired_noMultiple: GTF file for genome {wildcards.genome} not found in config or does not exist: {gtf}")
    in_dict = {
        "bams": bams,
        "gtf": gtf
    }
    return in_dict

rule featureCounts_paired_noMultiple:
    input:
        unpack(get_bams_for_featureCounts_paired)
    output:
        outfile = outdir + "/{genome}/{genome}_paired_featureCounts.tsv",
    log:
        logdir + "/{genome}/featureCounts_paired_noMultiple.log"
    conda:
        "../featureCounts.yaml"
    container:
        sif("../featureCounts.yaml")
    threads:
        10
    params:
        featureCounts = config.get('Procedure',{}).get('featureCounts') or 'featureCounts',
        M = config.get("Params",{}).get("featureCounts",{}).get("M") or False,
        O = config.get("Params",{}).get("featureCounts",{}).get("O") or False,
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("featureCounts_paired_noMultiple", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start featureCounts_paired_noMultiple at {current_time}")
            script = os.path.join(os.path.dirname(output.outfile), f"featureCounts_paired_noMultiple_{current_time}.sh")
            os.makedirs(os.path.dirname(output.outfile), exist_ok=True)
            bams_str = " ".join(input.bams)
            cmd = [
                params.featureCounts,
                "-T", str(threads),
                "-B", 
                "-p",
                "--countReadPairs",
                "-t", "exon",
                "-g", "gene_id",
                "-a", str(input.gtf),
                "-o", str(output.outfile)
            ]
            if params.M:
                cmd.append("-M")
            if params.O:
                cmd.append("-O")
            cmd.append(bams_str)
            with open(script, "w") as f:
                # for multiple -M -O
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'featureCounts_paired_noMultiple for {wildcards.genome} completed at {current_time}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"featureCounts_paired_noMultiple failed: {e}\n")
            logger.error(f"featureCounts_paired_noMultiple failed: {e}\n")
            raise e

def get_input_for_featureCounts_merge(wildcards):
    """Resolve merge inputs based on which sample layouts exist for this genome."""
    logger.info(f"[get_input_for_featureCounts_merge] called with wildcards: {wildcards}")
    inputs = {}
    if genome_paired_samples.get(wildcards.genome):
        inputs["paired"] = outdir + f"/{wildcards.genome}/{wildcards.genome}_paired_featureCounts.tsv"
    if genome_single_samples.get(wildcards.genome):
        inputs["single"] = outdir + f"/{wildcards.genome}/{wildcards.genome}_single_featureCounts.tsv"
    if not inputs:
        raise ValueError(f"No PE or SE samples for genome {wildcards.genome}, cannot build featureCounts input")
    return inputs

rule featureCounts_merge:
    """Merge PE and/or SE featureCounts output into a single count matrix.

    Input adapts dynamically via input function:
    - PE only → copy PE output as merged
    - SE only → copy SE output as merged
    - PE + SE → merge by Geneid column
    """
    input:
        unpack(get_input_for_featureCounts_merge)
    output:
        merged = outdir + "/{genome}/{genome}_featureCounts.tsv"
    log:
        logdir + "/{genome}/featureCounts_merge.log"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("featureCounts_merge", log_file=log_path)
            rule_logger.info(f"Start featureCounts merge for {wildcards.genome}")
            os.makedirs(os.path.dirname(output.merged), exist_ok=True)
            import pandas as pd
            has_paired = hasattr(input, 'paired') and input.paired
            has_single = hasattr(input, 'single') and input.single
            if has_paired and has_single:
                # Both PE + SE: merge by Geneid
                df_pe = pd.read_csv(input.paired, sep="\t", comment='#')
                df_pe.columns = [re.sub(r'.*\/([^\/]+)\.[^.]+$', r'\1', c) if '/' in c else c for c in df_pe.columns]
                df_se = pd.read_csv(input.single, sep="\t", comment='#')
                df_se.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'], inplace=True, errors='ignore')
                df_se.columns = [re.sub(r'.*\/([^\/]+)\.[^.]+$', r'\1', c) if '/' in c else c for c in df_se.columns]
                df_merged = pd.merge(df_pe, df_se, on='Geneid')
                df_merged.to_csv(output.merged, sep="\t", index=False)
                rule_logger.info(f"PE+SE merged: {df_merged.shape[0]} genes × {df_merged.shape[1] - 6} samples")
            elif has_paired:
                # PE only: copy
                import shutil
                shutil.copy2(input.paired, output.merged)
                rule_logger.info(f"PE-only: copied {input.paired} → {output.merged}")
            elif has_single:
                # SE only: copy
                import shutil
                shutil.copy2(input.single, output.merged)
                rule_logger.info(f"SE-only: copied {input.single} → {output.merged}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"featureCounts merge failed: {e}\n")
            raise e


rule featureCounts_results:
    input:
        merged = outdir + "/{genome}/{genome}_featureCounts.tsv"
