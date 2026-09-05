include: "../../common/common.smk"
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
logdir_combine = config.get("logdir_combine") or "log"
genome_samples = config.get("genome_samples") or {}
sample_groups = config.get("sample_groups") or {}
ROOT_DIR = config.get("ROOT_DIR", ".")
rule stringTie:
    input:
        bam = indir + "/{genome}/{sample_id}/{sample_id}.bam"
    output:
        gtf = outdir + "/{genome}/raw/{sample_id}/{sample_id}.gtf"
    log:
        logdir + "/{genome}/{sample_id}/stringTie.log"
    params:
        gtf = lambda wildcards: config.get('genomes', {}).get('reference', {}).get(wildcards.genome, {}).get('gtf'),
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    threads: 5
    conda:
        "../StringTie.yaml"
    container:
        sif("../StringTie.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="stringTie_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start stringTie run for sample {wildcards.sample_id} at {current_time}")
            script = f"{outdir}/{wildcards.genome}/raw/{wildcards.sample_id}/stringTie_{current_time}.sh"
            cmd = [params.stringtie, "-o", output.gtf, input.bam, "-G", params.gtf, "-p", str(threads)]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error occurred during stringTie run: {e}\n")
            raise f"Error occurred during stringTie run: {e}"


rule TEChimericTranscripts:
    input:
        gtf = outdir + "/{genome}/raw/{sample_id}/{sample_id}.gtf"
    output:
        txt = outdir + "/{genome}/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt"
    log:
        logdir + "/{genome}/{sample_id}/TEChimericTranscripts.log"
    params:
        te_gtf = lambda wildcards: config.get('genomes', {}).get('reference', {}).get(wildcards.genome, {}).get('TE_gtf'),
        TEChimericTranscripts = ROOT_DIR + "/modules/StringTie/bin/TEChimericTranscripts.py"
    threads: 5
    conda:
        "../StringTie.yaml"
    container:
        sif("../StringTie.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger(logger_name="TEChimericTranscripts_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            rule_logger.info(f"Start TEChimericTranscripts run for sample {wildcards.sample_id} at {current_time}")
            script = f"{outdir}/{wildcards.genome}/raw/{wildcards.sample_id}/TEChimericTranscripts.{current_time}.sh"
            cmd = f"python {params.TEChimericTranscripts} -s {input.gtf} -t {params.te_gtf} -o {output.txt} > {log} 2>&1"
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(cmd + "\n")
                f.write(f'echo "TEChimericTranscripts for sample {wildcards.sample_id} on genome {wildcards.genome} successfully completed!"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error occurred during TEChimericTranscripts: {e}\n")
            raise RuntimeError(f"Error occurred during TEChimericTranscripts: {e}\n")

def get_input_for_TEChimericPlot(wildcards):
    logger.info(f"[get_input_for_TEChimericPlot] called with wildcards: {wildcards}")
    if wildcards.genome not in genome_samples:
        raise ValueError(f"Genome {wildcards.genome} not found in genome_samples configuration.")
    txts = []
    for genome, samples in genome_samples.items():
        if genome == wildcards.genome:
            txts += [outdir + f"/{genome}/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt" for sample_id in samples]
    if len(txts) == 0:
        raise ValueError(f"No TE chimeric transcript files found for genome {wildcards.genome}.")
    return txts
rule TEChimericPlot:
    input:
        txts = get_input_for_TEChimericPlot
    output:
        group_stack = outdir + "/{genome}/TE_chimeric/TE_chimeric_group_stacked.png",
        type_top = outdir + "/{genome}/TE_chimeric/TE_chimeric_te_type_top.png",
        type_by_group = outdir + "/{genome}/TE_chimeric/TE_chimeric_te_type_by_group.png",
        sample_summary = outdir + "/{genome}/TE_chimeric/TE_chimeric_sample_summary.tsv",
        group_summary = outdir + "/{genome}/TE_chimeric/TE_chimeric_group_summary.tsv",
        te_type_counts = outdir + "/{genome}/TE_chimeric/TE_chimeric_te_type_counts.tsv"
    log:
        logdir_combine + "/{genome}/stringtie/TEChimericPlot.log"
    params:
        TEChimericPlot = ROOT_DIR + "/modules/StringTie/bin/TEChimericPlot.py"
    threads: 1
    conda:
        "../StringTie.yaml"
    container:
        sif("../StringTie.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger(logger_name="TEChimericPlot_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = f"{outdir}/{wildcards.genome}/TE_chimeric/TEChimericPlot.{current_time}.sh"
            rule_logger.info(f"Start TEChimericPlot run at {current_time}")
            group_tsv = outdir + f"/{wildcards.genome}/TE_chimeric/sample_groups.tsv"
            with open(group_tsv, 'w') as f:
                f.write("sample\tgroup\n")
                for group, sample_list in sample_groups.items():
                    for sample_id in sample_list:
                        f.write(f"{sample_id}\t{group}\n")
            cmd = [
                "python", params.TEChimericPlot,
                "-i", os.path.join(outdir, wildcards.genome, "raw"),
                "-g", group_tsv,
                "-o", f"{outdir}/{wildcards.genome}/TE_chimeric/TE_chimeric"
            ]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(' '.join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error occurred during TEChimericPlot: {e}\n")
            raise RuntimeError(f"Error occurred during TEChimericPlot: {e}\n")


def get_input_for_stringTieMerge(wildcards):
    logger.info(f"[get_input_for_stringTieMerge] called with wildcards: {wildcards}")
    gtfs = []
    for sample_id in genome_samples.get(wildcards.genome, []):
        gtfs.append(outdir + f"/{wildcards.genome}/raw/{sample_id}/{sample_id}.gtf")
    if len(gtfs) == 0:
        raise ValueError(f"No GTF files found for genome {wildcards.genome}.")
    return gtfs

rule stringTieMerge:
    input:
        gtfs = get_input_for_stringTieMerge
    output:
        gtf = outdir + "/{genome}/stringtie_merged.gtf"
    log:
        logdir_combine + "/{genome}/stringtie/stringTieMerge.log"
    params:
        gtf = lambda wildcards: config.get('genomes', {}).get('reference', {}).get(wildcards.genome, {}).get('gtf'),
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    conda:
        "../StringTie.yaml"
    container:
        sif("../StringTie.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger(logger_name="stringTieMerge_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start stringTieMerge run at {current_time}")
            script = os.path.join(outdir, wildcards.genome, f"stringTieMerge_{current_time}.sh")
            cmd = [params.stringtie, "--merge"] + list(input.gtfs) + ["-o", output.gtf, "-G", params.gtf]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error occurred during stringTieMerge: {e}\n")
            raise RuntimeError(f"Error occurred during stringTieMerge: {e}\n")
