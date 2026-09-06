include: "../../common/common.smk"
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
logdir_combine = config.get("logdir_combine") or "log"
genome_samples = config.get("genome_samples") or {}
sample_groups = config.get("sample_groups") or {}
ROOT_DIR = config.get("ROOT_DIR", ".")

def get_input_for_stringTie(wildcards):
    logger.info(f"[get_input_for_stringTie] called with wildcards: {wildcards}")
    bam_path = indir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    gtf = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"GTF file for genome {wildcards.genome} is not specified or does not exist in the configuration.")
    in_dict = {
        "bam": bam_path,
        "gtf": gtf
    }
    return in_dict

rule stringTie:
    input:
        unpack(get_input_for_stringTie)
    output:
        gtf = outdir + "/{genome}/raw/{sample_id}/{sample_id}.gtf"
    log:
        logdir + "/{sample_id}/{genome}/stringTie.log"
    params:
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
            sample_outdir = os.path.dirname(output.gtf)
            script = f"{sample_outdir}/stringTie_{current_time}.sh"
            cmd = [params.stringtie, "-o", output.gtf, input.bam, "-G", input.gtf, "-p", str(threads)]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "StringTie for sample {wildcards.sample_id} on genome {wildcards.genome} successfully completed!"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error occurred during stringTie run: {e}\n")
            logger.error(f"Error occurred during stringTie run: {e}")
            raise e

def get_input_for_TEChimericTranscripts(wildcards):
    logger.info(f"[get_input_for_TEChimericTranscripts] called with wildcards: {wildcards}")
    gtf_path = outdir + f"/{wildcards.genome}/raw/{wildcards.sample_id}/{wildcards.sample_id}.gtf"
    te_gtf = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('TE_gtf')
    if not te_gtf or not os.path.exists(te_gtf):
        raise ValueError(f"TE GTF file for genome {wildcards.genome} is not specified or does not exist in the configuration.")
    in_dict = {
        "gtf": gtf_path,
        "te_gtf": te_gtf
    }
    return in_dict

rule TEChimericTranscripts:
    input:
        unpack(get_input_for_TEChimericTranscripts)
    output:
        txt = outdir + "/{genome}/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt"
    log:
        logdir + "/{sample_id}/{genome}/TEChimericTranscripts.log"
    params:
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
            sample_outdir = os.path.dirname(output.txt)
            script = f"{sample_outdir}/TEChimericTranscripts.{current_time}.sh"
            cmd = [
                "python", params.TEChimericTranscripts,
                "-s", input.gtf,
                "-t", input.te_gtf,
                "-o", output.txt
            ]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "TEChimericTranscripts for sample {wildcards.sample_id} on genome {wildcards.genome} successfully completed!"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error occurred during TEChimericTranscripts: {e}\n")
            logger.error(f"Error occurred during TEChimericTranscripts: {e}")
            raise e

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
        logdir_combine + "/stringtie/{genome}/TEChimericPlot.log"
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
            sample_outdir = os.path.dirname(output.group_stack)
            script = f"{sample_outdir}/TEChimericPlot.{current_time}.sh"
            rule_logger.info(f"Start TEChimericPlot run at {current_time}")
            group_tsv = outdir + f"/{wildcards.genome}/TE_chimeric/sample_groups.tsv"
            if genome_samples.get(wildcards.genome) is None or len(genome_samples[wildcards.genome]) == 0:
                raise ValueError(f"No samples found for genome {wildcards.genome} in genome_samples configuration.")
            with open(group_tsv, 'w') as f:
                f.write("sample\tgroup\n")
                for group, sample_list in sample_groups.items():
                    for sample_id in genome_samples.get(wildcards.genome, []):
                        f.write(f"{sample_id}\t{group}\n")
            
            cmd = [
                "python", params.TEChimericPlot,
                "-i", os.path.join(outdir, wildcards.genome, "raw"),
                "-g", group_tsv,
                "-o", f"{outdir}/{wildcards.genome}/TE_chimeric/TE_chimeric"
            ]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(' '.join(cmd) + "\n")
                f.write(f'echo "TEChimericPlot for genome {wildcards.genome} successfully completed!"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error occurred during TEChimericPlot: {e}\n")
            logger.error(f"Error occurred during TEChimericPlot: {e}")
            raise e


def get_input_for_stringTieMerge(wildcards):
    logger.info(f"[get_input_for_stringTieMerge] called with wildcards: {wildcards}")
    gtfs = []
    for sample_id in genome_samples.get(wildcards.genome, []):
        gtfs.append(outdir + f"/{wildcards.genome}/raw/{sample_id}/{sample_id}.gtf")
    if len(gtfs) == 0:
        raise ValueError(f"No GTF files found for genome {wildcards.genome}.")
    gtf = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"GTF file for genome {wildcards.genome} is not specified or does not exist in the configuration.")
    in_dict = {
        "gtfs": gtfs,
        "gtf": gtf
    }
    return in_dict

rule stringTieMerge:
    input:
        unpack(get_input_for_stringTieMerge)
    output:
        gtf = outdir + "/{genome}/stringtie_merged.gtf"
    log:
        logdir_combine + "/stringtie/{genome}/stringTieMerge.log"
    params:
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
            sample_outdir = os.path.dirname(output.gtf)
            script = os.path.join(sample_outdir, f"stringTieMerge_{current_time}.sh")
            cmd = [params.stringtie, "--merge"] + list(input.gtfs) + ["-o", output.gtf, "-G", input.gtf]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "stringTieMerge for genome {wildcards.genome} successfully completed!"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error occurred during stringTieMerge: {e}\n")
            logger.error(f"Error occurred during stringTieMerge: {e}")
            raise e
