include: "../common/common.smk"

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
samples = config.get("samples") or []
sample_groups = config.get("sample_groups") or {}
ROOT_DIR = config.get("ROOT_DIR", ".")
rule stringTie:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        gtf = outdir + "/raw/{sample_id}/{sample_id}.gtf"
    log:
        logdir + "/{sample_id}/stringTie.log"
    params:
        gtf = config.get('genome', {}).get('gtf'), #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    threads: 5
    conda:
        "StringTie.yaml"
    container:
        sif("StringTie.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger(logger_name="stringTie_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start stringTie run for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"stringTie_{wildcards.sample_id}_{current_time}.sh")
            plain_script = f"...{wildcards.sample_id}/stringTie.sh"
            cmd = [params.stringtie, "-o", output.gtf, input.bam, "-G", params.gtf, "-p", str(threads)]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f"echo 'stringTie run completed for sample {wildcards.sample_id} was successful completed'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during stringTie run for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during stringTie run for sample {wildcards.sample_id}: {e}\n")

rule TEChimericTranscripts:
    input:
        gtf = outdir + "/raw/{sample_id}/{sample_id}.gtf"
    output:
        txt = outdir + "/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt"
    log:
        logdir + "/{sample_id}/TEChimericTranscripts.log"
    params:
        te_gtf = config.get('genome', {}).get('TE_gtf'),
        TEChimericTranscripts = ROOT_DIR + "/modules/StringTie/bin/TEChimericTranscripts.py"
    threads: 5
    conda:
        "StringTie.yaml"
    container:
        sif("StringTie.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger(logger_name="TEChimericTranscripts_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            rule_logger.info(f"Start TEChimericTranscripts run for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"TEChimericTranscripts_{wildcards.sample_id}_{current_time}.sh")
            plain_script = f"...{wildcards.sample_id}/TEChimericTranscripts.sh"
            cmd = f"python {params.TEChimericTranscripts} -s {input.gtf} -t {params.te_gtf} -o {output.txt} > {log} 2>&1"
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(cmd + "\n")
                f.write(f"echo 'TEChimericTranscripts run completed for sample {wildcards.sample_id} was successful completed'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during TEChimericTranscripts for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during TEChimericTranscripts for sample {wildcards.sample_id}: {e}\n")

rule TEChimericPlot:
    input:
        txts = expand(outdir + "/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt", sample_id=samples)
    output:
        group_stack = outdir + "/TE_chimeric/TE_chimeric_group_stacked.png",
        type_top = outdir + "/TE_chimeric/TE_chimeric_te_type_top.png",
        type_by_group = outdir + "/TE_chimeric/TE_chimeric_te_type_by_group.png",
        sample_summary = outdir + "/TE_chimeric/TE_chimeric_sample_summary.tsv",
        group_summary = outdir + "/TE_chimeric/TE_chimeric_group_summary.tsv",
        te_type_counts = outdir + "/TE_chimeric/TE_chimeric_te_type_counts.tsv"
    log:
        logdir + "/../group/stringtie/TEChimericPlot.log"
    params:
        TEChimericPlot = ROOT_DIR + "/modules/StringTie/bin/TEChimericPlot.py"
    threads: 1
    conda:
        "StringTie.yaml"
    container:
        sif("StringTie.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger(logger_name="TEChimericPlot_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = os.path.join(RECEIPT_DIR, f"TEChimericPlot_{current_time}.sh")
            rule_logger.info(f"Start TEChimericPlot run at {current_time}")
            group_tsv = outdir + "/TE_chimeric/sample_groups.tsv"
            with open(group_tsv, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write("sample\tgroup\n")
                for group, sample_list in sample_groups.items():
                    for sample_id in sample_list:
                        f.write(f"{sample_id}\t{group}\n")
                f.write(f'echo "script completed successfully"\n')
            cmd = [
                "python", params.TEChimericPlot,
                "-i", os.path.join(outdir,"raw"),
                "-g", group_tsv,
                "-o", f"{outdir}/TE_chimeric/TE_chimeric"
            ]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + "\n")
                f.write(f"echo 'TEChimericPlot run completed successfully'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during TEChimericPlot: {e}")
            raise RuntimeError(f"Error occurred during TEChimericPlot: {e}\n")

rule stringTieMerge:
    input:
        gtfs = expand(outdir + "/raw/{sample_id}/{sample_id}.gtf",sample_id=samples)
    output:
        gtf = outdir + "/stringtie_merged.gtf"
    log:
        logdir + "/../group/stringtie/stringTieMerge.log"
    params:
        gtf = config.get('genome', {}).get('gtf'), #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    conda:
        "StringTie.yaml"
    container:
        sif("StringTie.yaml")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger(logger_name="stringTieMerge_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start stringTieMerge run at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"stringTieMerge_{current_time}.sh")
            plain_script = os.path.join(outdir, f"stringTieMerge.sh")
            cmd = [params.stringtie, "--merge"] + list(input.gtfs) + ["-o", output.gtf, "-G", params.gtf]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f"echo 'stringTieMerge run completed successfully'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during stringTieMerge: {e}")
            raise RuntimeError(f"Error occurred during stringTieMerge: {e}\n")
