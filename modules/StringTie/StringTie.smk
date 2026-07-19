include: "../common/common.smk"
from snakemake.logging import logger
import time
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
        gtf = outdir + "/{sample_id}/{sample_id}.gtf"
    log:
        logdir + "/{sample_id}/stringTie.log"
    params:
        gtf = config.get('genome', {}).get('gtf'), #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    threads: 5
    conda:
        "StringTie.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="stringTie_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start stringTie run for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/stringTie_{current_time}.sh")
            cmd = [params.stringtie, "-o", output.gtf, input.bam, "-G", params.gtf, "-p", str(threads)]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")
rule TEChimericTranscripts:
    input:
        gtf = outdir + "/{sample_id}/{sample_id}.gtf"
    output:
        txt = outdir + "/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt"
    log:
        logdir + "/{sample_id}/TEChimericTranscripts.log"
    params:
        te_gtf = config.get('genome', {}).get('TE_gtf'),
        TEChimericTranscripts = ROOT_DIR + "/modules/StringTie/bin/TEChimericTranscripts.py"
    threads: 5
    conda:
        "StringTie.yaml"
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/{wildcards.sample_id}/TEChimericTranscripts.{current_time}.sh"
        cmd = f"python {params.TEChimericTranscripts} -s {input.gtf} -t {params.te_gtf} -o {output.txt} > {log} 2>&1"
        with open(script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(cmd + "\n")
        shell(f"bash {script}")

rule TEChimericPlot:
    input:
        txts = expand(outdir + "/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt", sample_id=samples)
    output:
        group_stack = outdir + "/result/TE_chimeric/TE_chimeric_group_stacked.png",
        type_top = outdir + "/result/TE_chimeric/TE_chimeric_te_type_top.png",
        type_by_group = outdir + "/result/TE_chimeric/TE_chimeric_te_type_by_group.png",
        sample_summary = outdir + "/result/TE_chimeric/TE_chimeric_sample_summary.tsv",
        group_summary = outdir + "/result/TE_chimeric/TE_chimeric_group_summary.tsv",
        te_type_counts = outdir + "/result/TE_chimeric/TE_chimeric_te_type_counts.tsv"
    log:
        logdir + "/all/stringtie/TEChimericPlot.log"
    params:
        TEChimericPlot = ROOT_DIR + "/modules/StringTie/bin/TEChimericPlot.py"
    threads: 1
    conda:
        "StringTie.yaml"
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/result/TE_chimeric/TEChimericPlot.{current_time}.sh"
        group_tsv = outdir + "/result/TE_chimeric/sample_groups.tsv"
        with open(group_tsv, 'w') as f:
            f.write("sample\tgroup\n")
            for group, sample_list in sample_groups.items():
                for sample_id in sample_list:
                    f.write(f"{sample_id}\t{group}\n")
        cmd = f"python {params.TEChimericPlot} -i {outdir} -g {group_tsv} -o {outdir}/result/TE_chimeric/TE_chimeric > {log} 2>&1"
        with open(script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(cmd + "\n")
        shell(f"bash {script}")


rule stringTieMerge:
    input:
        gtfs = expand(outdir + "/{sample_id}/{sample_id}.gtf",sample_id=samples)
    output:
        gtf = outdir + "/stringtie_merged.gtf"
    log:
        logdir + "/all/stringtie/stringTieMerge.log"
    params:
        gtf = config.get('genome', {}).get('gtf'), #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    conda:
        "StringTie.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="stringTieMerge_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start stringTieMerge run at {current_time}")
            script = os.path.join(outdir, f"stringTieMerge_{current_time}.sh")
            cmd = [params.stringtie, "--merge"] + list(input.gtfs) + ["-o", output.gtf, "-G", params.gtf]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")
