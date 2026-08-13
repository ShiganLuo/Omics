include: "../common/common.smk"
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
samples = config.get("samples",[])
bam_substring = config.get("bam_substring") or ""
FUSION_FIGS = [
    "fig1_per_sample_counts.png",
    "fig2_fusion_type_distribution.png",
    "fig3_type_by_sample.png",
    "fig4_reading_frame.png",
    "fig5_recurrent_heatmap.png",
    "fig6_inframe_support.png",
]

def get_input_for_arriba(wildcards):
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
    container:
        sif("arriba.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("arriba", log_file=log_path)
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = f"{outdir}/{wildcards.sample_id}/{wildcards.sample_id}_arriba.{current_time}.sh"
            if not os.path.exists(input.blacklist):
                rule_logger.error(f"Blacklist file not found: {input.blacklist}")
                raise ValueError(f"Blacklist file not found: {input.blacklist}")
            if not os.path.exists(input.known_fusions):
                rule_logger.error(f"Known fusions file not found: {input.known_fusions}")
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
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as fh:
                fh.write(f"Arriba fusion detection failed: {e}\n")
            raise f"Arriba fusion detection failed: {e}\n"

rule arriba_report:
    input:
        passed_fusions = expand(outdir + "/{sample_id}/{sample_id}_passed_fusions.tsv", sample_id=samples),
        discarded_fusions = expand(outdir + "/{sample_id}/{sample_id}_discarded_fusions.tsv", sample_id=samples)
    output:
        report = outdir + "/arriba_report/arriba_fusion_report.html",
        fusion_summary = outdir + "/arriba_report/per_sample_summary.tsv",
        recurrent_fusions = outdir + "/arriba_report/recurrent_fusions.tsv",
        high_medium_fusions = outdir + "/arriba_report/high_medium_confidence_fusions.tsv",
        inframe_fusions = outdir + "/arriba_report/inframe_fusions.tsv",
        fusion_figs = expand(outdir + "/arriba_report/figures/{fig}", fig=FUSION_FIGS),
    log:
        logdir + "/../group/arriba/arriba_report.log"
    params:
        summary_script = os.path.join(ROOT_DIR, "modules/arriba/bin/summarize_arriba_fusions.py")
    conda:
        "arriba.yaml"
    container:
        sif("arriba.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("arriba_report", log_file=log_path)
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = f"{outdir}/arriba_report.{current_time}.sh"
            cmd = [
                "python", params.summary_script,
                "-p", ",".join(input.passed_fusions),
                "-d", ",".join(input.discarded_fusions),
                "-o", outdir + "/arriba_report"
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as fh:
                fh.write(f"Arriba report generation failed: {e}\n")
            raise f"Arriba report generation failed: {e}\n"
rule arriba_result:
    input:
        passed_fusion_tsv = outdir + "/{sample_id}/{sample_id}_passed_fusions.tsv",
        discarded_fusion_tsv = outdir + "/{sample_id}/{sample_id}_discarded_fusions.tsv"
