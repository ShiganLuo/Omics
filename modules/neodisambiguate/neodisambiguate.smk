include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
# wait for test
def get_inputFile_for_neodisambiguate(wildcards):
    logger.info(f"[get_inputFile_for_neodisambiguate] called with wildcards: {wildcards}")

    return {
        "bam1": f"{indir}/{wildcards.genomeA}/{wildcards.sample_id}.bam",
        "bam2": f"{indir}/{wildcards.genomeB}/{wildcards.sample_id}.bam"
    }
rule neodisambiguate:
    input:
        **get_inputFile_for_neodisambiguate()
    output:
        clean_bam1 = outdir + "/{sample_id}/{sample_id}.{genomeA}.neodisambiguatedA.bam", # for wildcards, the order of genomeA and genomeB is not fixed, so we use .neodisambiguatedA and .neodisambiguatedB to represent the two output bam files
        clean_bam2 = outdir + "/{sample_id}/{sample_id}.{genomeB}.neodisambiguatedB.bam",
        ambiguous_bam1 = outdir + "/{sample_id}/ambiguous-alignments/{sample_id}.{genomeA}.ambiguousA.bam",
        ambiguous_bam2 = outdir + "/{sample_id}/ambiguous-alignments/{sample_id}.{genomeB}.ambiguousB.bam",
    params:
        prefix = lambda wildcards: f"{outdir}/{wildcards.sample_id}.neodisambiguated",
        neodisambiguate = config.get("Procedure", {}).get("neodisambiguate") or "neodisambiguate"
    threads: 8
    log:
        logdir + "/{sample_id}/neodisambiguate.log"
    conda:
        "neodisambiguate.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="neodisambiguate", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start neodisambiguate for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/neodisambiguate_{current_time}.sh")
            cmd = [
                str(params.neodisambiguate),
                "-s", str(params.prefix),
                "-o", str(outdir),
                str(input.bam1),
                str(input.bam2)
            ]
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