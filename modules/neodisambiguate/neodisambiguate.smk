include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

def get_inputFile_for_neodisambiguate(wildcards):
    logger.info(f"[get_inputFile_for_neodisambiguate] called with {wildcards}")

    return {
        "bam1": f"{indir}/{wildcards.genomeA}/{wildcards.sample_id}.bam",
        "bam2": f"{indir}/{wildcards.genomeB}/{wildcards.sample_id}.bam"
    }

rule neodisambiguate:
    input:
        unpack(get_inputFile_for_neodisambiguate)
    output:
        clean_bam1 = outdir + "/{sample_id}/{genomeA}_{genomeB}/{sample_id}.neodisambiguated.{genomeA}.bam",
        clean_bam2 = outdir + "/{sample_id}/{genomeA}_{genomeB}/{sample_id}.neodisambiguated.{genomeB}.bam",
        ambiguous_bam1 = outdir + "/{sample_id}/{genomeA}_{genomeB}/ambiguous-alignments/{sample_id}.{genomeA}.ambiguous.bam",
        ambiguous_bam2 = outdir + "/{sample_id}/{genomeA}_{genomeB}/ambiguous-alignments/{sample_id}.{genomeB}.ambiguous.bam",
    params:
        pair_dir = lambda wildcards: f"{outdir}/{wildcards.sample_id}/{wildcards.genomeA}_{wildcards.genomeB}",
        prefix = lambda wildcards: f"{outdir}/{wildcards.sample_id}/{wildcards.genomeA}_{wildcards.genomeB}/{wildcards.sample_id}.neodisambiguated",
        neodisambiguate = config.get("Procedure", {}).get("neodisambiguate") or "neodisambiguate"
    threads: rule_res("neodisambiguate").get("threads", 8)
    benchmark: logdir + "/{sample_id}/{genomeA}_{genomeB}/neodisambiguate.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "neodisambiguate", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    log:
        logdir + "/{sample_id}/{genomeA}_{genomeB}/neodisambiguate.log"
    conda:
        "neodisambiguate.yaml"
    container:
        sif("neodisambiguate.yaml", "neodisambiguate")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        logger = setup_logger(logger_name="neodisambiguate", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start neodisambiguate for sample {wildcards.sample_id} ({wildcards.genomeA} vs {wildcards.genomeB}) at {current_time}")
            script = receipt_script(f"neodisambiguate_{wildcards.sample_id}_{wildcards.genomeA}_{wildcards.genomeB}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            pair_dir = str(params.pair_dir)
            staged_dir = os.path.join(pair_dir, "inputs")
            os.makedirs(staged_dir, exist_ok=True)
            os.makedirs(os.path.join(pair_dir, "ambiguous-alignments"), exist_ok=True)
            staged = []
            for genome, src in ((wildcards.genomeA, str(input.bam1)), (wildcards.genomeB, str(input.bam2))):
                dst = os.path.join(staged_dir, f"{wildcards.sample_id}.{genome}.bam")
                if os.path.lexists(dst):
                    os.remove(dst)
                os.symlink(os.path.abspath(src), dst)
                staged.append(dst)
            plain_script = os.path.join(pair_dir, "neodisambiguate.sh")
            cmd = [
                "bash", f'"$(command -v {params.neodisambiguate})"',
                "-i", staged[0], staged[1],
                "-o", str(params.prefix),
                "-n", str(wildcards.genomeA), str(wildcards.genomeB)
            ]
            with open(script, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise RuntimeError(f"Error: {e}")
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")
