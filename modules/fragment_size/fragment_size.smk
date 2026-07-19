include: "../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

rule samtools_stats:
    """Get BAM statistics including fragment size distribution."""
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        stats = outdir + "/stats/{sample_id}/{sample_id}_stats.txt"
    log:
        logdir + "/{sample_id}/samtools_stats.log"
    conda:
        "fragment_size.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("samtools_stats", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start samtools_stats for sample {wildcards.sample_id} at {current_time}")
            outdir_sample = os.path.dirname(str(output.stats))
            os.makedirs(outdir_sample, exist_ok=True)
            script = os.path.join(outdir_sample, f"samtools_stats_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.samtools} stats {input.bam} > {output.stats}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during samtools_stats for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during samtools_stats for sample {wildcards.sample_id}: {e}")
            raise e

rule getFragmentSize:
    """Extract and summarize fragment size distribution from all samples."""
    input:
        stats = expand(outdir + "/stats/{sample_id}/{sample_id}_stats.txt", sample_id=samples)
    output:
        hist = outdir + "/fragment/FragmentSize.txt"
    log:
        logdir + "/getFragmentSize.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/fragment_size/bin/getFragmentSize.py")
    conda:
        "fragment_size.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("getFragmentSize", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start getFragmentSize at {current_time}")
            outdir_fragment = os.path.dirname(str(output.hist))
            os.makedirs(outdir_fragment, exist_ok=True)
            script = os.path.join(outdir_fragment, f"getFragmentSize_{current_time}.sh")
            stats_str = " ".join(input.stats)
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python {params.script} --input {stats_str} --out {output.hist}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during getFragmentSize: {e}\n")
            logger.error(f"Error occurred during getFragmentSize: {e}")
            raise e

rule plotFragmentSize:
    """Plot fragment size distribution."""
    input:
        hist = outdir + "/fragment/FragmentSize.txt"
    output:
        png = outdir + "/fragment/FragmentSize.png"
    log:
        logdir + "/plotFragmentSize.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/fragment_size/bin/plotFragmentSize.py")
    conda:
        "fragment_size.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("plotFragmentSize", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start plotFragmentSize at {current_time}")
            outdir_fragment = os.path.dirname(str(output.png))
            script = os.path.join(outdir_fragment, f"plotFragmentSize_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python {params.script} --input {input.hist} --output {output.png}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during plotFragmentSize: {e}\n")
            logger.error(f"Error occurred during plotFragmentSize: {e}")
            raise e

rule fragment_size_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        hist = outdir + "/fragment/FragmentSize.txt",
        png = outdir + "/fragment/FragmentSize.png"
