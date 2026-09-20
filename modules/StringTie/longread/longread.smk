include: "../../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
upstream_indir = config.get("upstream_indir", config.get("indir", "input"))
merged_gtf = config.get("merged_gtf", outdir + "/stringtie_merged.gtf")

rule stringtie_lr_assemble:
    input:
        bam = upstream_indir + "/{sample_id}/{sample_id}.sorted.bam",
        gtf = config.get("genome", {}).get("gtf") or "/dev/null"
    output:
        gtf = outdir + "/raw/{sample_id}/{sample_id}.gtf"
    log:
        logdir + "/{sample_id}/stringtie_lr_assemble.log"
    threads: 8
    conda:
        "../StringTie.yaml"
    container:
        sif("../StringTie.yaml")
    params:
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("stringtie_lr_assemble", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start StringTie long-read assembly for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.gtf))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"stringtie_lr_assemble_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                params.stringtie, "-L",
                "-p", str(threads),
                "-G", input.gtf,
                "-o", output.gtf,
                input.bam
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "StringTie LR assembly for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during StringTie LR assembly for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error during StringTie LR assembly for sample {wildcards.sample_id}: {e}")
            raise e
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"StringTie LR assembly finished for sample {wildcards.sample_id} at {current_time}")

rule stringtie_lr_merge:
    input:
        gtfs = expand(outdir + "/raw/{sample_id}/{sample_id}.gtf", sample_id=samples)
    output:
        gtf = outdir + "/stringtie_merged.gtf"
    log:
        logdir + "/stringtie_lr_merge.log"
    conda:
        "../StringTie.yaml"
    container:
        sif("../StringTie.yaml")
    params:
        gtf = config.get("genome", {}).get("gtf") or "/dev/null",
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("stringtie_lr_merge", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start StringTie merge at {current_time}")
            script = os.path.join(outdir, f"stringtie_lr_merge_{current_time}.sh")
            cmd = [params.stringtie, "--merge"] + list(input.gtfs) + ["-o", output.gtf, "-G", params.gtf]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "StringTie merge at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during StringTie merge: {e}\n")
            logger.error(f"Error during StringTie merge: {e}")
            raise e
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"StringTie merge finished at {current_time}")

rule stringtie_lr_quant:
    input:
        bam = upstream_indir + "/{sample_id}/{sample_id}.sorted.bam",
        gtf = config.get("genome", {}).get("gtf") or "/dev/null",
        merged_gtf = merged_gtf
    output:
        abund = outdir + "/{sample_id}/{sample_id}.gene_abund.tab",
        gtf = outdir + "/{sample_id}/{sample_id}.quant.gtf"
    log:
        logdir + "/{sample_id}/stringtie_lr_quant.log"
    threads: 4
    conda:
        "../StringTie.yaml"
    container:
        sif("../StringTie.yaml")
    params:
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("stringtie_lr_quant", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start StringTie quantification for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.gtf))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"stringtie_lr_quant_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                params.stringtie, "-e", "-B", "-L",
                "-p", str(threads),
                "-G", input.merged_gtf,
                "-o", output.gtf,
                "-A", output.abund,
                input.bam
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "StringTie quantification for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during StringTie quantification for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error during StringTie quantification for sample {wildcards.sample_id}: {e}")
            raise e
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"StringTie quantification finished for sample {wildcards.sample_id} at {current_time}")

rule stringtie_lr_result:
    input:
        merged_gtf = outdir + "/stringtie_merged.gtf",
        abund = expand(outdir + "/{sid}/{sid}.gene_abund.tab", sid=samples),
        quant_gtf = expand(outdir + "/{sid}/{sid}.quant.gtf", sid=samples)