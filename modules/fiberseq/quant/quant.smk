include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_quant")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
p = config.get("Params", {}).get("fiberseq_quant", {})

wildcard_constraints:
    hap = "hap[12]",


rule fiberseq_quant_run:
    """Per-region percent accessibility / occupancy / m6A / FIRE fraction."""
    input:
        peaks = indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
        m6a = indir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        fire = indir + "/{sample_id}/{sample_id}.fire.bed",
    output:
        quant = outdir + "/{sample_id}/{sample_id}.region_quant.tsv",
    log:
        logdir + "/{sample_id}/fiberseq_quant_run.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/quant_matrix.py",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_quant_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start quantification for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.quant))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_quant_run_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script, "run",
                "-p", str(input.peaks),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-a", str(input.m6a),
                "-f", str(input.fire),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_quant_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_quant_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_quant_run for {wildcards.sample_id}: {e}\n")


rule fiberseq_quant_hap:
    """Per-region quantification on a haplotype-split element set."""
    input:
        peaks = indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        msp = indir + "/{sample_id}/{sample_id}.{hap}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.{hap}.nuc.bed.gz",
        m6a = indir + "/{sample_id}/{sample_id}.{hap}.m6a.bed.gz",
    output:
        quant = outdir + "/{sample_id}/{sample_id}.{hap}.region_quant.tsv",
    log:
        logdir + "/{sample_id}/fiberseq_quant_hap_{hap}.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/quant_matrix.py",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(f"fiberseq_quant_hap_{wildcards.hap}", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start hap quantification for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.quant))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_quant_hap_{wildcards.hap}_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script, "run",
                "-p", str(input.peaks),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-a", str(input.m6a),
                "-o", sample_outdir,
                "--sample", f"{wildcards.sample_id}.{wildcards.hap}",
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_quant_hap for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_quant_hap for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_quant_hap for {wildcards.sample_id}: {e}\n")


rule fiberseq_quant_merge:
    """Merge per-sample quant tables into regions x samples matrices."""
    input:
        quants = expand(outdir + "/{sample_id}/{sample_id}.region_quant.tsv",
                        sample_id=samples),
    output:
        accessible = outdir + "/matrix/percent_accessible_matrix.tsv",
        nuc = outdir + "/matrix/nuc_occupancy_matrix.tsv",
        m6a = outdir + "/matrix/m6a_per_fiber_matrix.tsv",
        fire = outdir + "/matrix/fire_frac_matrix.tsv",
    log:
        logdir + "/group/fiberseq_quant_merge.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/quant_matrix.py",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_quant_merge", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start quant merge at {current_time}")
            group_outdir = outdir + "/matrix"
            os.makedirs(group_outdir, exist_ok=True)
            script = os.path.join(group_outdir, f"fiberseq_quant_merge_{current_time}.sh")
            cmd = [
                "python", params.script, "merge",
                "-q", ",".join(str(x) for x in input.quants),
                "-o", group_outdir,
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_quant_merge at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_quant_merge: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_quant_merge: {e}\n")


rule fiberseq_quant_result:
    input:
        quants = expand(outdir + "/{sample_id}/{sample_id}.region_quant.tsv",
                        sample_id=samples),
        matrix = outdir + "/matrix/percent_accessible_matrix.tsv",
