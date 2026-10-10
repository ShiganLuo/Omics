include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_censat")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
p = config.get("Params", {}).get("fiberseq_censat", {})


rule fiberseq_censat_run:
    """Region-class accessibility and methylation profiles (CDR/HOR/flank)."""
    input:
        regions = lambda wildcards: _require(genome_cfg.get("censat_bed"), "genome.censat_bed", "fiberseq_censat"),
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
        cpg = indir + "/{sample_id}/{sample_id}.cpg.bed.gz",
    output:
        regions_tsv = outdir + "/{sample_id}/{sample_id}.censat_regions.tsv",
        class_summary = outdir + "/{sample_id}/{sample_id}.censat_class_summary.tsv",
        profile = outdir + "/{sample_id}/{sample_id}.censat_profile.tsv",
        pdf = outdir + "/{sample_id}/{sample_id}.censat_profile.pdf",
    log:
        logdir + "/{sample_id}/fiberseq_censat_run.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/censat_profile.py",
        pileup = genome_cfg.get("mcpg_pileup"),
        msp_min_len = p.get("msp_min_len", 50),
        flank = p.get("flank", 10000),
        bin = p.get("bin", 500),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_censat_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start censat profile for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.regions_tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_censat_run_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-r", str(input.regions),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-c", str(input.cpg),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-m", str(params.msp_min_len),
                "-f", str(params.flank),
                "-b", str(params.bin),
                "-l", log_path,
            ]
            if params.pileup:
                cmd += ["-p", str(params.pileup)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_censat_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_censat_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_censat_run for {wildcards.sample_id}: {e}\n")


rule fiberseq_censat_result:
    input:
        class_summary = expand(outdir + "/{sample_id}/{sample_id}.censat_class_summary.tsv",
                               sample_id=samples),
