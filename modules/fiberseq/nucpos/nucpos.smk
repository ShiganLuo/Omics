include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_nucpos")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
p = config.get("Params", {}).get("fiberseq_nucpos", {})


rule fiberseq_nucpos_run:
    """Nucleosome positioning/phasing around TSS anchors."""
    input:
        gtf = lambda wildcards: _require(genome_cfg.get("gtf"), "genome.gtf", "fiberseq_nucpos"),
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        per_anchor = outdir + "/{sample_id}/{sample_id}.nuc_offset_per_anchor.tsv",
        summary = outdir + "/{sample_id}/{sample_id}.nuc_offset_summary.tsv",
        pdf = outdir + "/{sample_id}/{sample_id}.nucpos.pdf",
    log:
        logdir + "/{sample_id}/fiberseq_nucpos_run.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/nucpos.py",
        upstream = p.get("upstream_nucs", 5),
        downstream = p.get("downstream_nucs", 5),
        min_fibers = p.get("min_fibers", 10),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_nucpos_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start nucpos for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.per_anchor))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_nucpos_run_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-g", str(input.gtf),
                "--anchor", "tss",
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-5", str(params.upstream),
                "-3", str(params.downstream),
                "-i", str(params.min_fibers),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_nucpos_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_nucpos_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_nucpos_run for {wildcards.sample_id}: {e}\n")


rule fiberseq_nucpos_result:
    input:
        summary = expand(outdir + "/{sample_id}/{sample_id}.nuc_offset_summary.tsv",
                         sample_id=samples),
