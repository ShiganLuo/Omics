include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_haplotype")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
p = config.get("Params", {}).get("fiberseq_haplotype", {})

# upstream contract: per-hap quant tables from modules/fiberseq_quant
QUANT_DIR = config.get("quant_dir") or indir + "/fiberseq_quant"


rule fiberseq_haplotype_diff:
    """Haplotype-specific accessibility differences (Fisher per region)."""
    input:
        hap1 = QUANT_DIR + "/{sample_id}/{sample_id}.hap1.region_quant.tsv",
        hap2 = QUANT_DIR + "/{sample_id}/{sample_id}.hap2.region_quant.tsv",
    output:
        diff = outdir + "/{sample_id}/{sample_id}.haplotype_diff.tsv",
        summary = outdir + "/{sample_id}/{sample_id}.haplotype_summary.tsv",
    log:
        logdir + "/{sample_id}/fiberseq_haplotype_diff.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/haplotype_diff.py",
        min_fibers = p.get("min_fibers", 10),
        alpha = p.get("alpha", 0.05),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_haplotype_diff", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start haplotype diff for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.diff))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_haplotype_diff_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-1", str(input.hap1),
                "-2", str(input.hap2),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-i", str(params.min_fibers),
                "-a", str(params.alpha),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_haplotype_diff for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_haplotype_diff for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_haplotype_diff for {wildcards.sample_id}: {e}\n")


rule fiberseq_haplotype_result:
    input:
        diff = expand(outdir + "/{sample_id}/{sample_id}.haplotype_diff.tsv",
                      sample_id=samples),
