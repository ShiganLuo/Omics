include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_diff")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
p = config.get("Params", {}).get("fiberseq_diff", {})

# upstream contract: matrices from modules/fiberseq_quant, annotations from
# modules/fiberseq_annotate
QUANT_DIR = config.get("quant_dir") or indir + "/fiberseq_quant"
ANNOT_DIR = config.get("annot_dir") or indir + "/fiberseq_annotate"
ANNOT_SAMPLE = p.get("annot_sample") or (samples[0] if samples else "")


rule fiberseq_diff_run:
    """Per-region differential accessibility for one contrast."""
    input:
        matrix = QUANT_DIR + "/matrix/percent_accessible_matrix.tsv",
        design = lambda wildcards: _require(p.get("design_tsv"),
                                            "Params.diff.design_tsv",
                                            "fiberseq_diff"),
        annot = ANNOT_DIR + "/" + ANNOT_SAMPLE + "/" + ANNOT_SAMPLE + ".peaks_annotated.tsv",
    output:
        diff = outdir + "/{contrast}/diff_{contrast}.tsv",
        volcano = outdir + "/{contrast}/volcano_{contrast}.pdf",
        class_summary = outdir + "/{contrast}/class_summary_{contrast}.tsv",
        gene_table = outdir + "/{contrast}/{contrast}.TEcount_Gene.name.tsv",
    log:
        logdir + "/group/fiberseq_diff_run_{contrast}.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/diff_accessibility.py",
        alpha = p.get("alpha", 0.05),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_diff_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start diff analysis for {wildcards.contrast} at {current_time}")
            contrast_outdir = os.path.dirname(str(output.diff))
            os.makedirs(contrast_outdir, exist_ok=True)
            script = os.path.join(contrast_outdir, f"fiberseq_diff_run_{wildcards.contrast}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.matrix),
                "-d", str(input.design),
                "-a", str(input.annot),
                "-C", str(wildcards.contrast),
                "-o", contrast_outdir,
                "-p", str(params.alpha),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_diff_run for {wildcards.contrast} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_diff_run for {wildcards.contrast}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_diff_run for {wildcards.contrast}: {e}\n")


rule fiberseq_diff_result:
    input:
        matrix = QUANT_DIR + "/matrix/percent_accessible_matrix.tsv",
