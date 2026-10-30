include: "../../common/common.smk"

import os

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
smallrna_types = config.get("Params", {}).get("smallrna_types",
    ["miRNA", "snRNA", "snoRNA", "rRNA", "misc_RNA", "scRNA", "scaRNA", "vaultRNA"])
flank = config.get("Params", {}).get("smallrna_flank", 50)
ROOT_DIR = config.get("ROOT_DIR", os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

rule chromosome_sizes:
    input:
        fasta = config.get('genome', {}).get(wildcards.genome, {}).get('fasta')
    output:
        chrom_sizes = outdir + "/{genome}/chrom.sizes"
    log:
        logdir + "/{genome}/chromosome_sizes.log"
    threads: 1
    conda:
        "../genome.yaml"
    container:
        sif("../genome.yaml")
    params:
        samtools = config.get('Procedure', {}).get('samtools') or 'samtools'
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("chromosome_sizes", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start chromosome_sizes at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"chromosome_sizes_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(output.chrom_sizes), f"chromosome_sizes.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"{params.samtools} faidx {input.fasta}\n")
                f.write(f"cut -f1,2 {input.fasta}.fai > {output.chrom_sizes}\n")
                f.write(f"echo 'chromosome_sizes completed for genome {wildcards.genome} at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"chromosome_sizes completed")
        except Exception as e:
            rule_logger.error(f"Error occurred during chromosome_sizes: {e}")
            raise RuntimeError(f"Error occurred during chromosome_sizes: {e}, you can check the log file {log_path} for more details.")

        # ── Extract smallRNA BED + FASTA from GENCODE GTF ───────────────────────
rule extract_smallrna:
    input:
        gtf = config.get("genome", {}).get("gtf"),
        fasta = config.get("genome", {}).get("fasta"),
        chrom_sizes = outdir + "/{genome}/chrom.sizes"
    output:
        bed = outdir + "/genome/smallrna/smallrna_genes.bed",
        fasta = outdir + "/genome/smallrna/smallrna_genes_flank.fa",
    log:
        logdir + "/genome/extract_smallrna.log"
    threads: 1
    conda:
        "../genome.yaml"
    container:
        sif("../genome.yaml")
    params:
        script = os.path.join(ROOT_DIR, "modules", "genome", "bin", "extract_smallrna.py"),
        outdir = outdir + "/genome/smallrna",
        flank = flank,
        types = " ".join(smallrna_types),
        bedtools = config.get("Procedure", {}).get("bedtools") or "bedtools",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("extract_smallrna", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start extract_smallrna at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"extract_smallrna_{current_time}.sh")
            plain_script = os.path.join(params.outdir, f"extract_smallrna.sh")
            os.makedirs(params.outdir, exist_ok=True)
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"python {params.script} \\\n")
                f.write(f"    --gtf {input.gtf} \\\n")
                f.write(f"    --fasta {input.fasta} \\\n")
                f.write(f"    --chrom-sizes {input.chrom_sizes} \\\n")
                f.write(f"    --outdir {params.outdir} \\\n")
                f.write(f"    --flank {params.flank} \\\n")
                f.write(f"    --types {params.types} \\\n")
                f.write(f"    --bedtools {params.bedtools} \\\n")
                f.write(f"echo 'extract_smallrna completed successfully at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during extract_smallrna: {e}")
            raise RuntimeError(f"Error occurred during extract_smallrna: {e}, you can check the log file {log_path} for more details.")