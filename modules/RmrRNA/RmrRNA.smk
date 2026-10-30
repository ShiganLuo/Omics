include: "../common/common.smk"

import time
import os
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
rule extract_rRNA:
    input:
        fasta = config.get('genome', {}).get('fasta'),
        gtf = config.get('genome', {}).get('gtf')
    output:
        rRNA_fasta = outdir + "/rRNA.fasta"
    log:
        logdir + "/../group/RmrRNA/extract_rRNA.log"
    threads: 2
    params:
        extract_rRNA_script = os.path.join(ROOT_DIR, "modules/RmrRNA/bin/extract_rRNA.py")
    conda:
        "RmrRNA.yaml"
    container:
        sif("RmrRNA.yaml")
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = os.path.join(RECEIPT_DIR, f"extract_rRNA_{current_time}.sh")
        plain_script = f"{outdir}/extract_rRNA.sh"
        cmd = ["python", params.extract_rRNA_script, 
                "--fasta", input.fasta, 
                "--gtf", input.gtf, 
                "--output", output.rRNA_fasta,
                "--threads", str(threads)
                ]
        with open(script, "w") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(" ".join(cmd) + "\n")
            f.write(f'echo "script completed successfully"\n')
        shell("bash {plain_script} > {log} 2>&1")

