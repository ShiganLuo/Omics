from snakemake.logging import logger
import time
import os
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
ROOT_DIR = config.get("ROOT_DIR", ".")

def get_input_for_somatic_spectrum(wildcards):
    vcf_files = []
    for normal_sample_id in config.get("normal_samples", []):
        for experimental_sample_id in config.get("experimental_samples", []):
            vcf_path = f"{outdir}/mutect2-vcf/{normal_sample_id}_vs_{experimental_sample_id}/{normal_sample_id}_vs_{experimental_sample_id}.vcf.gz"
            vcf_files.append(vcf_path)
    return vcf_files

rule somatic_spectrum:
    input: 
        get_input_for_somatic_spectrum
    output:
        spectrum_png = outdir + "/mutect2-vcf/spectrum/somatic_spectrum.png"
    log:
        logdir + "/all/gatk/somatic_spectrum.log"
    conda:
        "../gatk.yaml"
    params:
        spectrum_script = os.path.join(ROOT_DIR, "modules/spectrum/bin/spectrum.py")
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/somatic_spectrum.{current_time}.sh"
        cmd = ["python", params.spectrum_script, 
                "--vcf_files"] + input + [
                "--output", output.spectrum_png
                ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")
