include: "../common/common.smk"

from snakemake.logging import logger
import time
import os
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
ROOT_DIR = config.get("ROOT_DIR", ".")
experiment_somatic_vcf_dict = config.get("sample_somatic_vcf_dict", {})
somatic_mutation_vcf_files = experiment_somatic_vcf_dict.values()
experiment_group_dict = config.get("sample_group_dict", {})

rule somatic_spectrum:
    input: 
        somatic_mutation_vcf_files
    output:
        spectrum_png = outdir + "/somatic_spectrum_stacked_bar.png"
    log:
        logdir + "/all/spectrum/somatic_spectrum.log"
    params:
        spectrum_script = os.path.join(ROOT_DIR, "modules/spectrum/bin/spectrum.py"),
        outprefix = outdir + "/somatic_spectrum"
    conda:
        "spectrum.yaml"
    run:
        try:
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = f"{outdir}/somatic_spectrum.{current_time}.sh"
            experiment_vcf_map_file = os.path.join(outdir, "experiment_vcf_map.tsv")
            with open(experiment_vcf_map_file, "w") as f:
                f.write("experiment_sample_id\tvcf\n")
                for experiment_id, vcf_file in experiment_somatic_vcf_dict.items():
                    f.write(f"{experiment_id}\t{vcf_file}\n")
            experiment_group_map_file = os.path.join(outdir, "experiment_group_map.tsv")
            with open(experiment_group_map_file, "w") as f:
                f.write("experiment_sample_id\tgroup\n")
                for experiment_id, group_id in experiment_group_dict.items():
                    f.write(f"{experiment_id}\t{group_id}\n")
            fasta = config.get('genome', {}).get('fasta')
            if not os.path.exists(fasta):
                raise ValueError(f"Reference fasta file not found: {fasta}")
            cmd = ["python", params.spectrum_script, 
                    "--vcf-map", experiment_vcf_map_file,
                    "--group-map", experiment_group_map_file,
                    "--ref-fasta", fasta,
                    "--output-prefix", params.outprefix
                    ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell("bash {script} > {log} 2>&1")
        except Exception as e:
            with open(log) as f:
                f.write(f"Error in somatic_spectrum rule: {str(e)}\n")
            raise e
