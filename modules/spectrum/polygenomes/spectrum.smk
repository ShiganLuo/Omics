include: "../../common/common.smk"

import time
import os

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
ROOT_DIR = config.get("ROOT_DIR", ".")
experiment_somatic_vcf_dict = config.get("sample_somatic_vcf_dict", {})
somatic_mutation_vcf_files = list(experiment_somatic_vcf_dict.values())
experiment_group_dict = config.get("sample_group_dict", {})
genomes = sorted((config.get("genome", {}) or {}).keys())

rule somatic_spectrum:
    input:
        somatic_mutation_vcf_files
    output:
        spectrum_png = outdir + "/{genome}/somatic_spectrum_stacked_bar.png",
        spectrum_heatmap = outdir + "/{genome}/somatic_spectrum_heatmap.png",
        spectrum_csv = outdir + "/{genome}/somatic_spectrum.csv"
    log:
        logdir + "/{genome}/spectrum/somatic_spectrum.log"
    threads: rule_res("somatic_spectrum").get("threads", 1)
    resources:
        mem_mb = rule_res("somatic_spectrum").get("mem_mb", 8000),
    params:
        spectrum_script = os.path.join(ROOT_DIR, "modules/spectrum/bin/spectrum.py"),
        outprefix = lambda wildcards: outdir + f"/{wildcards.genome}/somatic_spectrum"
    conda:
        "../spectrum.yaml"
    container:
        sif("../spectrum.yaml", "somatic_spectrum")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("somatic_spectrum", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start somatic spectrum for genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.dirname(str(output.spectrum_png))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"somatic_spectrum_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "somatic_spectrum.sh")
            experiment_vcf_map_file = os.path.join(sample_outdir, "experiment_vcf_map.tsv")
            with open(experiment_vcf_map_file, "w") as f:
                f.write("experiment_sample_id\tvcf\n")
                for experiment_id, vcf_file in experiment_somatic_vcf_dict.items():
                    f.write(f"{experiment_id}\t{vcf_file}\n")
            experiment_group_map_file = os.path.join(sample_outdir, "experiment_group_map.tsv")
            with open(experiment_group_map_file, "w") as f:
                f.write("experiment_sample_id\tgroup\n")
                for experiment_id, group_id in experiment_group_dict.items():
                    f.write(f"{experiment_id}\t{group_id}\n")
            fasta = config.get('genome', {}).get(wildcards.genome, {}).get('fasta')
            if not fasta:
                raise ValueError(f"Missing reference fasta: set 'genome.{wildcards.genome}.fasta' in the workflow config JSON, then rerun.")
            if not os.path.exists(fasta):
                raise ValueError(f"Reference fasta file not found: {fasta}")
            cmd = ["python", params.spectrum_script,
                   "--vcf-map", experiment_vcf_map_file,
                   "--group-map", experiment_group_map_file,
                   "--ref-fasta", fasta,
                   "--output-prefix", params.outprefix
                   ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write('echo "somatic spectrum completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"Somatic spectrum completed for genome {wildcards.genome}")
        except Exception as e:
            rule_logger.error(f"Error occurred during somatic spectrum for genome {wildcards.genome}: {e}\n")
            raise RuntimeError(f"Error occurred during somatic spectrum for genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.\n")

rule spectrum_result:
    input:
        spectrum_png = expand(outdir + "/{genome}/somatic_spectrum_stacked_bar.png", genome=genomes),
        spectrum_heatmap = expand(outdir + "/{genome}/somatic_spectrum_heatmap.png", genome=genomes),
        spectrum_csv = expand(outdir + "/{genome}/somatic_spectrum.csv", genome=genomes)
