include: "../../common/common.smk"

import time
import os

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
ROOT_DIR = config.get("ROOT_DIR", ".")
genomes = sorted((config.get("genome", {}) or {}).keys())


def _spectrum_cfg():
    """Normalized spectrum config block (single source for input fns + params).

    Keys: vcf_map (sample -> VCF), group_map (sample -> group),
    signature (analysis parameters, shared across genomes).
    """
    cfg = {
        "vcf_map": dict(config.get("sample_somatic_vcf_dict") or {}),
        "group_map": dict(config.get("sample_group_dict") or {}),
        "signature": dict(config.get("signature") or {}),
    }
    cfg["signature"].setdefault("mode", "both")
    cfg["signature"].setdefault("min_mutations", 200)
    cfg["signature"].setdefault("min_exposure_frac", 0.02)
    cfg["signature"].setdefault("family_map", None)
    cfg["signature"].setdefault("max_rank", 8)
    cfg["signature"].setdefault("restarts", 20)
    cfg["signature"].setdefault("seed", 42)
    cfg["signature"].setdefault("rank", None)
    return cfg


def _genome_fasta(genome):
    """Resolve genome.<name>.fasta, failing at build time when unset."""
    fasta = (config.get("genome", {}).get(genome) or {}).get("fasta")
    if not fasta:
        raise ValueError(
            f"Missing reference fasta: set 'genome.{genome}.fasta' in the "
            f"workflow config JSON, then rerun."
        )
    return fasta


def _genome_signature_matrix(genome):
    """Resolve genome.<name>.signature_matrix, failing at build time when unset."""
    matrix = (config.get("genome", {}).get(genome) or {}).get("signature_matrix")
    if not matrix:
        raise ValueError(
            f"Missing COSMIC signature matrix: set 'genome.{genome}.signature_matrix' "
            f"in the workflow config JSON (e.g. /data/16T/Database/cosmic_signatures/"
            f"mm39/COSMIC_v3.6_SBS_mm39.txt), then rerun."
        )
    return matrix


def _spectrum_inputs(wildcards):
    """Declare all inputs of somatic_spectrum; missing key fails at build time."""
    cfg = _spectrum_cfg()
    return {
        "fasta": _genome_fasta(wildcards.genome),
        "vcfs": list(cfg["vcf_map"].values()),
    }


def _signature_inputs(wildcards):
    """Declare inputs of the signature rules; matrix key is required."""
    return {
        "counts": outdir + f"/{wildcards.genome}/somatic_spectrum_per_sample.csv",
        "signatures": _genome_signature_matrix(wildcards.genome),
    }


rule somatic_spectrum:
    input:
        unpack(_spectrum_inputs)
    output:
        spectrum_png = outdir + "/{genome}/somatic_spectrum_stacked_bar.png",
        spectrum_heatmap = outdir + "/{genome}/somatic_spectrum_heatmap.png",
        spectrum_csv = outdir + "/{genome}/somatic_spectrum.csv",
        per_sample = outdir + "/{genome}/somatic_spectrum_per_sample.csv",
        sample_stats = outdir + "/{genome}/somatic_spectrum_sample_stats.csv"
    log:
        logdir + "/{genome}/spectrum/somatic_spectrum.log"
    threads: rule_res("somatic_spectrum").get("threads", 1)
    resources:
        mem_mb = rule_res("somatic_spectrum").get("mem_mb", 8000),
    params:
        spectrum_script = os.path.join(ROOT_DIR, "modules/spectrum/bin/spectrum.py"),
        outprefix = lambda wildcards: outdir + f"/{wildcards.genome}/somatic_spectrum",
        cfg = _spectrum_cfg()
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
                for experiment_id, vcf_file in params.cfg["vcf_map"].items():
                    f.write(f"{experiment_id}\t{vcf_file}\n")
            experiment_group_map_file = os.path.join(sample_outdir, "experiment_group_map.tsv")
            with open(experiment_group_map_file, "w") as f:
                f.write("experiment_sample_id\tgroup\n")
                for experiment_id, group_id in params.cfg["group_map"].items():
                    f.write(f"{experiment_id}\t{group_id}\n")
            cmd = ["python", params.spectrum_script,
                   "--vcf-map", experiment_vcf_map_file,
                   "--group-map", experiment_group_map_file,
                   "--ref-fasta", str(input.fasta),
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


rule somatic_signature_refit:
    input:
        unpack(_signature_inputs)
    output:
        exposures = outdir + "/{genome}/signature_refit_exposures.csv",
        fractions = outdir + "/{genome}/signature_refit_exposure_fractions.csv",
        families = outdir + "/{genome}/signature_refit_family_exposures.csv",
        summary = outdir + "/{genome}/signature_refit_summary.csv",
        stacked_bar = outdir + "/{genome}/signature_refit_exposure_stacked_bar.png",
        heatmap = outdir + "/{genome}/signature_refit_exposure_heatmap.png"
    log:
        logdir + "/{genome}/spectrum/somatic_signature_refit.log"
    threads: rule_res("somatic_signature_refit").get("threads", 1)
    resources:
        mem_mb = rule_res("somatic_signature_refit").get("mem_mb", 8000),
    params:
        refit_script = os.path.join(ROOT_DIR, "modules/spectrum/bin/signature_refit.py"),
        outprefix = lambda wildcards: outdir + f"/{wildcards.genome}/signature_refit",
        cfg = _spectrum_cfg()
    conda:
        "../spectrum.yaml"
    container:
        sif("../spectrum.yaml", "somatic_signature_refit")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("somatic_signature_refit", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start signature refit for genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.dirname(str(output.exposures))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"somatic_signature_refit_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "somatic_signature_refit.sh")
            sig = params.cfg["signature"]
            cmd = ["python", params.refit_script,
                   "--counts", str(input.counts),
                   "--signatures", str(input.signatures),
                   "--output-prefix", params.outprefix,
                   "--min-mutations", str(sig["min_mutations"]),
                   "--min-exposure-frac", str(sig["min_exposure_frac"]),
                   ]
            if sig["family_map"]:
                cmd += ["--family-map", sig["family_map"]]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write('echo "somatic signature refit completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"Signature refit completed for genome {wildcards.genome}")
        except Exception as e:
            rule_logger.error(f"Error occurred during signature refit for genome {wildcards.genome}: {e}\n")
            raise RuntimeError(f"Error occurred during signature refit for genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.\n")


rule somatic_signature_nmf:
    input:
        unpack(_signature_inputs)
    output:
        signatures = outdir + "/{genome}/signature_nmf_signatures.csv",
        exposures = outdir + "/{genome}/signature_nmf_exposures.csv",
        rank_diagnostics = outdir + "/{genome}/signature_nmf_rank_diagnostics.csv",
        rank_plot = outdir + "/{genome}/signature_nmf_rank_diagnostics.png",
        profiles = outdir + "/{genome}/signature_nmf_signatures.png",
        stacked_bar = outdir + "/{genome}/signature_nmf_exposures_stacked_bar.png",
        cosmic_similarity = outdir + "/{genome}/signature_nmf_cosmic_similarity.csv",
        cosmic_best = outdir + "/{genome}/signature_nmf_cosmic_best_match.csv"
    log:
        logdir + "/{genome}/spectrum/somatic_signature_nmf.log"
    threads: rule_res("somatic_signature_nmf").get("threads", 1)
    resources:
        mem_mb = rule_res("somatic_signature_nmf").get("mem_mb", 8000),
    params:
        nmf_script = os.path.join(ROOT_DIR, "modules/spectrum/bin/signature_nmf.py"),
        outprefix = lambda wildcards: outdir + f"/{wildcards.genome}/signature_nmf",
        cfg = _spectrum_cfg()
    conda:
        "../spectrum.yaml"
    container:
        sif("../spectrum.yaml", "somatic_signature_nmf")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("somatic_signature_nmf", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start signature NMF for genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.dirname(str(output.signatures))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"somatic_signature_nmf_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "somatic_signature_nmf.sh")
            sig = params.cfg["signature"]
            cmd = ["python", params.nmf_script,
                   "--counts", str(input.counts),
                   "--signatures", str(input.signatures),
                   "--output-prefix", params.outprefix,
                   "--max-rank", str(sig["max_rank"]),
                   "--restarts", str(sig["restarts"]),
                   "--seed", str(sig["seed"]),
                   ]
            if sig["rank"]:
                cmd += ["--rank", str(sig["rank"])]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write('echo "somatic signature NMF completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"Signature NMF completed for genome {wildcards.genome}")
        except Exception as e:
            rule_logger.error(f"Error occurred during signature NMF for genome {wildcards.genome}: {e}\n")
            raise RuntimeError(f"Error occurred during signature NMF for genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.\n")


rule spectrum_result:
    input:
        spectrum_png = expand(outdir + "/{genome}/somatic_spectrum_stacked_bar.png", genome=genomes),
        spectrum_heatmap = expand(outdir + "/{genome}/somatic_spectrum_heatmap.png", genome=genomes),
        spectrum_csv = expand(outdir + "/{genome}/somatic_spectrum.csv", genome=genomes),
        per_sample = expand(outdir + "/{genome}/somatic_spectrum_per_sample.csv", genome=genomes),
        sample_stats = expand(outdir + "/{genome}/somatic_spectrum_sample_stats.csv", genome=genomes),
        refit_exposures = expand(outdir + "/{genome}/signature_refit_exposures.csv", genome=genomes),
        refit_summary = expand(outdir + "/{genome}/signature_refit_summary.csv", genome=genomes),
        nmf_signatures = expand(outdir + "/{genome}/signature_nmf_signatures.csv", genome=genomes),
        nmf_rank_diagnostics = expand(outdir + "/{genome}/signature_nmf_rank_diagnostics.csv", genome=genomes)
