include: "../common/common.smk"

from snakemake.logging import logger
import os
import time

outdir = config.get("outdir", "results/ldsc")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})

ldsc_cfg = config.get("Params", {}).get("ldsc", {})
gwas = ldsc_cfg.get("gwas", {}) or {}
traits = sorted(gwas.keys())
BIN = ROOT_DIR + "/modules/ldsc/bin"

wildcard_constraints:
    trait = ".+",


rule ldsc_annot:
    """Split FIRE peak clusters into LDSC annotation BEDs."""
    input:
        peaks = lambda wildcards: _require(ldsc_cfg.get("annot_bed"),
                                           "Params.ldsc.annot_bed", "LDSC"),
    output:
        manifest = outdir + "/annot/" + (ldsc_cfg.get("annot_prefix") or "peaks") + ".annot_manifest.tsv",
    log:
        logdir + "/ldsc_annot.log"
    threads: 1
    conda:
        "ldsc.yaml"
    container:
        sif("ldsc.yaml")
    params:
        script = BIN + "/ldsc_annot.py",
        prefix = ldsc_cfg.get("annot_prefix") or "peaks",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("ldsc_annot", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ldsc_annot at {current_time}")
            group_outdir = outdir + "/annot"
            os.makedirs(group_outdir, exist_ok=True)
            script = os.path.join(group_outdir, f"ldsc_annot_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-p", str(input.peaks),
                "-o", group_outdir,
                "--prefix", params.prefix,
                "-l", log_path,
            ]
            if ldsc_cfg.get("annot_cluster_from_name"):
                cmd.append("-c")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "ldsc_annot at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during ldsc_annot: {e}\n")
            raise RuntimeError(f"Error occurred during ldsc_annot: {e}\n")


rule ldsc_munge:
    """Munge one GWAS sumstats file with munge_sumstats.py."""
    input:
        sumstats = lambda wildcards: _require(
            gwas.get(wildcards.trait), f"Params.ldsc.gwas.{wildcards.trait}", "LDSC"),
    output:
        munged = outdir + "/{trait}/{trait}.sumstats.gz",
    log:
        logdir + "/ldsc_munge_{trait}.log"
    threads: 1
    conda:
        "ldsc.yaml"
    container:
        sif("ldsc.yaml")
    params:
        munge = config.get("Procedure", {}).get("munge_sumstats") or "munge_sumstats.py",
        extra = ldsc_cfg.get("munge_args", ""),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("ldsc_munge", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ldsc_munge for {wildcards.trait} at {current_time}")
            trait_outdir = outdir + f"/{wildcards.trait}"
            os.makedirs(trait_outdir, exist_ok=True)
            script = os.path.join(trait_outdir, f"ldsc_munge_{current_time}.sh")
            cmd = [
                params.munge,
                "--sumstats", str(input.sumstats),
                "--out", trait_outdir + f"/{wildcards.trait}",
                "--merge-alleles", _require(ldsc_cfg.get("merge_alleles"),
                                            "Params.ldsc.merge_alleles", "LDSC"),
            ]
            if params.extra:
                cmd += str(params.extra).split()
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "ldsc_munge for {wildcards.trait} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during ldsc_munge for {wildcards.trait}: {e}\n")
            raise RuntimeError(f"Error occurred during ldsc_munge for {wildcards.trait}: {e}\n")


rule ldsc_h2_cts:
    """Partitioned heritability (ldsc.py --h2-cts) for one GWAS trait."""
    input:
        munged = outdir + "/{trait}/{trait}.sumstats.gz",
        manifest = outdir + "/annot/" + (ldsc_cfg.get("annot_prefix") or "peaks") + ".annot_manifest.tsv",
    output:
        results = outdir + "/{trait}/{trait}.cts.results",
    log:
        logdir + "/ldsc_h2_cts_{trait}.log"
    threads: 1
    conda:
        "ldsc.yaml"
    container:
        sif("ldsc.yaml")
    params:
        ldsc = config.get("Procedure", {}).get("ldsc") or "ldsc.py",
        ref_ld = ldsc_cfg.get("ref_ld_chr"),
        w_ld = ldsc_cfg.get("w_ld_chr"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("ldsc_h2_cts", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ldsc_h2_cts for {wildcards.trait} at {current_time}")
            trait_outdir = outdir + f"/{wildcards.trait}"
            os.makedirs(trait_outdir, exist_ok=True)
            script = os.path.join(trait_outdir, f"ldsc_h2_cts_{current_time}.sh")
            cmd = [
                params.ldsc,
                "--h2-cts", str(input.munged),
                "--ref-ld-chr", str(_require(params.ref_ld, "Params.ldsc.ref_ld_chr", "LDSC")),
                "--w-ld-chr", str(_require(params.w_ld, "Params.ldsc.w_ld_chr", "LDSC")),
                "--out", trait_outdir + f"/{wildcards.trait}",
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "ldsc_h2_cts for {wildcards.trait} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during ldsc_h2_cts for {wildcards.trait}: {e}\n")
            raise RuntimeError(f"Error occurred during ldsc_h2_cts for {wildcards.trait}: {e}\n")


rule ldsc_result:
    """Panel of per-trait LDSC results."""
    input:
        results = expand(outdir + "/{trait}/{trait}.cts.results", trait=traits),
