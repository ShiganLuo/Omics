from snakemake.logging import logger
include: "../common/common.smk"
import time
import os

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])

# ========================
# mimseq main run
# ========================

rule mimseq_run:
    """Run mim-tRNAseq pipeline: alignment, coverage, modification quantification, DESeq2."""
    input:
        sample_data = config.get("sample_data", ""),
        trnas = config.get("genome", {}).get("trnas", ""),
        trnaout = config.get("genome", {}).get("trnaout", ""),
    output:
        outdir = directory(outdir + "/mimseq"),
    log:
        logdir + "/mimseq/mimseq_run.log"
    conda:
        "mimseq.yaml"
    params:
        script = "mimseq",
        species = config.get("Params", {}).get("mimseq", {}).get("species", ""),
        name = config.get("Params", {}).get("mimseq", {}).get("name", "tRNAseq"),
        control_cond = config.get("Params", {}).get("mimseq", {}).get("control_cond", ""),
        cluster_id = config.get("Params", {}).get("mimseq", {}).get("cluster_id", 0.97),
        threads = config.get("Params", {}).get("mimseq", {}).get("threads", 8),
        min_cov = config.get("Params", {}).get("mimseq", {}).get("min_cov", 0.0005),
        max_mismatches = config.get("Params", {}).get("mimseq", {}).get("max_mismatches", 0.075),
        max_multi = config.get("Params", {}).get("mimseq", {}).get("max_multi", 4),
        misinc_thresh = config.get("Params", {}).get("mimseq", {}).get("misinc_thresh", 0.1),
        remap_mismatches = config.get("Params", {}).get("mimseq", {}).get("remap_mismatches", 0.05),
        p_adj = config.get("Params", {}).get("mimseq", {}).get("p_adj", 0.05),
        mito_trnas = config.get("genome", {}).get("mito_trnas", ""),
        plastid_trnas = config.get("genome", {}).get("plastid_trnas", ""),
        no_cluster = config.get("Params", {}).get("mimseq", {}).get("no_cluster", False),
        no_cca = config.get("Params", {}).get("mimseq", {}).get("no_cca", False),
        double_cca = config.get("Params", {}).get("mimseq", {}).get("double_cca", False),
        remap = config.get("Params", {}).get("mimseq", {}).get("remap", True),
        snp_tolerance = config.get("Params", {}).get("mimseq", {}).get("snp_tolerance", True),
        keep_temp = config.get("Params", {}).get("mimseq", {}).get("keep_temp", False),
        crosstalks = config.get("Params", {}).get("mimseq", {}).get("crosstalks", False),
        pretRNAs = config.get("Params", {}).get("mimseq", {}).get("pretRNAs", False),
        posttrans = config.get("Params", {}).get("mimseq", {}).get("posttrans_mod_off", False),
        cov_diff = config.get("Params", {}).get("mimseq", {}).get("cov_diff", 0.5),
        local_mod = config.get("Params", {}).get("mimseq", {}).get("local_modomics", False),
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            mimseq_outdir = f"{outdir}/mimseq"
            script_path = os.path.join(outdir, f"mimseq_run_{current_time}.sh")

            # Build command
            cmd = [params.script]

            # Species or custom tRNA references
            if params.species:
                cmd += ["-s", params.species]
            else:
                cmd += ["-t", input.trnas, "-o", input.trnaout]
                if params.mito_trnas:
                    cmd += ["-m", params.mito_trnas]
                if params.plastid_trnas:
                    cmd += ["-p", params.plastid_trnas]

            # Required options
            cmd += [
                "-n", params.name,
                "--out-dir", mimseq_outdir,
                "--control-condition", params.control_cond,
                "--threads", str(params.threads),
                "--cluster-id", str(params.cluster_id),
                "--min-cov", str(params.min_cov),
                "--max-mismatches", str(params.max_mismatches),
                "--max-multi", str(params.max_multi),
                "--misinc-thresh", str(params.misinc_thresh),
                "--remap-mismatches", str(params.remap_mismatches),
                "--p-adj", str(params.p_adj),
                "--deconv-cov-ratio", str(params.cov_diff),
            ]

            # Boolean flags
            if params.no_cluster:
                cmd += ["--no-cluster"]
            if params.no_cca:
                cmd += ["--no-cca-analysis"]
            if params.double_cca:
                cmd += ["--double-cca"]
            if params.remap:
                cmd += ["--remap"]
            if not params.snp_tolerance:
                cmd += ["--no-snp-tolerance"]
            if params.keep_temp:
                cmd += ["--keep-temp"]
            if params.crosstalks:
                cmd += ["--crosstalks"]
            if params.pretRNAs:
                cmd += ["--pretRNAs"]
            if params.posttrans:
                cmd += ["--posttrans-mod-off"]
            if params.local_mod:
                cmd += ["--local-modomics"]

            # Sample data (positional, must be last)
            cmd += [input.sample_data]

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"mimseq_run failed: {e}\n")
            raise

# ========================
# Result aggregation
# ========================

rule mimseq_result:
    """Aggregate mimseq outputs as dependency endpoint."""
    input:
        outdir = outdir + "/mimseq",
