from snakemake.logging import logger
import time
import os
import sys

ROOT_DIR = config.get("ROOT_DIR", ".")
data_dir = config.get("Params", {}).get("mimseq", {}).get("data_dir", "") or os.path.join(ROOT_DIR, "modules", "mimseq", "mimseq", "data")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

# Import serialization utilities
_bin_dir = os.path.join(ROOT_DIR, "modules", "mimseq", "bin")
if _bin_dir not in sys.path:
    sys.path.insert(0, _bin_dir)
from serialize import load_tRNA_state

def get_tRNAtools_input(wildcards):
    """Return input files for tRNAtools based on species configuration."""
    species = config.get("Params", {}).get("mimseq", {}).get("species", "")
    if species:
        # Using built-in species, no external files needed
        return {}
    else:
        # Using custom tRNA files
        return {
            "trnas": config.get("genome", {}).get("trnas", ""),
            "trnaout": config.get("genome", {}).get("trnaout", ""),
        }

rule mimseq_tRNAtools:
    """Parse tRNA sequences, generate SNP index and GSNAP indices."""
    input:
        unpack(get_tRNAtools_input)
    output:
        state_dir = directory(outdir + "/state"),
    log:
        logdir + "/mimseq_tRNAtools.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "tRNAtools", "run.py"),
        mimseq = config.get("Procedure", {}).get("mimseq") or "mimseq",
        species = config.get("Params", {}).get("mimseq", {}).get("species", ""),
        name = config.get("Params", {}).get("mimseq", {}).get("name", "tRNAseq"),
        cluster_id = config.get("Params", {}).get("mimseq", {}).get("cluster_id", 0.97),
        double_cca = config.get("Params", {}).get("mimseq", {}).get("double_cca", False),
        posttrans = config.get("Params", {}).get("mimseq", {}).get("posttrans_mod_off", False),
        pretrnas = config.get("Params", {}).get("mimseq", {}).get("pretRNAs", False),
        local_mod = config.get("Params", {}).get("mimseq", {}).get("local_modomics", False),
        mito_trnas = config.get("genome", {}).get("mito_trnas", ""),
        plastid_trnas = config.get("genome", {}).get("plastid_trnas", ""),
        modifications = os.path.join(ROOT_DIR, "modules", "mimseq", "mimseq", "modifications"),
        data_dir = data_dir,
    threads: 4
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"mimseq_tRNAtools_{current_time}.sh")

            cmd = [
                "/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_/bin/python", params.script,
                "--modifications", params.modifications,
                "--name", params.name,
                "--out", outdir + "/",
                "--species", params.species,
                "--data-dir", params.data_dir,
                "--threads", str(threads),
                "--cluster-id", str(params.cluster_id),
            ]
            # Add trnas and trnaout if using custom files
            if hasattr(input, "trnas") and input.trnas:
                cmd += ["--trnas", input.trnas]
            if hasattr(input, "trnaout") and input.trnaout:
                cmd += ["--trnaout", input.trnaout]
            if params.mito_trnas:
                cmd += ["--mito-trnas", params.mito_trnas]
            if params.plastid_trnas:
                cmd += ["--plastid-trnas", params.plastid_trnas]
            if params.double_cca:
                cmd += ["--double-cca"]
            if params.posttrans:
                cmd += ["--posttrans"]
            if params.pretrnas:
                cmd += ["--pretrnas"]
            if params.local_mod:
                cmd += ["--local-mod"]
            if config.get("Params", {}).get("mimseq", {}).get("snp_tolerance", True):
                cmd += ["--snp-tolerance"]
            if config.get("Params", {}).get("mimseq", {}).get("cluster", True):
                cmd += ["--cluster"]

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("source /home/zhousha/miniforge3/bin/activate /data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"mimseq_tRNAtools failed: {e}\n")
            raise

rule mimseq_tRNAtools_result:
    """Result rule for tRNAtools module."""
    input:
        state_dir = outdir + "/state",
