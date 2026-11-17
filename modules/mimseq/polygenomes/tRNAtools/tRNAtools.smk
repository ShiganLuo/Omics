include: "../../../common/common.smk"


ROOT_DIR = config.get("ROOT_DIR", ".")

# 参数精确到 rule:Params.mimseq.mimseq_tRNAtools.<key>(rule_param 见 common.smk;
# 无 module 级共享回退;data_dir 归 mimseq_tRNAtools 所有)。
_mimseq_params = config.get("Params", {}).get("mimseq", {})
_data_dir_mimseq_tRNAtools = rule_param(_mimseq_params, "mimseq_tRNAtools", "data_dir", "") or os.path.join(
    ROOT_DIR, "modules", "mimseq", "mimseq", "data")
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
    species = rule_param(_mimseq_params, "mimseq_tRNAtools", "species", "")
    if species:
        # Using built-in species, no external files needed
        return {}
    else:
        # Using custom tRNA files
        return {
            "trnas": config["genome"][wildcards.genome].get("trnas", ""),
            "trnaout": config["genome"][wildcards.genome].get("trnaout", ""),
        }

rule mimseq_tRNAtools:
    """Parse tRNA sequences, generate SNP index and GSNAP indices."""
    input:
        unpack(get_tRNAtools_input)
    output:
        state_dir = directory(outdir + "/{genome}/state"),
    log:
        logdir + "/{genome}/mimseq_tRNAtools.log"
    params:
        script = os.path.join(ROOT_DIR, "modules", "mimseq", "bin", "tRNAtools", "run.py"),
        mimseq = config.get("Procedure", {}).get("mimseq") or "mimseq",
        species = rule_param(_mimseq_params, "mimseq_tRNAtools", "species", ""),
        name = rule_param(_mimseq_params, "mimseq_tRNAtools", "name", "tRNAseq"),
        cluster_id = rule_param(_mimseq_params, "mimseq_tRNAtools", "cluster_id", 0.97),
        double_cca = rule_param(_mimseq_params, "mimseq_tRNAtools", "double_cca", False),
        posttrans = rule_param(_mimseq_params, "mimseq_tRNAtools", "posttrans_mod_off", False),
        pretrnas = rule_param(_mimseq_params, "mimseq_tRNAtools", "pretRNAs", False),
        local_mod = rule_param(_mimseq_params, "mimseq_tRNAtools", "local_modomics", False),
        mito_trnas = lambda wildcards: config["genome"][wildcards.genome].get("mito_trnas", ""),
        plastid_trnas = lambda wildcards: config["genome"][wildcards.genome].get("plastid_trnas", ""),
        modifications = os.path.join(ROOT_DIR, "modules", "mimseq", "mimseq", "modifications"),
        data_dir = _data_dir_mimseq_tRNAtools,
    threads: rule_res("mimseq_tRNAtools").get("threads", 4)
    benchmark: logdir + "/{genome}/mimseq_tRNAtools.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_tRNAtools", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../../mimseq.yaml"
    container:
        sif("../../mimseq.yaml", "mimseq_tRNAtools")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(RECEIPT_DIR, f"mimseq_tRNAtools_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(outdir, wildcards.genome, f"mimseq_tRNAtools.sh")

            cmd = [
                "python", params.script,
                "--modifications", params.modifications,
                "--name", params.name,
                "--out", outdir + f"/{wildcards.genome}/",
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
            if rule_param(_mimseq_params, "mimseq_tRNAtools", "snp_tolerance", True):
                cmd += ["--snp-tolerance"]
            if rule_param(_mimseq_params, "mimseq_tRNAtools", "cluster", True):
                cmd += ["--cluster"]

            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule mimseq_tRNAtools_result:
    """Result rule for tRNAtools module."""
    input:
        state_dir = outdir + "/{genome}/state",

    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "mimseq_tRNAtools_result", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
