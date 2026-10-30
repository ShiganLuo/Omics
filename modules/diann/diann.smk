include: "../common/common.smk"

import os
import time

indir = config.get("indir", "common/2_raw2mzml")
outdir = config.get("outdir", "common/3_diann")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
p = config.get("Params", {}).get("diann", {})

def _require(value, key):
    if not value:
        raise ValueError(
            f"Missing file for DIA-NN: set '{key}' in the workflow config JSON, "
            f"then rerun."
        )
    return value

def _mzml_inputs(wildcards=None):
    return [f"{indir}/{s}/{s}.mzML" for s in samples]

rule diann_run:
    """DIA-NN library-free (or library) analysis over all DIA runs jointly.

    One job for the whole cohort: DIA-NN aligns runs and controls FDR across
    them, which is the recommended usage. Outputs the main report plus stats.
    """
    input:
        mzmls = _mzml_inputs,
        fasta = lambda wildcards: _require(genome_cfg.get("protein_fasta"),
                                           "genome.protein_fasta"),
    output:
        report = outdir + "/diann.report.tsv",
        stats = outdir + "/diann.report.stats.tsv",
    log:
        logdir + "/group/diann_run.log"
    threads: 16
    conda:
        "diann.yaml"
    container:
        sif("diann.yaml")
    params:
        diann = config.get("Procedure", {}).get("diann") or "diann",
        qvalue = p.get("qvalue", 0.01),
        matrix_qvalue = p.get("matrix_qvalue", 0.01),
        mass_acc_ms1 = p.get("mass_acc_ms1", 20),
        mass_acc = p.get("mass_acc", 20),
        min_fr_mz = p.get("min_fr_mz", 200),
        max_fr_mz = p.get("max_fr_mz", 1800),
        min_pept_len = p.get("min_pept_len", 7),
        max_pept_len = p.get("max_pept_len", 30),
        fixed_mod = p.get("fixed_mod") or "Carbamidomethyl (C),57.0215",
        var_mods = p.get("var_mods") or "Oxidation (M),15.9949",
        library = p.get("library"),
        gen_spec_lib = p.get("gen_spec_lib", True),
        predictor = p.get("predictor", True),
        verbose = p.get("verbose", True),
        extra = p.get("extra") or "",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("diann_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DIA-NN over {len(samples)} runs at {current_time}")
            group_outdir = outdir
            os.makedirs(group_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"diann_run_{current_time}.sh")
            plain_script = os.path.join(group_outdir, f"diann_run.sh")
            cmd = [
                params.diann,
                "--out", str(output.report),
                "--fasta", str(input.fasta),
                "--threads", str(threads),
                "--qvalue", str(params.qvalue),
                "--matrix-qvalue", str(params.matrix_qvalue),
                "--mass-acc-ms1", str(params.mass_acc_ms1),
                "--mass-acc", str(params.mass_acc),
                "--min-fr-mz", str(params.min_fr_mz),
                "--max-fr-mz", str(params.max_fr_mz),
                "--min-pept-len", str(params.min_pept_len),
                "--max-pept-len", str(params.max_pept_len),
                "--fixed-mod", str(params.fixed_mod),
                "--var-mods", str(params.var_mods),
            ]
            if params.library:
                cmd += ["--lib", str(params.library)]
            else:
                if params.gen_spec_lib:
                    cmd.append("--gen-spec-lib")
                if params.predictor:
                    cmd.append("--predictor")
            if params.verbose:
                cmd.append("--verbose")
            if params.extra:
                cmd += str(params.extra).split()
            cmd += [str(m) for m in input.mzmls]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "diann_run at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during diann_run: {e}\n")
            raise RuntimeError(f"Error occurred during diann_run: {e}\n")

rule diann_summary:
    """Per-run QC summary and protein matrix from the DIA-NN report."""
    input:
        report = outdir + "/diann.report.tsv",
    output:
        qc = outdir + "/diann_qc.tsv",
        matrix = outdir + "/diann_protein_matrix.tsv",
    log:
        logdir + "/group/diann_summary.log"
    threads: 1
    conda:
        "diann.yaml"
    container:
        sif("diann.yaml")
    params:
        script = ROOT_DIR + "/modules/diann/bin/diann_summary.py",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("diann_summary", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DIA-NN summary at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"diann_summary_{current_time}.sh")
            plain_script = os.path.join(outdir, f"diann_summary.sh")
            cmd = [
                "python", params.script,
                "-i", str(input.report),
                "-o", str(outdir),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "diann_summary at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during diann_summary: {e}\n")
            raise RuntimeError(f"Error occurred during diann_summary: {e}\n")

rule diann_result:
    input:
        report = outdir + "/diann.report.tsv",
        stats = outdir + "/diann.report.stats.tsv",
        qc = outdir + "/diann_qc.tsv",
        matrix = outdir + "/diann_protein_matrix.tsv",
