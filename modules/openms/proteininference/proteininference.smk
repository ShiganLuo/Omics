include: "../../common/common.smk"

import shutil

indir = config.get("indir", "data/psm_fdr")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

# Get parameters
method = config.get("Params", {}).get("protein_inference", {}).get("method", "epifany")
greedy_group_resolution = config.get("Params", {}).get("protein_inference", {}).get("greedy_group_resolution", "none")

# Get executables
epifany = config.get("Procedure", {}).get("epifany") or "Epifany"

rule protein_inference:
    input:
        idxml = indir + "/{sample_id}/{sample_id}_filtered.idXML"
    output:
        protein_idxml = outdir + "/{sample_id}/{sample_id}_protein.idXML"
    log:
        logdir + "/{sample_id}/protein_inference.log"
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        epifany = epifany,
        method = method,
        greedy_group_resolution = greedy_group_resolution
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("protein_inference", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start protein inference for sample {wildcards.sample_id} at {current_time}")
            # Epifany segfaults on idXML without peptide hits (e.g. when no PSM
            # survived FDR); fail early with an actionable message instead.
            n_hits = 0
            with open(input.idxml) as fh:
                for line in fh:
                    if "<PeptideHit" in line:
                        n_hits += 1
                        break
            if n_hits == 0:
                raise RuntimeError(
                    f"no peptide hits in {input.idxml}: zero PSMs survived FDR "
                    f"for sample {wildcards.sample_id}; check the psm_fdr step "
                    "(search/identification quality) before protein inference"
                )
            os.makedirs(os.path.join(logdir, "cmd", str(wildcards.sample_id)), exist_ok=True)
            script = os.path.join(logdir, "cmd", str(wildcards.sample_id), f"{wildcards.sample_id}/protein_inference_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/protein_inference.sh")
            # Epifany (OpenMS 3.1) has no -algorithm:method: it only runs the
            # Bayesian inference; greedy_group_resolution is a TOPP-level
            # option (-greedy_group_resolution), not an algorithm subsection one.
            if params.method not in ("", "epifany"):
                raise ValueError(
                    f"Params.protein_inference.method={params.method!r} is not supported: "
                    "this Epifany build only provides Bayesian 'epifany' inference"
                )
            cmd = [
                params.epifany,
                "-in", input.idxml,
                "-out", output.protein_idxml,
                "-greedy_group_resolution", params.greedy_group_resolution,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'Protein inference for {wildcards.sample_id} was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule protein_inference was call failed,error: {e}")
            raise RuntimeError(f"rule protein_inference was call failed,error: {e}")

rule protein_inference_result:
    input:
        protein_idxml = outdir + "/{sample_id}/{sample_id}_protein.idXML"
