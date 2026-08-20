include: "../../common/common.smk"

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
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("protein_inference", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start protein inference for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/protein_inference_{current_time}.sh")
            cmd = [
                params.epifany,
                "-in", input.idxml,
                "-out", output.protein_idxml,
                "-algorithm:method", params.method,
                "-algorithm:greedy_group_resolution", params.greedy_group_resolution
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"rule protein_inference was call failed,error: {e}")
            raise RuntimeError(f"rule protein_inference was call failed,error: {e}")

rule protein_inference_result:
    input:
        protein_idxml = outdir + "/{sample_id}/{sample_id}_protein.idXML"
