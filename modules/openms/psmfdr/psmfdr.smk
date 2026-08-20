include: "../../common/common.smk"

indir = config.get("indir", "data/psm_rescoring")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

# Get parameters
fdr = config.get("Params", {}).get("psm_fdr_control", {}).get("fdr", 0.01)
method = config.get("Params", {}).get("psm_fdr_control", {}).get("method", "BH")

# Get executables
openms = config.get("Procedure", {}).get("openms") or "FalseDiscoveryRate"

rule psm_fdr:
    input:
        idxml = indir + "/{sample_id}/{sample_id}_scored.idXML"
    output:
        filtered_idxml = outdir + "/{sample_id}/{sample_id}_filtered.idXML"
    log:
        logdir + "/{sample_id}/psm_fdr.log"
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        openms = openms,
        fdr = fdr,
        method = method
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("psm_fdr", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start PSM FDR control for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/psm_fdr_{current_time}.sh")
            cmd = [
                params.openms,
                "-in", input.idxml,
                "-out", output.filtered_idxml,
                "-FDR:PSM", str(params.fdr),
                "-method", params.method
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"rule psm_fdr was call failed,error: {e}")
            raise RuntimeError(f"rule psm_fdr was call failed,error: {e}")
        

rule psm_fdr_result:
    input:
        filtered_idxml = outdir + "/{sample_id}/{sample_id}_filtered.idXML"
