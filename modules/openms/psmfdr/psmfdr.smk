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
        open(log_path, 'w').close()
        rule_logger = setup_logger("psm_fdr", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start PSM FDR control for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"psm_fdr_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/psm_fdr.sh")
            cmd = [
                params.openms,
                "-in", input.idxml,
                "-out", output.filtered_idxml,
                "-FDR:PSM", str(params.fdr),
                "-protein", "false",
            ]
            # FalseDiscoveryRate is the single FDR authority of the chain: it
            # recomputes target-decoy q-values from the Percolator scores and
            # overwrites the percolator q-value annotations, so filtering and
            # the reported q-values always come from one estimator.
            # method: "BH" (or "qvalue") = Benjamini-Hochberg q-values
            # (default); "strict" = plain FDR without q-value smoothing.
            if str(params.method).lower() in ("bh", "qvalue"):
                pass
            elif str(params.method).lower() == "strict":
                cmd.append("-algorithm:no_qvalues")
            else:
                raise ValueError(
                    f"Params.psm_fdr_control.method={params.method!r} is not "
                    "supported: use 'BH' (q-values) or 'strict'"
                )
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "PSM FDR control for {wildcards.sample_id} was completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule psm_fdr was call failed,error: {e}")
            raise RuntimeError(f"rule psm_fdr was call failed,error: {e}")
        

rule psm_fdr_result:
    input:
        filtered_idxml = outdir + "/{sample_id}/{sample_id}_filtered.idXML"
