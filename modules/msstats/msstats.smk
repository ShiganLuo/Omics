include: "../common/common.smk"


indir = config.get("indir", "data/quantification")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
quantification_method = config.get("quantification_method", "lfq")

# Get parameters
reference_condition = config.get("Params", {}).get("msstats", {}).get("reference_condition", "")
reference_samples = config.get("Params", {}).get("msstats", {}).get("reference_samples", [])

# Get executables
msstats = config.get("Procedure", {}).get("msstats") or "MSstatsConverter"

rule msstats:
    input:
        mztab = indir + "/{quant_method}_quantification.mzTab".format(quant_method=quantification_method)
    output:
        csv = outdir + "/msstats_results.csv"
    log:
        logdir + "/msstats/msstats.log"
    conda:
        "msstats.yaml"
    container:
        sif("msstats.yaml")
    params:
        msstats = msstats,
        reference_condition = reference_condition,
        reference_samples = reference_samples
    run:
        log_path = str(log)
        try:
            open(log_path,"w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger = setup_logger("msstats",log_file=log_path)
            rule_logger.info(f"Start MSstats analysis at {current_time}")
            script = os.path.join(outdir, f"msstats_{current_time}.sh")
            cmd = [
                params.msstats,
                "-in", input.mztab,
                "-out", output.csv
            ]
            if params.reference_condition:
                cmd.extend(["-reference_condition", params.reference_condition])
            if params.reference_samples:
                cmd.extend(["-reference_samples", ",".join(params.reference_samples)])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'MSstats analysis completed at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}'\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during MSstats analysis: {str(e)}\n")
            raise RuntimeError(f"Error during MSstats analysis: {str(e)}\n")

rule msstats_result:
    input:
        outdir + "/msstats_results.csv"
