include: "../../common/common.smk"
msstats_converter = config.get("Procedure", {}).get("msstats") or "MSstatsConverter"

indir = config.get("indir", "data/quantification")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
quantification_method = (config.get("Params", {}).get("workflow", {}).get("quantification_method") or "lfq").lower()

expdesign = config.get("expdesign", "")

_MSSTATS_SOURCE = {
    "tmt": indir + "/tmt_quantification.consensusXML",
    "lfq": indir + "/lfq_msstats.csv",
    "dia": indir + "/dia_msstats.csv",
}

rule msstats:
    input:
        **({"source": _MSSTATS_SOURCE.get(quantification_method, indir + "/lfq_msstats.csv")}
           if (expdesign or quantification_method != "tmt") else {})
    output:
        csv = outdir + "/msstats_input.csv"
    log:
        logdir + "/msstats/msstats.log"
    threads: rule_res("msstats").get("threads", 1)
    benchmark: logdir + "/msstats/msstats.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "msstats", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml", "msstats")
    params:
        msstats_converter = msstats_converter,
        expdesign = expdesign,
        method = quantification_method,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("msstats", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start MSstats input preparation at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            plain_script = os.path.join(outdir, f"msstats.sh")
            if params.method == "tmt":
                if not params.expdesign:
                    rule_logger.info(
                        "Params.workflow.expdesign is empty: skipping TMT "
                        "MSstats input preparation (needs a plex->channel->"
                        "condition design). Protein-level deliverable is "
                        "7_quantification/tmt_protein_abundances.csv."
                    )
                    os.makedirs(str(outdir), exist_ok=True)
                    open(str(output.csv), "w").close()
                else:
                    cmd = [
                        params.msstats_converter,
                        "-in", str(input.source),
                        "-in_design", str(params.expdesign),
                        "-method", "ISO",
                        "-out", output.csv,
                    ]
                    script = receipt_script(f"msstats_{current_time}.sh")
                    with open(script, "w") as f:
                        f.write("#!/bin/bash\nset -euo pipefail\n")
                        f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                        f.write(f"echo 'MSstats input preparation completed at {current_time}'\n")
                    shutil.copyfile(script, plain_script)
                    if SCRIPT_ONLY:
                        script_only_touch(output)
                    else:
                        shell(f"bash {plain_script} >> {log_path} 2>&1")
            else:
                rule_logger.info(f"Staging {input.source} as {output.csv}")
                shutil.copy(str(input.source), str(output.csv))
                rule_logger.info(f"MSstats input preparation completed at {current_time}")
        except Exception as e:
            rule_logger.error(f"Error during MSstats input preparation: {str(e)}\n")
            raise RuntimeError(f"Error during MSstats input preparation: {str(e)}\n")

rule msstats_result:
    input:
        outdir + "/msstats_input.csv"
