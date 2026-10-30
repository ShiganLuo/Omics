include: "../../common/common.smk"

import shutil


indir = config.get("indir", "data/quantification")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
quantification_method = (config.get("quantification_method") or "lfq").lower()

# Experimental design TSV (OpenMS design format); wired from
# Params.workflow.expdesign at subworkflow level. MSstatsConverter REQUIRES
# it for the isobaric path, ProteomicsLFQ uses it to fill Run/Condition.
expdesign = config.get("expdesign", "")

# Get executables
msstats_converter = config.get("Procedure", {}).get("msstats") or "MSstatsConverter"

# Per-method product of the quantification module that carries the peptide
# intensities for the MSstats input table:
#   lfq/dia -> ProteomicsLFQ -out_msstats CSV (already MSstats input format)
#   tmt     -> merged isobaric consensusXML (only produced with an expdesign),
#              converted here by MSstatsConverter
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
    threads: 1
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        msstats_converter = msstats_converter,
        expdesign = expdesign,
        method = quantification_method,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("msstats", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start MSstats input preparation at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"msstats_{current_time}.sh")
            plain_script = os.path.join(outdir, f"msstats.sh")
            if params.method == "tmt":
                # Isobaric consensusXML -> MSstats (MSstatsTMT) input CSV.
                # Requires an experimental design AND a merged consensusXML
                # with intact per-plex sample identity; per the quantification
                # module, merged-design mzTab/MSstats inputs are only produced
                # when Params.workflow.expdesign is set. Without a design the
                # deliverable is the long abundance table instead, so skip.
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
                    with open(script, "w") as f:
                        f.write("#!/bin/bash\nset -euo pipefail\n")
                        f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                        f.write(f"echo 'MSstats input preparation completed at {current_time}'\n")
                    shutil.copyfile(script, plain_script)
                    shell(f"bash {plain_script} >> {log_path} 2>&1")
            else:
                # ProteomicsLFQ -out_msstats already writes the MSstats input
                # CSV; stage it under the uniform final path.
                rule_logger.info(f"Staging {input.source} as {output.csv}")
                shutil.copy(str(input.source), str(output.csv))
                rule_logger.info(f"MSstats input preparation completed at {current_time}")
        except Exception as e:
            rule_logger.error(f"Error during MSstats input preparation: {str(e)}\n")
            raise RuntimeError(f"Error during MSstats input preparation: {str(e)}\n")


rule msstats_result:
    input:
        outdir + "/msstats_input.csv"
