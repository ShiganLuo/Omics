shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
import time

indir = config.get("indir", "data/protein_inference")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
quantification_method = config.get("quantification_method", "lfq")

# Get parameters
protein_quant_params = config.get("Params", {}).get("protein_quant", {})
tmt_params = config.get("Params", {}).get("tmt", {})
lfq_params = config.get("Params", {}).get("lfq", {})
dia_params = config.get("Params", {}).get("dia", {})

# Get executables
proteomicslfq = config.get("Procedure", {}).get("proteomicslfq") or "ProteomicsLFQ"
proteinquantifier = config.get("Procedure", {}).get("proteinquantifier") or "ProteinQuantifier"

# Get experimental design file
expdesign = config.get("expdesign", "")

rule quantification_lfq:
    input:
        idxml = expand(indir + "/{sid}/{sid}_protein.idXML", sid=samples)
    output:
        mztab = outdir + "/lfq_quantification.mzTab"
    log:
        logdir + "/quantification_lfq.log"
    threads: 8
    conda:
        "openms.yaml"
    params:
        proteomicslfq = proteomicslfq,
        expdesign = expdesign,
        protein_inference_method = protein_quant_params.get("method", "feature_intensity"),
        top = protein_quant_params.get("top", 3),
        average = protein_quant_params.get("average", "median"),
        best_charge_and_score = protein_quant_params.get("best_charge_and_score", True)
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start LFQ quantification at {current_time}")
        script = os.path.join(outdir, f"quantification_lfq_{current_time}.sh")
        cmd = [
            params.proteomicslfq,
            "-in", " ".join(input.idxml),
            "-out", output.mztab,
            "-threads", str(threads),
            "-protein_inference", params.protein_inference_method,
            "-quantification_method", "feature_intensity",
            "-top", str(params.top),
            "-average", params.average
        ]
        if params.best_charge_and_score:
            cmd.append("-best_charge_and_score")
        if params.expdesign:
            cmd.extend(["-design", params.expdesign])
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule quantification_tmt:
    input:
        idxml = expand(indir + "/{sid}/{sid}_protein.idXML", sid=samples)
    output:
        mztab = outdir + "/tmt_quantification.mzTab"
    log:
        logdir + "/quantification_tmt.log"
    threads: 8
    conda:
        "openms.yaml"
    params:
        proteinquantifier = proteinquantifier,
        channel_mass_tolerance = tmt_params.get("channel_mass_tolerance", 0.003),
        channel_annotation = tmt_params.get("channel_annotation", ""),
        isotope_correction = tmt_params.get("isotope_correction", True)
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start TMT quantification at {current_time}")
        script = os.path.join(outdir, f"quantification_tmt_{current_time}.sh")
        cmd = [
            params.proteinquantifier,
            "-in", " ".join(input.idxml),
            "-out", output.mztab,
            "-threads", str(threads),
            "-channel_mass_tolerance", str(params.channel_mass_tolerance)
        ]
        if params.channel_annotation:
            cmd.extend(["-channel_annotation", params.channel_annotation])
        if params.isotope_correction:
            cmd.append("-isotope_correction")
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule quantification_dia:
    input:
        idxml = expand(indir + "/{sid}/{sid}_protein.idXML", sid=samples)
    output:
        mztab = outdir + "/dia_quantification.mzTab"
    log:
        logdir + "/quantification_dia.log"
    threads: 8
    conda:
        "openms.yaml"
    params:
        proteomicslfq = proteomicslfq,
        library = dia_params.get("library", ""),
        dia_window = dia_params.get("dia_window", "")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start DIA quantification at {current_time}")
        script = os.path.join(outdir, f"quantification_dia_{current_time}.sh")
        cmd = [
            params.proteomicslfq,
            "-in", " ".join(input.idxml),
            "-out", output.mztab,
            "-threads", str(threads),
            "-quantification_method", "feature_intensity"
        ]
        if params.library:
            cmd.extend(["-library", params.library])
        if params.dia_window:
            cmd.extend(["-dia_window", params.dia_window])
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

# Select quantification method
if quantification_method == "tmt":
    rule quantification_result:
        input:
            outdir + "/tmt_quantification.mzTab"
elif quantification_method == "lfq":
    rule quantification_result:
        input:
            outdir + "/lfq_quantification.mzTab"
elif quantification_method == "dia":
    rule quantification_result:
        input:
            outdir + "/dia_quantification.mzTab"
else:
    rule quantification_result:
        input:
            outdir + "/lfq_quantification.mzTab"
