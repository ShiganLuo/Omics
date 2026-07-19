include: "../../common/common.smk"

shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
import time

indir = config.get("indir", "data/mzml")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
mzml_files = config.get("mzml_files", [])

# Get parameters
search_engine_params = config.get("Params", {}).get("search_engine", {})
search_engines = config.get("search_engines", "comet").split(",")

# Get executables
comet = config.get("Procedure", {}).get("comet") or "CometAdapter"
msgf = config.get("Procedure", {}).get("msgf") or "MSGFPlusAdapter"
sage = config.get("Procedure", {}).get("sage") or "SageAdapter"

# Get decoy database
decoy_fasta = config.get("genome", {}).get("decoy_fasta")

def get_mzml_file(wildcards):
    """Get mzML file for a sample."""
    for i, sample_id in enumerate(samples):
        if sample_id == wildcards.sample_id:
            return mzml_files[i]
    raise ValueError(f"Sample {wildcards.sample_id} not found in samples list")

rule search_engine_comet:
    input:
        mzml = get_mzml_file,
        fasta = decoy_fasta
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_comet.idXML"
    log:
        logdir + "/{sample_id}/search_engine_comet.log"
    threads: 4
    conda:
        "openms.yaml"
    params:
        comet = comet,
        precursor_mass_tolerance = search_engine_params.get("comet", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = search_engine_params.get("comet", {}).get("fragment_mass_tolerance", 0.02),
        fragment_bin_tolerance = search_engine_params.get("comet", {}).get("fragment_bin_tolerance", 0.02),
        fragment_bin_offset = search_engine_params.get("comet", {}).get("fragment_bin_offset", 0)
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start Comet search for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"{wildcards.sample_id}/comet_{current_time}.sh")
        cmd = [
            params.comet,
            "-in", input.mzml,
            "-out", output.idxml,
            "-database", input.fasta,
            "-threads", str(threads),
            "-precursor_mass_tolerance", str(params.precursor_mass_tolerance),
            "-fragment_mass_tolerance", str(params.fragment_mass_tolerance),
            "-fragment_bin_tolerance", str(params.fragment_bin_tolerance),
            "-fragment_bin_offset", str(params.fragment_bin_offset)
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule search_engine_msgf:
    input:
        mzml = get_mzml_file,
        fasta = decoy_fasta
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_msgf.idXML"
    log:
        logdir + "/{sample_id}/search_engine_msgf.log"
    threads: 4
    conda:
        "openms.yaml"
    params:
        msgf = msgf,
        precursor_mass_tolerance = search_engine_params.get("msgf", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = search_engine_params.get("msgf", {}).get("fragment_mass_tolerance", 0.02),
        isotope_error_range = search_engine_params.get("msgf", {}).get("isotope_error_range", "0,1")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start MSGF+ search for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"{wildcards.sample_id}/msgf_{current_time}.sh")
        cmd = [
            params.msgf,
            "-in", input.mzml,
            "-out", output.idxml,
            "-database", input.fasta,
            "-threads", str(threads),
            "-precursor_mass_tolerance", str(params.precursor_mass_tolerance),
            "-fragment_mass_tolerance", str(params.fragment_mass_tolerance),
            "-isotope_error_range", params.isotope_error_range
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule search_engine_sage:
    input:
        mzml = get_mzml_file,
        fasta = decoy_fasta
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_sage.idXML"
    log:
        logdir + "/{sample_id}/search_engine_sage.log"
    threads: 4
    conda:
        "openms.yaml"
    params:
        sage = sage,
        precursor_mass_tolerance = search_engine_params.get("sage", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = search_engine_params.get("sage", {}).get("fragment_mass_tolerance", 0.02)
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start Sage search for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"{wildcards.sample_id}/sage_{current_time}.sh")
        cmd = [
            params.sage,
            "-in", input.mzml,
            "-out", output.idxml,
            "-database", input.fasta,
            "-threads", str(threads),
            "-precursor_mass_tolerance", str(params.precursor_mass_tolerance),
            "-fragment_mass_tolerance", str(params.fragment_mass_tolerance)
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

# Select search engine based on configuration
if "comet" in search_engines:
    rule search_engine_result:
        input:
            expand(outdir + "/{sid}/{sid}_comet.idXML", sid=samples)
elif "msgf" in search_engines:
    rule search_engine_result:
        input:
            expand(outdir + "/{sid}/{sid}_msgf.idXML", sid=samples)
elif "sage" in search_engines:
    rule search_engine_result:
        input:
            expand(outdir + "/{sid}/{sid}_sage.idXML", sid=samples)
else:
    rule search_engine_result:
        input:
            expand(outdir + "/{sid}/{sid}_comet.idXML", sid=samples)
