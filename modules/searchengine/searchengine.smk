include: "../common/common.smk"
indir = config.get("indir", "data/mzml")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

search_engine = config.get("Params", {}).get("search_engine") or "comet"

# Get executables
comet = config.get("Procedure", {}).get("comet") or "CometAdapter"
msgf = config.get("Procedure", {}).get("msgf") or "MSGFPlusAdapter"
sage = config.get("Procedure", {}).get("sage") or "SageAdapter"

# Get decoy database
decoy_fasta = config.get("genome", {}).get("decoy_fasta")

def get_input_for_search_engine(wildcards):
    in_dict = {}
    in_dict["mzml"] = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.mzML"
    if os.path.exists(decoy_fasta):
        in_dict["fasta"] = decoy_fasta
    else:
        in_dict["fasta"] = indir + "/genome_decoy.fasta"
    return in_dict
        

rule search_engine_comet:
    input:
        unpack(get_input_for_search_engine)
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_comet.idXML"
    log:
        logdir + "/{sample_id}/search_engine_comet.log"
    threads: 4
    conda:
        "searchengine.yaml"
    container:
        sif("searchengine.yaml")
    params:
        comet = comet,
        precursor_mass_tolerance = config.get("Params").get("search_engine",{}).get("comet", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = config.get("Params").get("search_engine",{}).get("comet", {}).get("fragment_mass_tolerance", 0.02),
        fragment_bin_tolerance = config.get("Params").get("search_engine",{}).get("comet", {}).get("fragment_bin_tolerance", 0.02),
        fragment_bin_offset = config.get("Params").get("search_engine",{}).get("comet", {}).get("fragment_bin_offset", 0)
    run:
        log_path = str(log)
        try:
            
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
    unpack(get_input_for_search_engine)
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_msgf.idXML"
    log:
        logdir + "/{sample_id}/search_engine_msgf.log"
    threads: 4
    conda:
        "searchengine.yaml"
    container:
        sif("searchengine.yaml")
    params:
        msgf = msgf,
        precursor_mass_tolerance = config.get("Params").get("search_engine",{}).get("msgf", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = config.get("Params").get("search_engine",{}).get("msgf", {}).get("fragment_mass_tolerance", 0.02),
        isotope_error_range = config.get("Params").get("search_engine",{}).get("msgf", {}).get("isotope_error_range", "0,1")
    run:
        log_path = str(log)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger = setup_logger("search_engine_msgf",log_file = log_path)
            rule_logger.info(f"Start MSGF+ search for sample {wildcards.sample_id} at {current_time}")
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
                f.write(f'echo "rule search_engine_msgf called for {wildcards.sample_id} was successfully completed"\n')
            shell(f"bash {script} >> {log} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"rule search_engine_msgf was call failed,error: {e}")
            raise RuntimeError(f"rule search_engine_msgf was call failed,error: {e}")

rule search_engine_sage:
    input:
        unpack(get_input_for_search_engine)
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_sage.idXML"
    log:
        logdir + "/{sample_id}/search_engine_sage.log"
    threads: 4
    conda:
        "searchengine.yaml"
    container:
        sif("searchengine.yaml")
    params:
        sage = sage,
        precursor_mass_tolerance = config.get("Params").get("search_engine",{}).get("sage", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = config.get("Params").get("search_engine",{}).get("sage", {}).get("fragment_mass_tolerance", 0.02)
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
if search_engine == "comet":
    rule search_engine_result:
        input:
            expand(outdir + "/{sid}/{sid}_comet.idXML", sid=samples)
elif search_engine == "msgf":
    rule search_engine_result:
        input:
            expand(outdir + "/{sid}/{sid}_msgf.idXML", sid=samples)
elif search_engine == "sage":
    rule search_engine_result:
        input:
            expand(outdir + "/{sid}/{sid}_sage.idXML", sid=samples)
else:
    rule search_engine_result:
        input:
            expand(outdir + "/{sid}/{sid}_comet.idXML", sid=samples)
