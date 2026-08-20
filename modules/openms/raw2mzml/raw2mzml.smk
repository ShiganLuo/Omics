include: "../../common/common.smk"

indir = config.get("indir", "data/raw")
outdir = config.get("outdir", "output/raw2mzml")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
raw_files = config.get("raw_files", [])

converter = config.get("Params", {}).get("raw_to_mzml", {}).get("converter") or "thermorawfileparser"
converter_mode = config.get("Params", {}).get("raw_to_mzml", {}).get("mode") or "ThermoRawFileParser"
converter_args = config.get("Params", {}).get("raw_to_mzml", {}).get("args") or ""
extra_filter = config.get("Params", {}).get("raw_to_mzml", {}).get("peak_picking", True)


def get_raw_input(wildcards):
    if raw_files and wildcards.sample_id in samples:
        idx = samples.index(wildcards.sample_id)
        if idx < len(raw_files):
            return os.path.realpath(raw_files[idx])
    candidates = [
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.raw",
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.RAW",
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.raw.gz",
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.RAW.gz",
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.mzML",
    ]
    for path in candidates:
        if os.path.exists(path):
            return os.path.realpath(path)
    return candidates[0]


rule raw2mzml:
    input:
        infile = get_raw_input
    output:
        mzml = outdir + "/{sample_id}/{sample_id}.mzML"
    log:
        logdir + "/{sample_id}/raw2mzml.log"
    threads: 2
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        converter = converter,
        converter_mode = converter_mode,
        converter_args = converter_args,
        peak_picking = extra_filter
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger = setup_logger("raw2mzml",log_file=log_path)
            rule_logger.info(f"Start raw2mzml for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"raw2mzml_{current_time}.sh")
            output_dir = os.path.dirname(output.mzml)
            output_name = os.path.basename(output.mzml)
            if input.infile.endswith(".raw") or input.infile.endswith(".RAW"):
                cmd = [params.converter]
                if params.converter_args:
                    cmd.extend(shlex.split(str(params.converter_args)))
                if params.converter_mode == "msconvert":
                    cmd.extend([str(input.infile), "--mzML"])
                    if params.peak_picking:
                        cmd.extend(["--filter", "peakPicking true 1-"])
                    cmd.extend(["-o", output_dir, "--outfile", output_name])
                elif params.converter_mode == "thermorawfileparser":
                    cmd.extend(["-i", str(input.infile), "-b", str(output.mzml), "-f", "1"])
                    if not params.peak_picking:
                        cmd.append("-p")
                else:
                    cmd.extend([str(input.infile), str(output.mzml)])
            elif input.infile.endswith(".mzML"):
                rule_logger.info(f"Input file {input.infile} is already in mzML format. Copying to output directory.")
                cmd = ["ln", "-s", str(input.infile), str(output.mzml)]
            else:
                raise ValueError(f"Unsupported input file format: {input.infile}")
            
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(x) for x in cmd) + "\n")
                f.write(f'echo "raw2mzml completed for sample {wildcards.sample_id} at $(date)"\n')
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"rule raw2mzml was call failed,error: {e}")
            raise RuntimeError(f"rule raw2mzml was call failed,error: {e}")

rule raw2mzml_result:
    input:
        expand(outdir + "/{sample_id}/{sample_id}.mzML", sample_id=samples)
