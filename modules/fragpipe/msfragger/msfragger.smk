include: "../../common/common.smk"

indir = config.get("indir", "data/raw2mzml")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
search_fasta = config.get("genome", {}).get("search_fasta", "")

# Get parameters
ms_params = config.get("Params", {}).get("msfragger", {})
threads = ms_params.get("threads", 16)
label = ms_params.get("label", "tmt16")
license_key = ms_params.get("license_key") or ""
if label not in ("tmt16", "label-free"):
    raise ValueError(
        f"Params.msfragger.label={label!r} is not supported: use 'tmt16' or 'label-free'"
    )

# Get executables
msfragger = config.get("Procedure", {}).get("msfragger") or "MSFragger"
make_params = ROOT_DIR + "/modules/fragpipe/bin/make_fragger_params.py"


rule msfragger_run:
    input:
        mzml = indir + "/{sample_id}/{sample_id}.mzML",
        fasta = lambda wildcards: _require(search_fasta, "genome.search_fasta", "msfragger"),
    output:
        pepxml = outdir + "/{sample_id}/{sample_id}.pepXML",
        tsv = outdir + "/{sample_id}/{sample_id}.tsv",
    log:
        logdir + "/{sample_id}/msfragger.log"
    threads: 16
    conda:
        "../fragpipe.yaml"
    container:
        sif("../fragpipe.yaml")
    params:
        msfragger = msfragger,
        make_params = make_params,
        label = label,
        license_key = license_key,
        workdir = outdir + "/{sample_id}",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("msfragger", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start MSFragger search for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)

            # MSFragger writes results next to its input file: hardlink the
            # mzML into the work dir so outputs never pollute raw2mzml.
            local_mzml = os.path.join(sample_outdir, f"{wildcards.sample_id}.mzML")
            if not os.path.exists(local_mzml):
                try:
                    os.link(str(input.mzml), local_mzml)
                except OSError:
                    shutil.copy(str(input.mzml), local_mzml)

            params_file = os.path.join(sample_outdir, f"{wildcards.sample_id}.fragger.params")
            subprocess.run(
                [sys.executable, str(params.make_params),
                 "--fasta", str(input.fasta),
                 "-o", params_file,
                 "-t", str(threads),
                 "--label", str(params.label)],
                check=True,
            )
            script = os.path.join(sample_outdir, f"msfragger_{current_time}.sh")
            cmd = [params.msfragger]
            if params.license_key:
                cmd += ["--key", str(params.license_key)]
            else:
                raise ValueError(
                    "MSFragger >= 4.x requires a free academic license key: "
                    "set 'Params.msfragger.license_key' in the workflow config "
                    "JSON (obtain at msfragger-upgrader.nesvilab.org/upgrader/), "
                    "then rerun."
                )
            cmd += [params_file, local_mzml]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'MSFragger search for {wildcards.sample_id} completed at {current_time}'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during MSFragger search: {str(e)}\n")
            raise RuntimeError(f"Error during MSFragger search: {str(e)}\n")


rule msfragger_result:
    input:
        expand(outdir + "/{sample_id}/{sample_id}.pepXML", sample_id=samples)
