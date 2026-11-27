include: "../../common/common.smk"

indir = config.get("indir", "data/raw2mzml")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
search_fasta = config.get("genome", {}).get("search_fasta", "")

ms_params = config.get("Params", {}).get("msfragger", {})
threads = ms_params.get("threads", 16)
label = ms_params.get("label", "tmt16")
license_key = ms_params.get("license_key") or ""
msfragger_jar = ms_params.get("jar") or ""
java = ms_params.get("java") or "java"
java_xmx = ms_params.get("java_xmx", "32G")
if label not in ("tmt16", "label-free"):
    raise ValueError(
        f"Params.msfragger.label={label!r} is not supported: use 'tmt16' or 'label-free'"
    )

msfragger = config.get("Procedure", {}).get("msfragger") or "msfragger"
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
    threads: rule_res("msfragger_run").get("threads", 16)
    benchmark: logdir + "/{sample_id}/msfragger.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "msfragger_run", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../fragpipe.yaml"
    container:
        sif("../fragpipe.yaml", "msfragger_run")
    params:
        msfragger = msfragger,
        make_params = make_params,
        label = label,
        license_key = license_key,
        jar = msfragger_jar,
        java = java,
        java_xmx = java_xmx,
        workdir = outdir + "/{sample_id}",
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("msfragger", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start MSFragger search for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)

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
            script = receipt_script(f"msfragger_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"msfragger.sh")
            if params.jar:
                cmd = [params.java, f"-Xmx{params.java_xmx}", "-jar", str(params.jar)]
            else:
                cmd = [params.msfragger, f"-Xmx{params.java_xmx}"]
                if params.license_key:
                    cmd += ["--key", str(params.license_key)]
                else:
                    raise ValueError(
                        "MSFragger >= 4.x requires a free academic license key: "
                        "set 'Params.msfragger.license_key' in the workflow config "
                        "JSON (obtain at msfragger-upgrader.nesvilab.org/upgrader/), "
                        "or point 'Params.msfragger.jar' at an MSFragger-3.x jar "
                        "(no key required), then rerun."
                    )
            cmd += [params_file, local_mzml]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'MSFragger search for {wildcards.sample_id} completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during MSFragger search: {str(e)}\n")
            raise RuntimeError(f"Error during MSFragger search: {str(e)}\n")

rule msfragger_idxml:
    input:
        pepxml = outdir + "/{sample_id}/{sample_id}.pepXML",
        mzml = indir + "/{sample_id}/{sample_id}.mzML",
        fasta = lambda wildcards: _require(search_fasta, "genome.search_fasta", "msfragger"),
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_msfragger.idXML",
    log:
        logdir + "/{sample_id}/msfragger_idxml.log"
    threads: rule_res("msfragger_idxml").get("threads", 1)
    benchmark: logdir + "/{sample_id}/msfragger_idxml.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "msfragger_idxml", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../../openms/openms.yaml"
    container:
        sif("../../openms/openms.yaml", "msfragger_idxml")
    params:
        converter = config.get("Procedure", {}).get("id_converter") or "IDFileConverter",
        indexer = config.get("Procedure", {}).get("peptide_indexer") or "PeptideIndexer",
        unmatched_action = config.get("Params", {}).get("msfragger", {}).get("unmatched_action", "remove"),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("msfragger_idxml", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start pepXML -> idXML conversion for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"msfragger_idxml_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"msfragger_idxml.sh")
            tmp_pepxml = os.path.join(sample_outdir, f"{wildcards.sample_id}.sanitized.pepXML")
            tmp_idxml = str(output.idxml) + ".tmp"
            sanitize_cmd = (
                "sed 's/name=\"stricttrypsin\"/name=\"Trypsin\"/' {} > {}".format(
                    shlex.quote(str(input.pepxml)), shlex.quote(tmp_pepxml)
                )
            )
            convert_cmd = [
                params.converter,
                "-in", tmp_pepxml,
                "-out", tmp_idxml,
                "-out_type", "idXML",
            ]
            annotate_expr = (
                "s|^\\([[:space:]]*<ProteinIdentification [^>]*>\\)$|\\1\\n"
                "\\t\\t\\t<UserParam type=\"stringList\" name=\"spectra_data\" "
                "value=\"[" + str(input.mzml) + "]\"/>|"
            )
            annotate_cmd = (
                "sed -i '/name=\"spectra_data\"/d' {f}\nsed -i {expr} {f}".format(
                    f=shlex.quote(str(output.idxml)),
                    expr=shlex.quote(annotate_expr),
                )
            )
            index_cmd = [
                params.indexer,
                "-in", tmp_idxml,
                "-fasta", str(input.fasta),
                "-out", str(output.idxml),
            ]
            if params.unmatched_action:
                index_cmd += ["-unmatched_action", str(params.unmatched_action)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(sanitize_cmd + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in convert_cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in index_cmd) + "\n")
                f.write(annotate_cmd + "\n")
                f.write(f"rm -f {shlex.quote(tmp_idxml)} {shlex.quote(tmp_pepxml)}\n")
                f.write(f"echo 'pepXML -> idXML conversion for {wildcards.sample_id} completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during pepXML -> idXML conversion: {str(e)}\n")
            raise RuntimeError(f"Error during pepXML -> idXML conversion: {str(e)}\n")

rule msfragger_result:
    input:
        expand(outdir + "/{sample_id}/{sample_id}.pepXML", sample_id=samples),
        expand(outdir + "/{sample_id}/{sample_id}_msfragger.idXML", sample_id=samples),
