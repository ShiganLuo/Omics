include: "../../common/common.smk"

import shutil

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
# Direct jar mode: run MSFragger-3.x.jar via java (paper version, no license
# key needed). Takes precedence over the Procedure.msfragger wrapper.
msfragger_jar = ms_params.get("jar") or ""
java = ms_params.get("java") or "java"
java_xmx = ms_params.get("java_xmx", "32G")
if label not in ("tmt16", "label-free"):
    raise ValueError(
        f"Params.msfragger.label={label!r} is not supported: use 'tmt16' or 'label-free'"
    )

# Get executables (bioconda installs the launcher in lower case: 'msfragger')
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
        jar = msfragger_jar,
        java = java,
        java_xmx = java_xmx,
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
            os.makedirs(os.path.join(logdir, "cmd", "group"), exist_ok=True)
            script = os.path.join(logdir, "cmd", "group", f"msfragger_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"msfragger.sh")
            if params.jar:
                cmd = [params.java, f"-Xmx{params.java_xmx}", "-jar", str(params.jar)]
            else:
                # The conda launcher defaults to -Xmx1g — far too small for a
                # ~500k-protein search DB; pass -Xmx through (it forwards
                # -Xm* args to the JVM).
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
    threads: 1
    conda:
        "../../openms/openms.yaml"
    container:
        sif("../../openms/openms.yaml")
    params:
        converter = config.get("Procedure", {}).get("id_converter") or "IDFileConverter",
        indexer = config.get("Procedure", {}).get("peptide_indexer") or "PeptideIndexer",
        unmatched_action = config.get("Params", {}).get("msfragger", {}).get("unmatched_action", "remove"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("msfragger_idxml", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start pepXML -> idXML conversion for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(os.path.join(logdir, "cmd", "group"), exist_ok=True)
            script = os.path.join(logdir, "cmd", "group", f"msfragger_idxml_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"msfragger_idxml.sh")
            # Staged copy must keep a .pepXML suffix: IDFileConverter has no
            # -in_type in OpenMS 3.x and infers the INPUT type from the
            # extension (it only offers -out_type).
            tmp_pepxml = os.path.join(sample_outdir, f"{wildcards.sample_id}.sanitized.pepXML")
            tmp_idxml = str(output.idxml) + ".tmp"
            # MSFragger writes <sample_enzyme name="stricttrypsin">, which the
            # OpenMS EnzymeDB cannot resolve ("the element 'stricttrypsin' could
            # not be found"); rewrite it to the equivalent standard enzyme name
            # on a staged copy -- the <specificity> element keeps the actual
            # cleavage semantics. The upstream pepXML stays untouched.
            sanitize_cmd = (
                "sed 's/name=\"stricttrypsin\"/name=\"Trypsin\"/' {} > {}".format(
                    shlex.quote(str(input.pepxml)), shlex.quote(tmp_pepxml)
                )
            )
            # -out_type is required: OpenMS would otherwise infer the output
            # type from the .tmp staging extension and abort with
            # "Could not determine output file type!".
            convert_cmd = [
                params.converter,
                "-in", tmp_pepxml,
                "-out", tmp_idxml,
                "-out_type", "idXML",
            ]
            # IDMapper at quantification requires an annotated primary MS run
            # path (idXML meta 'spectra_data'); pepXML -> idXML conversion
            # leaves it empty (even with -mz_file) and later OpenMS tools
            # serialize it as "[]", which IDMapper rejects with
            # "No MS run path annotated in ProteinIdentification.". Inject the
            # mzML run path once here; downstream tools round-trip the value
            # unchanged (verified PSMFeatureExtractor -> PercolatorAdapter ->
            # FalseDiscoveryRate -> Epifany -> IDMapper). Delete-before-insert
            # keeps it single: duplicate 'spectra_data' keys silently overwrite
            # (last one wins), which re-introduces the empty value.
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
            # PeptideIndexer re-annotates target/decoy + protein mapping so the
            # downstream OpenMS chain (PercolatorAdapter/FalseDiscoveryRate)
            # sees the same annotations the CometAdapter path produces.
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
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during pepXML -> idXML conversion: {str(e)}\n")
            raise RuntimeError(f"Error during pepXML -> idXML conversion: {str(e)}\n")


rule msfragger_result:
    input:
        expand(outdir + "/{sample_id}/{sample_id}.pepXML", sample_id=samples),
        expand(outdir + "/{sample_id}/{sample_id}_msfragger.idXML", sample_id=samples),
