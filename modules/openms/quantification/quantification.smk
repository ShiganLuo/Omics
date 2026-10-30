include: "../../common/common.smk"

indir = config.get("indir", "data/protein_inference")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
quantification_method = (config.get("quantification_method") or "lfq").lower()
mzml_dir = config.get("mzml_dir", "")

# PSM-level filtered idXML (5_psm_fdr). ProteomicsLFQ documents its -ids as
# PSM-level identifications and performs protein inference itself; feeding it
# pre-inferred proteins (indir, 6_protein_inference) makes it warn that the
# inference will be invalidated.
psm_fdr_indir = config.get("psm_fdr_indir", "") or indir

# Search database (target+decoy fasta): ProteomicsLFQ needs protein sequences
# for coverage. Wired from the decoy_database output / user decoy_fasta.
search_fasta = config.get("genome", {}).get("search_fasta", "")

# Experimental design TSV (OpenMS design format) shared by ProteomicsLFQ and
# ProteinQuantifier; wired from Params.workflow.expdesign at subworkflow level.
expdesign = config.get("expdesign", "")

# Get parameters
tmt_params = config.get("Params", {}).get("tmt", {})
protein_quant_params = config.get("Params", {}).get("protein_quant", {})
add_triqler_output = config.get("Params", {}).get("add_triqler_output", False)
quantify_decoys = config.get("Params", {}).get("quantify_decoys", False)
protein_fdr_params = config.get("Params", {}).get("protein_fdr", {})

# Get executables
proteomicslfq = config.get("Procedure", {}).get("proteomicslfq") or "ProteomicsLFQ"
proteinquantifier = config.get("Procedure", {}).get("proteinquantifier") or "ProteinQuantifier"
isobaric_analyzer = config.get("Procedure", {}).get("isobaric_analyzer") or "IsobaricAnalyzer"
id_mapper = config.get("Procedure", {}).get("id_mapper") or "IDMapper"
file_merger = config.get("Procedure", {}).get("file_merger") or "FileMerger"

# ProteomicsLFQ protein summarization options. The legacy config keys
# top/average/best_charge_and_score map onto the ProteinQuantification
# subsection of ProteomicsLFQ / the top-N options of ProteinQuantifier.
_pq = protein_quant_params
quant_method_opt = _pq.get("quantification_method") or _pq.get("method") or "feature_intensity"
protein_inference_opt = _pq.get("protein_inference") or "aggregation"
top_n = _pq.get("top", 3)
top_aggregate = _pq.get("average") or _pq.get("top_aggregate") or "median"
best_charge_and_fraction = _pq.get("best_charge_and_fraction",
                                   _pq.get("best_charge_and_score", True))
if protein_inference_opt not in ("aggregation", "bayesian"):
    raise ValueError(
        f"Params.protein_quant.protein_inference={protein_inference_opt!r} is not "
        "supported: use 'aggregation' or 'bayesian'"
    )

# TMT / isobaric extraction options (IsobaricAnalyzer). The legacy
# channel_mass_tolerance/channel_annotation keys map onto
# extraction:reporter_mass_shift and the per-channel description options.
tmt_type = tmt_params.get("type") or "tmt16plex"
reporter_mass_shift = tmt_params.get("reporter_mass_shift",
                                     tmt_params.get("channel_mass_tolerance", 0.003))
isotope_correction = tmt_params.get("isotope_correction", True)
tmt_normalization = tmt_params.get("normalization", False)
reference_channel = tmt_params.get("reference_channel", "")
correction_matrix = tmt_params.get("correction_matrix", []) or []
channel_descriptions = tmt_params.get("channel_descriptions", {}) or {}
if not channel_descriptions and tmt_params.get("channel_annotation"):
    # legacy flat string form: "126:S1,127N:S2"
    for pair in str(tmt_params["channel_annotation"]).split(","):
        if ":" in pair:
            ch, desc = pair.split(":", 1)
            channel_descriptions[ch.strip()] = desc.strip()

def _plfq_summarization_cmd():
    """ProteomicsLFQ arguments shared by the LFQ and DIA group rules."""
    cmd = [
        "-quantification_method", str(quant_method_opt),
        "-protein_inference", str(protein_inference_opt),
        "-ProteinQuantification:top:N", str(top_n),
        "-ProteinQuantification:top:aggregate", str(top_aggregate),
    ]
    if best_charge_and_fraction:
        cmd.append("-ProteinQuantification:best_charge_and_fraction")
    if protein_fdr_params.get("fdr") is not None:
        cmd += ["-proteinFDR", str(protein_fdr_params["fdr"])]
    if protein_fdr_params.get("picked"):
        cmd += ["-picked_proteinFDR", "true"]
    if quantify_decoys:
        cmd.append("-PeptideQuantification:quantify_decoys")
    return cmd

rule quantification_lfq:
    input:
        mzml = expand(mzml_dir + "/{sample_id}/{sample_id}.mzML", sample_id=samples),
        idxml = expand(psm_fdr_indir + "/{sample_id}/{sample_id}_filtered.idXML", sample_id=samples),
        fasta = search_fasta
    output:
        mztab = outdir + "/lfq_quantification.mzTab",
        msstats_csv = outdir + "/lfq_msstats.csv",
        **({"triqler": outdir + "/lfq_quantification.triqler.tsv"} if add_triqler_output else {})
    log:
        logdir + "/quantification/quantification_lfq.log"
    threads: 8
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        proteomicslfq = proteomicslfq,
        expdesign = expdesign,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("quantification_lfq", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start LFQ quantification at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"quantification_lfq_{current_time}.sh")
            plain_script = os.path.join(outdir, f"quantification_lfq.sh")
            cmd = [
                params.proteomicslfq,
                "-in", *[str(x) for x in input.mzml],
                "-ids", *[str(x) for x in input.idxml],
                "-out", output.mztab,
                "-out_msstats", output.msstats_csv,
                "-fasta", str(input.fasta),
            ]
            if add_triqler_output:
                cmd += ["-out_triqler", output.triqler]
            if params.expdesign:
                cmd += ["-design", str(params.expdesign)]
            cmd += ["-threads", str(threads)]
            cmd += _plfq_summarization_cmd()
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'LFQ quantification completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during LFQ quantification: {str(e)}\n")
            raise RuntimeError(f"Error during LFQ quantification: {str(e)}\n")

rule quantification_tmt_reporter:
    input:
        mzml = mzml_dir + "/{sample_id}/{sample_id}.mzML",
        idxml = indir + "/{sample_id}/{sample_id}_protein.idXML",
    output:
        reporter = outdir + "/{sample_id}/{sample_id}.reporter.consensusXML",
        mapped = outdir + "/{sample_id}/{sample_id}.mapped.consensusXML",
    log:
        logdir + "/{sample_id}/quantification_tmt_reporter.log"
    threads: 2
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        isobaric_analyzer = isobaric_analyzer,
        id_mapper = id_mapper,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("quantification_tmt_reporter", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start isobaric reporter extraction for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"quantification_tmt_reporter_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"quantification_tmt_reporter.sh")
            extract_cmd = [
                params.isobaric_analyzer,
                "-in", str(input.mzml),
                "-out", output.reporter,
                "-type", str(tmt_type),
                "-extraction:reporter_mass_shift", str(reporter_mass_shift),
                "-quantification:isotope_correction",
                "true" if isotope_correction else "false",
            ]
            if tmt_normalization:
                extract_cmd.append("-quantification:normalization")
            if reference_channel:
                extract_cmd += [f"-{tmt_type}:reference_channel", str(reference_channel)]
            if correction_matrix:
                extract_cmd += [f"-{tmt_type}:correction_matrix",
                                *[str(x) for x in correction_matrix]]
            for channel, description in channel_descriptions.items():
                extract_cmd += [f"-{tmt_type}:channel_{channel}_description", str(description)]
            map_cmd = [
                params.id_mapper,
                "-in", output.reporter,
                "-id", str(input.idxml),
                "-out", output.mapped,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in extract_cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in map_cmd) + "\n")
                f.write(f"echo 'Isobaric reporter extraction for {wildcards.sample_id} completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during isobaric reporter extraction: {str(e)}\n")
            raise RuntimeError(f"Error during isobaric reporter extraction: {str(e)}\n")

rule quantification_tmt:
    input:
        mapped = expand(outdir + "/{sample_id}/{sample_id}.mapped.consensusXML",
                        sample_id=samples),
        protein_groups = expand(indir + "/{sample_id}/{sample_id}_protein.idXML",
                                sample_id=samples)
    output:
        # Per-plex protein tables stay SEPARATE: FileMerger folds the 16
        # channel maps of every plex into 16 shared maps named
        # "mergedConsensusXMLFile", which destroys sample identity (pooled
        # abundances, mzTab PSM->file mapping and MSstatsConverter design
        # matching all break on it). Each plex is quantified on its own
        # mapped consensusXML + its own protein inference; the per-plex CSVs
        # are then stacked into one long table (sample x channel x protein).
        proteins_long = outdir + "/tmt_protein_abundances.csv",
        consensus_tables = expand(outdir + "/{sample_id}/{sample_id}_proteins.csv",
                                  sample_id=samples),
        **({"mztab": outdir + "/tmt_quantification.mzTab"}
           if config.get("expdesign", "") else {}),
    log:
        logdir + "/quantification/quantification_tmt.log"
    threads: 8
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        proteinquantifier = proteinquantifier,
        expdesign = config.get("expdesign", ""),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("quantification_tmt", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start TMT quantification at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"quantification_tmt_{current_time}.sh")
            plain_script = os.path.join(outdir, f"quantification_tmt.sh")
            quant_cmds = []
            csvs = []
            for sample_id, mapped_cxml, protein_idxml in zip(
                    samples, input.mapped, input.protein_groups):
                per_csv = str(outdir) + f"/{sample_id}/{sample_id}_proteins.csv"
                cmd = [
                    params.proteinquantifier,
                    "-in", str(mapped_cxml),
                    "-protein_groups", str(protein_idxml),
                    "-out", per_csv,
                    "-top:N", str(top_n),
                    "-top:aggregate", str(top_aggregate),
                ]
                if best_charge_and_fraction:
                    cmd.append("-best_charge_and_fraction")
                if params.expdesign:
                    per_mztab = str(outdir) + f"/{sample_id}/{sample_id}_quant.mzTab"
                    cmd += ["-mztab", per_mztab]
                quant_cmds.append(" ".join(shlex.quote(str(x)) for x in cmd))
                csvs.append(per_csv)
            # Stack the per-plex tables into one long table:
            # plex,sample,channel,protein,n_proteins,protein_score,n_peptides,abundance
            # abundance_<k> maps to TMT channel k in IsobaricAnalyzer's fixed
            # order (126, 127N, 127C, ..., 134N for tmt16plex).
            channels = {
                "tmt16plex": ["126", "127N", "127C", "128N", "128C", "129N",
                              "129C", "130N", "130C", "131N", "131C", "132N",
                              "132C", "133N", "133C", "134N"],
                "tmt10plex": ["126", "127N", "127C", "128N", "128C", "129N",
                              "129C", "130N", "131N"],
                "tmt6plex": ["126", "127N", "127C", "128N", "128C", "131N"],
            }.get(str(tmt_type), [str(i + 1) for i in range(16)])
            stack_py = (
                "import csv, sys, glob, os\n"
                "out_path, channels_arg = sys.argv[1], sys.argv[2]\n"
                "plex_dirs = sys.argv[3:]\n"
                "channels = channels_arg.split(',')\n"
                "out = open(out_path, 'w', newline='')\n"
                "w = csv.writer(out, delimiter='\\t')\n"
                "w.writerow(['plex','sample','channel','protein','n_proteins',"
                "'protein_score','n_peptides','abundance'])\n"
                "for plex_dir in sorted(plex_dirs):\n"
                "    plex = os.path.basename(plex_dir)\n"
                "    sample = plex  # one plex == one multiplexed run == one sample file\n"
                "    f = glob.glob(os.path.join(plex_dir, '*_proteins.csv'))\n"
                "    if not f:\n"
                "        sys.exit(f'no quant table in {plex_dir}')\n"
                "    rows = list(csv.reader(open(f[0]), delimiter='\\t'))\n"
                "    hdr = None\n"
                "    for i, r in enumerate(rows):\n"
                "        if r and r[0] == 'protein':\n"
                "            hdr = (i, r)\n"
                "    if hdr is None:\n"
                "        sys.exit(f'no header in {f[0]}')\n"
                "    cols = {n: k for k, n in enumerate(hdr[1])}\n"
                "    for r in rows[hdr[0] + 1:]:\n"
                "        if not r or not r[cols['protein']]:\n"
                "            continue\n"
                "        for ch_i, ch in enumerate(channels):\n"
                "            key = f'abundance_{ch_i + 1}'\n"
                "            if key not in cols or cols[key] >= len(r):\n"
                "                continue\n"
                "            val = r[cols[key]]\n"
                "            if val == '' or float(val) == 0.0:\n"
                "                continue\n"
                "            w.writerow([plex, sample, ch, r[cols['protein']],\n"
                "                        r[cols['n_proteins']], r[cols['protein_score']],\n"
                "                        r[cols['n_peptides']], val])\n"
                "out.close()\n"
            )
            stack_cmd = (
                # arg order MUST match stack_py's parsing:
                #   argv[1]=out, argv[2]=channels(csv), argv[3:]=plex dirs
                "python3 -c {} {} {} {}".format(
                    shlex.quote(stack_py),
                    shlex.quote(str(output.proteins_long)),
                    shlex.quote(",".join(channels)),
                    " ".join(shlex.quote(str(outdir) + "/" + s) for s in samples),
                )
            )
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                for cmd_line in quant_cmds:
                    f.write(cmd_line + "\n")
                f.write(stack_cmd + "\n")
                f.write(f"echo 'TMT quantification completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during TMT quantification: {str(e)}\n")
            raise RuntimeError(f"Error during TMT quantification: {str(e)}\n")

rule quantification_dia:
    input:
        mzml = expand(mzml_dir + "/{sample_id}/{sample_id}.mzML", sample_id=samples),
        idxml = expand(psm_fdr_indir + "/{sample_id}/{sample_id}_filtered.idXML", sample_id=samples),
        fasta = search_fasta
    output:
        mztab = outdir + "/dia_quantification.mzTab",
        msstats_csv = outdir + "/dia_msstats.csv",
        **({"triqler": outdir + "/dia_quantification.triqler.tsv"} if add_triqler_output else {})
    log:
        logdir + "/quantification/quantification_dia.log"
    threads: 8
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        proteomicslfq = proteomicslfq,
        expdesign = expdesign,
    run:
        log_path = str(log)
        rule_logger = setup_logger("quantification_dia", log_file=log_path)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DIA quantification at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"quantification_dia_{current_time}.sh")
            plain_script = os.path.join(outdir, f"quantification_dia.sh")
            cmd = [
                params.proteomicslfq,
                "-in", *[str(x) for x in input.mzml],
                "-ids", *[str(x) for x in input.idxml],
                "-out", output.mztab,
                "-out_msstats", output.msstats_csv,
                "-fasta", str(input.fasta),
            ]
            if add_triqler_output:
                cmd += ["-out_triqler", output.triqler]
            if params.expdesign:
                cmd += ["-design", str(params.expdesign)]
            cmd += ["-threads", str(threads)]
            cmd += _plfq_summarization_cmd()
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'DIA quantification completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during DIA quantification: {str(e)}\n")
            raise RuntimeError(f"Error during DIA quantification: {str(e)}\n")

# Panel rule: which quantification products the workflow registers is decided
# by node.py through the output paths it writes into outfiles; this rule only
# aggregates the products of the configured quantification method.
# tmt: the merged-design mzTab needs an experimental design (plex -> channel
# -> condition mapping); without one the per-plex protein tables and the long
# abundance table are the deliverables.
_RESULT_MZTAB = {
    "tmt": outdir + "/tmt_quantification.mzTab" if config.get("expdesign", "") else outdir + "/tmt_protein_abundances.csv",
    "lfq": outdir + "/lfq_quantification.mzTab",
    "dia": outdir + "/dia_quantification.mzTab",
}

rule quantification_result:
    input:
        mztab = _RESULT_MZTAB.get(quantification_method, outdir + "/lfq_quantification.mzTab")
