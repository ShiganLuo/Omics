include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
mzML_dir = config.get("mzML_dir", "data/mzml")
decoy_dir = config.get("decoy_dir", "data/decoy_database")

search_engine = config.get("Params", {}).get("search_engine", {}).get("engine") or "comet"

# Get executables
comet = config.get("Procedure", {}).get("comet") or "CometAdapter"
msgf = config.get("Procedure", {}).get("msgf") or "MSGFPlusAdapter"
sage = config.get("Procedure", {}).get("sage") or "SageAdapter"

# Get decoy database
decoy_fasta = config.get("genome", {}).get("decoy_fasta")

def get_input_for_search_engine(wildcards):
    in_dict = {}
    in_dict["mzml"] = mzML_dir + f"/{wildcards.sample_id}/{wildcards.sample_id}.mzML"
    if decoy_fasta and os.path.exists(decoy_fasta):
        in_dict["fasta"] = decoy_fasta
    else:
        in_dict["fasta"] = decoy_dir + "/genome_decoy.fasta"
    return in_dict
        

# Modification defaults follow the sample chemistry implied by the
# quantification method: labelled samples need the label as modification,
# label-free (lfq/dia) samples only need the standard tryptic set. Explicit
# config values always win.
#
# CAUTION (verified against real TMT16 data + pepXML ground truth): CometAdapter
# SILENTLY IGNORES isobaric label mods configured as FIXED modifications — with
# TMTpro fixed the search output was identical to a no-TMT control (0/40 known
# PSMs recovered), while TMTpro as a VARIABLE modification recovered 13/40.
# N-terminal labelling is also frequently incomplete, which a variable mod
# handles naturally. Label mods therefore always end up in variable_modifications.
_MOD_DEFAULTS = {
    "tmt": ["Carbamidomethyl (C)"],
    "lfq": ["Carbamidomethyl (C)"],
    "dia": ["Carbamidomethyl (C)"],
}
_VAR_MOD_DEFAULTS = {
    "tmt": ["TMTpro (N-term)", "TMTpro (K)", "Oxidation (M)"],
    "lfq": ["Oxidation (M)"],
    "dia": ["Oxidation (M)"],
}
_LABEL_TOKENS = ("tmt", "itraq")

def _normalize_label_mods(fixed, variable, method):
    """Force isobaric label mods into the variable-modification list, first.

    Ordering matters: when another variable mod (e.g. Oxidation (M)) is
    listed BEFORE the N-terminal label mod, CometAdapter silently drops the
    label mod and the search degrades to a no-label control (verified:
    13/40 known PSMs recovered with the label mod first vs 0/40 otherwise).
    """
    fixed = [str(m) for m in (fixed or [])]
    variable = [str(m) for m in (variable or [])]
    moved = [m for m in fixed if any(t in m.lower() for t in _LABEL_TOKENS)]
    if moved:
        fixed = [m for m in fixed if m not in moved]
        variable += [m for m in moved if m not in variable]
        logger.warning(
            "[search_engine] isobaric label mods must be VARIABLE for Comet "
            f"(fixed ones are silently ignored): moved {moved} to "
            "variable_modifications"
        )
    if method == "tmt" and not any(
        t in m.lower() for t in _LABEL_TOKENS for m in fixed + variable
    ):
        variable += ["TMTpro (N-term)", "TMTpro (K)"]
        logger.warning(
            "[search_engine] quantification_method='tmt' but no label mods "
            "configured: adding TMTpro (N-term)/(K) as variable modifications"
        )
    label = [m for m in variable if any(t in m.lower() for t in _LABEL_TOKENS)]
    other = [m for m in variable if m not in label]
    return fixed, label + other

_quant_method = (config.get("quantification_method") or "").lower()
_comet_params = config.get("Params", {}).get("search_engine", {}).get("comet", {})
_fixed_mods, _variable_mods = _normalize_label_mods(
    _comet_params.get("fixed_modifications") or _MOD_DEFAULTS.get(
        _quant_method, ["Carbamidomethyl (C)"]),
    _comet_params.get("variable_modifications") or _VAR_MOD_DEFAULTS.get(
        _quant_method, ["Oxidation (M)"]),
    _quant_method,
)
# Label mods imply isobaric reporter ions must be present in the spectra;
# verify at runtime so a label/data mismatch fails loudly at the search step
# instead of silently producing noise PSMs (see check_reporter_ions).
_expect_reporters = any(
    "tmt" in str(m).lower() or "itraq" in str(m).lower()
    for m in list(_fixed_mods) + list(_variable_mods)
)

def _assert_label_matches_data(mzml_path):
    """Fail when isobaric label mods are configured but spectra carry no reporters."""
    if not _expect_reporters:
        return
    has_reporters, detail = check_reporter_ions(str(mzml_path))
    if not has_reporters:
        raise ValueError(
            f"isobaric label modifications ({_fixed_mods + _variable_mods}) are "
            f"configured but {mzml_path} shows no TMT/iTRAQ reporter ions "
            f"({detail}). The spectra are NOT isobarically labelled: set "
            f"Params.workflow.quantification_method to 'lfq' (or provide "
            f"fixed_modifications without label mods) — searching labelled "
            f"peptide masses against unlabelled spectra yields only noise."
        )

rule search_engine_comet:
    input:
        unpack(get_input_for_search_engine)
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_comet.idXML"
    log:
        logdir + "/{sample_id}/search_engine_comet.log"
    threads: 4
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        comet = comet,
        precursor_mass_tolerance = config.get("Params").get("search_engine",{}).get("comet", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = config.get("Params").get("search_engine",{}).get("comet", {}).get("fragment_mass_tolerance", 0.02),
        fragment_bin_offset = config.get("Params").get("search_engine",{}).get("comet", {}).get("fragment_bin_offset", 0),
        unmatched_action = config.get("Params").get("search_engine",{}).get("comet", {}).get("unmatched_action", "remove"),
        fixed_modifications = _fixed_mods,
        variable_modifications = _variable_mods
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("search_engine_comet",log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start Comet search for sample {wildcards.sample_id} at {current_time}")
            _assert_label_matches_data(input.mzml)
            script = os.path.join(RECEIPT_DIR, f"comet_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/comet.sh")
            cmd = [
                params.comet,
                "-in", input.mzml,
                "-out", output.idxml,
                "-database", input.fasta,
                "-threads", str(threads),
                "-precursor_mass_tolerance", str(params.precursor_mass_tolerance),
                "-fragment_mass_tolerance", str(params.fragment_mass_tolerance),
                "-fragment_bin_offset", str(params.fragment_bin_offset)
            ]
            # PeptideIndexer aborts (exit 13) on unmatched peptides by default;
            # they cannot be used for FDR/quant anyway, so drop them.
            if params.unmatched_action:
                cmd += ["-PeptideIndexing:unmatched_action", str(params.unmatched_action)]
            # Labelled samples (e.g. TMT) need their label modification set as
            # fixed, otherwise real peptides never match (mass shift lost).
            for mod in params.fixed_modifications:
                cmd += ["-fixed_modifications", str(mod)]
            for mod in params.variable_modifications:
                cmd += ["-variable_modifications", str(mod)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "rule search_engine_comet called for {wildcards.sample_id} was successfully completed"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule search_engine_comet was call failed,error: {e}")
            raise RuntimeError(f"rule search_engine_comet was call failed,error: {e}")
        

rule search_engine_msgf:
    input:
        unpack(get_input_for_search_engine)
    output:
        idxml = outdir + "/{sample_id}/{sample_id}_msgf.idXML"
    log:
        logdir + "/{sample_id}/search_engine_msgf.log"
    threads: 4
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        msgf = msgf,
        precursor_mass_tolerance = config.get("Params").get("search_engine",{}).get("msgf", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = config.get("Params").get("search_engine",{}).get("msgf", {}).get("fragment_mass_tolerance", 0.02),
        isotope_error_range = config.get("Params").get("search_engine",{}).get("msgf", {}).get("isotope_error_range", "0,1")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("search_engine_msgf",log_file = log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start MSGF+ search for sample {wildcards.sample_id} at {current_time}")
            _assert_label_matches_data(input.mzml)
            script = os.path.join(RECEIPT_DIR, f"msgf_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/msgf.sh")
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
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "rule search_engine_msgf called for {wildcards.sample_id} was successfully completed"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule search_engine_msgf was call failed,error: {e}")
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
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        sage = sage,
        precursor_mass_tolerance = config.get("Params").get("search_engine",{}).get("sage", {}).get("precursor_mass_tolerance", 20),
        fragment_mass_tolerance = config.get("Params").get("search_engine",{}).get("sage", {}).get("fragment_mass_tolerance", 0.02)
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("search_engine_sage",log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start Sage search for sample {wildcards.sample_id} at {current_time}")
            _assert_label_matches_data(input.mzml)
            script = os.path.join(RECEIPT_DIR, f"sage_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/sage.sh")
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
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "rule search_engine_sage called for {wildcards.sample_id} was successfully completed"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule search_engine_sage was call failed,error: {e}")
            raise RuntimeError(f"rule search_engine_sage was call failed,error: {e}")
