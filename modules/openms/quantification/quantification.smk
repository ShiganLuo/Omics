include: "../../common/common.smk"
proteomicslfq = config.get("Procedure", {}).get("proteomicslfq") or "ProteomicsLFQ"
proteinquantifier = config.get("Procedure", {}).get("proteinquantifier") or "ProteinQuantifier"
isobaric_analyzer = config.get("Procedure", {}).get("isobaric_analyzer") or "IsobaricAnalyzer"
id_mapper = config.get("Procedure", {}).get("id_mapper") or "IDMapper"
file_merger = config.get("Procedure", {}).get("file_merger") or "FileMerger"

indir = config.get("indir", "data/protein_inference")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
# 指令派生名(input: 按方法取聚合目标;config 原名 quantification_method),
# 赋值写 config 全量路径 Params.workflow.quantification_method。
quantification_method = (config.get("Params", {}).get("workflow", {}).get("quantification_method") or "lfq").lower()
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
# 参数精确到 rule(2026-11 规范):config 全量路径字面链,共享读点拆成
# 每 rule 一个 parse 期变量,键挂在消费 rule 名下,无 module 级共享回退:
#   Params.tmt.<rule_name>.<key>            # IsobaricAnalyzer 提取/通道映射
#   Params.protein_quant.<rule_name>.<key>  # 汇总参数(lfq/dia/tmt 各自独立)
#   Params.protein_fdr.<rule_name>.<key>    # protein FDR(lfq/dia)
#   Params.<rule_name>.add_triqler_output   # triqler 输出开关(lfq/dia)
#   Params.<rule_name>.quantify_decoys      # decoy 定量(lfq/dia)
# 旧键 method/average/best_charge_and_score/channel_mass_tolerance/
# channel_annotation 为兼容回退键,随 module 级共享读取一并删除。

# add_triqler_output(Params.<rule_name>.add_triqler_output)、quantify_decoys
# (Params.<rule_name>.quantify_decoys):lfq/dia 两 rule 共享 → 每 rule 独立键,
# output 展开与 CLI 转发共用 parse 期变量。
_add_triqler_lfq = config.get("Params", {}).get("quantification_lfq", {}).get("add_triqler_output", False)
_add_triqler_dia = config.get("Params", {}).get("quantification_dia", {}).get("add_triqler_output", False)
_quantify_decoys_lfq = config.get("Params", {}).get("quantification_lfq", {}).get("quantify_decoys", False)
_quantify_decoys_dia = config.get("Params", {}).get("quantification_dia", {}).get("quantify_decoys", False)

# ProteomicsLFQ protein summarization options (ProteinQuantification subsection
# of ProteomicsLFQ / top-N options of ProteinQuantifier):
# Params.protein_quant.<rule_name>.<key>,每 rule 一个 parse 期变量。lfq/dia 经
# _plfq_summarization_cmd 转发;tmt 只消费 top/top_aggregate/best_charge_and_fraction。
_quant_method_lfq = config.get("Params", {}).get("protein_quant", {}).get("quantification_lfq", {}).get("quantification_method", "feature_intensity")
_inference_lfq = config.get("Params", {}).get("protein_quant", {}).get("quantification_lfq", {}).get("protein_inference", "aggregation")
_top_n_lfq = config.get("Params", {}).get("protein_quant", {}).get("quantification_lfq", {}).get("top", 3)
_top_aggregate_lfq = config.get("Params", {}).get("protein_quant", {}).get("quantification_lfq", {}).get("top_aggregate", "median")
_best_caf_lfq = config.get("Params", {}).get("protein_quant", {}).get("quantification_lfq", {}).get("best_charge_and_fraction", True)
_quant_method_dia = config.get("Params", {}).get("protein_quant", {}).get("quantification_dia", {}).get("quantification_method", "feature_intensity")
_inference_dia = config.get("Params", {}).get("protein_quant", {}).get("quantification_dia", {}).get("protein_inference", "aggregation")
_top_n_dia = config.get("Params", {}).get("protein_quant", {}).get("quantification_dia", {}).get("top", 3)
_top_aggregate_dia = config.get("Params", {}).get("protein_quant", {}).get("quantification_dia", {}).get("top_aggregate", "median")
_best_caf_dia = config.get("Params", {}).get("protein_quant", {}).get("quantification_dia", {}).get("best_charge_and_fraction", True)
_top_n_tmt = config.get("Params", {}).get("protein_quant", {}).get("quantification_tmt", {}).get("top", 3)
_top_aggregate_tmt = config.get("Params", {}).get("protein_quant", {}).get("quantification_tmt", {}).get("top_aggregate", "median")
_best_caf_tmt = config.get("Params", {}).get("protein_quant", {}).get("quantification_tmt", {}).get("best_charge_and_fraction", True)
for _inf_rule, _inf_val in (("quantification_lfq", _inference_lfq),
                           ("quantification_dia", _inference_dia)):
    if _inf_val not in ("aggregation", "bayesian"):
        raise ValueError(
            f"Params.protein_quant.{_inf_rule}.protein_inference={_inf_val!r} is not "
            "supported: use 'aggregation' or 'bayesian'"
        )

# Protein FDR options (ProteomicsLFQ -proteinFDR/-picked_proteinFDR):
# Params.protein_fdr.<rule_name>.fdr/picked,lfq/dia 两 rule 共享 → 每 rule 独立。
_protein_fdr_lfq = config.get("Params", {}).get("protein_fdr", {}).get("quantification_lfq", {}).get("fdr", None)
_protein_fdr_picked_lfq = config.get("Params", {}).get("protein_fdr", {}).get("quantification_lfq", {}).get("picked", False)
_protein_fdr_dia = config.get("Params", {}).get("protein_fdr", {}).get("quantification_dia", {}).get("fdr", None)
_protein_fdr_picked_dia = config.get("Params", {}).get("protein_fdr", {}).get("quantification_dia", {}).get("picked", False)

# TMT / isobaric extraction options (IsobaricAnalyzer):
# Params.tmt.<rule_name>.<key>。type 被 quantification_tmt_reporter、
# quantification_tmt 两 rule 共享 → 每 rule 独立键;其余仅
# quantification_tmt_reporter 消费,同样挂 rule 名下。旧键
# channel_mass_tolerance/channel_annotation 为兼容回退键,已删除。
_tmt_type_reporter = config.get("Params", {}).get("tmt", {}).get("quantification_tmt_reporter", {}).get("type", "tmt16plex")
_reporter_mass_shift = config.get("Params", {}).get("tmt", {}).get("quantification_tmt_reporter", {}).get("reporter_mass_shift", 0.003)
_isotope_correction = config.get("Params", {}).get("tmt", {}).get("quantification_tmt_reporter", {}).get("isotope_correction", True)
_tmt_normalization = config.get("Params", {}).get("tmt", {}).get("quantification_tmt_reporter", {}).get("normalization", False)
_reference_channel = config.get("Params", {}).get("tmt", {}).get("quantification_tmt_reporter", {}).get("reference_channel", "")
_correction_matrix = config.get("Params", {}).get("tmt", {}).get("quantification_tmt_reporter", {}).get("correction_matrix", [])
_channel_descriptions = config.get("Params", {}).get("tmt", {}).get("quantification_tmt_reporter", {}).get("channel_descriptions", {})
_tmt_type_tmt = config.get("Params", {}).get("tmt", {}).get("quantification_tmt", {}).get("type", "tmt16plex")

def _plfq_summarization_cmd(quant_method, protein_inference, top_n, top_aggregate,
                            best_charge_and_fraction, fdr, picked, quantify_decoys):
    """ProteomicsLFQ arguments shared by the LFQ and DIA group rules.

    参数已按 rule 独立解析后逐 key 传入(无共享回退)。
    """
    cmd = [
        "-quantification_method", str(quant_method),
        "-protein_inference", str(protein_inference),
        "-ProteinQuantification:top:N", str(top_n),
        "-ProteinQuantification:top:aggregate", str(top_aggregate),
    ]
    if best_charge_and_fraction:
        cmd.append("-ProteinQuantification:best_charge_and_fraction")
    if fdr is not None:
        cmd += ["-proteinFDR", str(fdr)]
    if picked:
        cmd += ["-picked_proteinFDR", "true"]
    if quantify_decoys:
        cmd.append("-PeptideQuantification:quantify_decoys")
    return cmd

# 每 rule 一个 parse 期变量(lfq/dia 共享构造拆分)
_plfq_cmd_lfq = _plfq_summarization_cmd(
    _quant_method_lfq, _inference_lfq, _top_n_lfq, _top_aggregate_lfq,
    _best_caf_lfq, _protein_fdr_lfq, _protein_fdr_picked_lfq, _quantify_decoys_lfq)
_plfq_cmd_dia = _plfq_summarization_cmd(
    _quant_method_dia, _inference_dia, _top_n_dia, _top_aggregate_dia,
    _best_caf_dia, _protein_fdr_dia, _protein_fdr_picked_dia, _quantify_decoys_dia)

rule quantification_lfq:
    input:
        mzml = expand(mzml_dir + "/{sample_id}/{sample_id}.mzML", sample_id=samples),
        idxml = expand(psm_fdr_indir + "/{sample_id}/{sample_id}_filtered.idXML", sample_id=samples),
        fasta = search_fasta
    output:
        mztab = outdir + "/lfq_quantification.mzTab",
        msstats_csv = outdir + "/lfq_msstats.csv",
        **({"triqler": outdir + "/lfq_quantification.triqler.tsv"} if _add_triqler_lfq else {})
    log:
        logdir + "/quantification/quantification_lfq.log"
    threads: rule_res("quantification_lfq").get("threads", 8)
    benchmark: logdir + "/quantification/quantification_lfq.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "quantification_lfq", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml", "quantification_lfq")
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
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"quantification_lfq.sh")
            cmd = [
                params.proteomicslfq,
                "-in", *[str(x) for x in input.mzml],
                "-ids", *[str(x) for x in input.idxml],
                "-out", output.mztab,
                "-out_msstats", output.msstats_csv,
                "-fasta", str(input.fasta),
            ]
            if _add_triqler_lfq:
                cmd += ["-out_triqler", output.triqler]
            if params.expdesign:
                cmd += ["-design", str(params.expdesign)]
            cmd += ["-threads", str(threads)]
            cmd += _plfq_cmd_lfq
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'LFQ quantification completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
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
    threads: rule_res("quantification_tmt_reporter").get("threads", 2)
    benchmark: logdir + "/{sample_id}/quantification_tmt_reporter.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "quantification_tmt_reporter", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml", "quantification_tmt_reporter")
    params:
        proteinquantifier = proteinquantifier,
        isobaric_analyzer = isobaric_analyzer,
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
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"quantification_tmt_reporter.sh")
            extract_cmd = [
                params.isobaric_analyzer,
                "-in", str(input.mzml),
                "-out", output.reporter,
                "-type", str(_tmt_type_reporter),
                "-extraction:reporter_mass_shift", str(_reporter_mass_shift),
                "-quantification:isotope_correction",
                "true" if _isotope_correction else "false",
            ]
            if _tmt_normalization:
                extract_cmd.append("-quantification:normalization")
            if _reference_channel:
                extract_cmd += [f"-{_tmt_type_reporter}:reference_channel", str(_reference_channel)]
            if _correction_matrix:
                extract_cmd += [f"-{_tmt_type_reporter}:correction_matrix",
                                *[str(x) for x in _correction_matrix]]
            for channel, description in _channel_descriptions.items():
                extract_cmd += [f"-{_tmt_type_reporter}:channel_{channel}_description", str(description)]
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
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
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
    threads: rule_res("quantification_tmt").get("threads", 8)
    benchmark: logdir + "/quantification/quantification_tmt.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "quantification_tmt", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml", "quantification_tmt")
    params:
        id_mapper = id_mapper,
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
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
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
                    "-top:N", str(_top_n_tmt),
                    "-top:aggregate", str(_top_aggregate_tmt),
                ]
                if _best_caf_tmt:
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
            }.get(str(_tmt_type_tmt), [str(i + 1) for i in range(16)])
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
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
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
        **({"triqler": outdir + "/dia_quantification.triqler.tsv"} if _add_triqler_dia else {})
    log:
        logdir + "/quantification/quantification_dia.log"
    threads: rule_res("quantification_dia").get("threads", 8)
    benchmark: logdir + "/quantification/quantification_dia.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "quantification_dia", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml", "quantification_dia")
    params:
        file_merger = file_merger,
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
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"quantification_dia.sh")
            cmd = [
                params.proteomicslfq,
                "-in", *[str(x) for x in input.mzml],
                "-ids", *[str(x) for x in input.idxml],
                "-out", output.mztab,
                "-out_msstats", output.msstats_csv,
                "-fasta", str(input.fasta),
            ]
            if _add_triqler_dia:
                cmd += ["-out_triqler", output.triqler]
            if params.expdesign:
                cmd += ["-design", str(params.expdesign)]
            cmd += ["-threads", str(threads)]
            cmd += _plfq_cmd_dia
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'DIA quantification completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
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
