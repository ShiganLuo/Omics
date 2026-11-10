from snakemake.logging import logger
include: "../common/common.smk"

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
fasta = config.get("genome",{}).get("fasta")
if not fasta:
    raise ValueError(
        "gatk_prepare module requires 'genome.fasta' in config. "
        "Please provide a valid reference genome FASTA path."
    )
input_bam_substring = config.get("input_bam_substring") or ""
rule gatk_index:
    input:
        fasta = fasta
    output:
        fai_index = outdir + "/index/genome.fa.fai",
        dict_index = outdir + "/index/genome.dict",
        fasta_link = outdir + "/index/genome.fa"
    log:
        logdir + "/index/gatk_index.log"
    threads: rule_res("gatk_index").get("threads", 4)
    resources:
        mem_mb = rule_res("gatk_index").get("mem_mb", 8000),
    conda: "gatk.yaml"
    container:
        sif("gatk.yaml", "gatk_index")
    params:
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="gatk_index", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start gatk_index at {current_time}")
            index_outdir = os.path.join(outdir, "index")
            os.makedirs(index_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"gatk_index_{current_time}.sh")
            plain_script = os.path.join(index_outdir, f"gatk_index.sh")
            cmd1 = ["ln", "-s", input.fasta, output.fasta_link]
            cmd2 = [
                params.gatk, "CreateSequenceDictionary",
                "--java-options", params.javaOptions,
                "-R", input.fasta,
                "-O", output.dict_index
            ]
            if params.tmp_dir:
                cmd2.extend(["--TMP_DIR", params.tmp_dir])
            cmd3 = [
                "samtools", "faidx",
                input.fasta,
                "-o", output.fai_index
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                f.write(f"echo 'gatk_index was completed successfully for {input.fasta} !'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during gatk_index execution: {str(e)}\n")
            raise RuntimeError(f"Error occurred while running gatk_index: {e}, you can check the log file {log_path} for more details.")

def get_input_for_addReadsGroup(wildcards):
    logger.debug(f"addReadsGroup called with {wildcards}")
    in_dict = {}
    if input_bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + input_bam_substring + ".bam")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
    return in_dict

rule addReadsGroup:
    input:
        unpack(get_input_for_addReadsGroup)
    output:
        bam = temp(outdir + "/{sample_id}/{sample_id}.addReadsGroup.bam"),
        bai = temp(outdir + "/{sample_id}/{sample_id}.addReadsGroup.bai")
    log:
        logdir + "/{sample_id}/addReadsGroup.log"
    conda: "gatk.yaml"
    container:
        sif("gatk.yaml", "addReadsGroup")
    threads: rule_res("addReadsGroup").get("threads", 8)
    resources:
        mem_mb = rule_res("addReadsGroup").get("mem_mb", 16000),
    params:
        id = "{sample_id}",
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        RGLB = config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGLB") or "lib1",
        RGPL = config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGPL") or "illumina",
        RGPU = config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGPU") or "unit1",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="addReadsGroup", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start addReadsGroup for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"addReadsGroup_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"addReadsGroup.sh")
            cmd1 = [
                params.gatk, "AddOrReplaceReadGroups", "--java-options", f'"{params.javaOptions}"',
                "--INPUT", input.bam,
                "--OUTPUT", output.bam,
                "-SO", "coordinate",
                "--RGLB", params.RGLB,
                "--RGPL", params.RGPL,
                "--RGPU", params.RGPU,
                "--RGSM", params.id
            ]
            if params.tmp_dir:
                cmd1.extend(["--TMP_DIR", params.tmp_dir])
            cmd2 = [
                params.samtools, "index",
                "-@", str(threads),
                output.bam,
                "-o", output.bai
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(f"echo 'addReadsGroup was completed successfully for {wildcards.sample_id} !'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during addReadsGroup execution for sample {wildcards.sample_id}: {str(e)}\n")
            raise RuntimeError(f"Error occurred while running addReadsGroup for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

rule MarkDuplicates:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.addReadsGroup.bam",
        bai = outdir + "/{sample_id}/{sample_id}.addReadsGroup.bai"
    output:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bai",
        metrics = outdir + "/{sample_id}/{sample_id}.Markdup-metrics.txt"
    log:
        logdir + "/{sample_id}/MarkDuplicates.log"
    conda: "gatk.yaml"
    container:
        sif("gatk.yaml", "MarkDuplicates")
    threads: rule_res("MarkDuplicates").get("threads", 8)
    resources:
        mem_mb = rule_res("MarkDuplicates").get("mem_mb", 32000),
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="MarkDuplicates", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start MarkDuplicates for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"MarkDuplicates_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"MarkDuplicates.sh")
            cmd = [
                params.gatk, "MarkDuplicates", "--java-options", f'"{params.javaOptions}"',
                "--INPUT", input.bam,
                "--OUTPUT", output.bam,
                "--CREATE_INDEX", "true",
                "--VALIDATION_STRINGENCY", "SILENT",
                "--METRICS_FILE", output.metrics
            ]
            if params.tmp_dir:
                cmd.extend(["--TMP_DIR", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'MarkDuplicates was completed successfully for {wildcards.sample_id} !'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during MarkDuplicates execution for sample {wildcards.sample_id}: {str(e)}\n")
            raise RuntimeError(f"Error occurred while running MarkDuplicates for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

_targets_bed = config.get("genome", {}).get("targets")
_fai_index = config.get("genome", {}).get("fai_index")
_dict_index = config.get("genome", {}).get("dict_index")


def get_input_for_collectHsMetrics(wildcards):
    """Resolve collectHsMetrics inputs; missing targets fails at DAG build.

    The capture region (BED or Picard interval_list) is mandatory: WGS runs
    have no bait design and simply do not request this rule's outputs, so a
    missing ``genome.targets`` here means a wiring error and must stop the
    DAG build with a self-identifying message.
    """
    if not _targets_bed:
        raise ValueError(
            "gatk_prepare/collectHsMetrics requires 'genome.targets' (capture "
            "BED or interval list) in the config JSON, then rerun."
        )
    in_dict = {
        "targets": _targets_bed,
        "bam": outdir + f"/{wildcards.sample_id}/{wildcards.sample_id}.sorted_markdup.bam",
        "bai": outdir + f"/{wildcards.sample_id}/{wildcards.sample_id}.sorted_markdup.bai",
    }
    if _fai_index and _dict_index:
        in_dict["fai"] = _fai_index
        in_dict["dict"] = _dict_index
        in_dict["fasta"] = fasta
    else:
        in_dict["fai"] = outdir + "/index/genome.fa.fai"
        in_dict["dict"] = outdir + "/index/genome.dict"
        in_dict["fasta"] = outdir + "/index/genome.fa"
    return in_dict


rule collectHsMetrics:
    """Per-sample hybrid-selection QC metrics (GATK CollectHsMetrics).

    Runs on the analysis-ready markdup BAM (the exact file the variant
    callers consume), so the reported on-target depth and coverage
    uniformity describe the effective callable territory. The targets BED is
    converted to a Picard interval_list with BedToIntervalList (bait =
    target: the region is reconstructed from actual coverage, there is no
    separate bait design); a user-supplied ``.interval_list`` is staged
    as-is. Both outputs are declared unconditionally so the DAG shape never
    depends on the input format.
    """
    input:
        unpack(get_input_for_collectHsMetrics),
    output:
        metrics = outdir + "/{sample_id}/{sample_id}.hs_metrics.txt",
        interval_list = outdir + "/{sample_id}/{sample_id}.targets.interval_list",
    log:
        logdir + "/{sample_id}/collectHsMetrics.log",
    conda:
        "gatk.yaml"
    container:
        sif("gatk.yaml", "collectHsMetrics")
    threads: rule_res("collectHsMetrics").get("threads", 4)
    resources:
        mem_mb = rule_res("collectHsMetrics").get("mem_mb", 16000),
    params:
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        javaOptions = config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx8g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="collectHsMetrics", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start collectHsMetrics for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"collectHsMetrics_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "collectHsMetrics.sh")
            cmds = []
            if str(input.targets).endswith(".interval_list"):
                cmds.append(["cp", str(input.targets), str(output.interval_list)])
            else:
                cmd_bed = [
                    params.gatk, "BedToIntervalList",
                    "--java-options", params.javaOptions,
                    "--INPUT", str(input.targets),
                    "--OUTPUT", str(output.interval_list),
                    "--SEQUENCE_DICTIONARY", str(input.dict),
                ]
                if params.tmp_dir:
                    cmd_bed.extend(["--TMP_DIR", params.tmp_dir])
                cmds.append(cmd_bed)
            cmd_hsm = [
                params.gatk, "CollectHsMetrics",
                "--java-options", params.javaOptions,
                "--INPUT", str(input.bam),
                "--OUTPUT", str(output.metrics),
                "--REFERENCE_SEQUENCE", str(input.fasta),
                "--BAIT_INTERVALS", str(output.interval_list),
                "--TARGET_INTERVALS", str(output.interval_list),
            ]
            if params.tmp_dir:
                cmd_hsm.extend(["--TMP_DIR", params.tmp_dir])
            cmds.append(cmd_hsm)
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                for cmd in cmds:
                    f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'collectHsMetrics was completed successfully for {wildcards.sample_id} !'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during collectHsMetrics execution for {wildcards.sample_id}: {str(e)}")
            raise RuntimeError(f"Error occurred while running collectHsMetrics for {wildcards.sample_id}: {e}, you can check the log file {log_path} for more information.")


rule gatk_prepare_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bai",
        metrics = outdir + "/{sample_id}/{sample_id}.Markdup-metrics.txt"
    resources:
        mem_mb = rule_res("gatk_prepare_result").get("mem_mb", 8000),

