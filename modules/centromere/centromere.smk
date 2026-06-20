"""Centromere analysis module.

This module performs:
1. Assembly with hifiasm
2. RepeatMasker for satellite DNA annotation
3. Centromere statistics extraction

Note: This module includes common.smk for shared utilities.
"""

# Include common utilities
include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""


def get_input_for_hifiasm(wildcards):
    """Get input files for hifiasm assembly."""
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}." + bam_substring + ".bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}." + bam_substring + ".bam.pbi")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}.bam.pbi")
    return in_dict


rule hifiasm_assemble:
    """Assemble HiFi reads with hifiasm."""
    input:
        unpack(get_input_for_hifiasm)
    output:
        fasta = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.fa"
    log:
        logdir + "/{sample_id}/hifiasm.log"
    threads: 12
    conda:
        "centromere.yaml"
    params:
        prefix = outdir + "/{sample_id}/assembly/asm"
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="hifiasm_assemble", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start hifiasm assembly for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/hifiasm_{current_time}.sh")
            fq_gz = os.path.join(outdir, f"{wildcards.sample_id}/assembly/{wildcards.sample_id}.hifi.fq.gz")
            cmd0 = [
                "samtools", "fastq", "-@", str(threads), input.bam, "|", "gzip", "-c", ">", fq_gz
            ]
            cmd1 = [
                "hifiasm", "-o", params.prefix, "-t", str(threads), fq_gz
            ]
            cmd2 = [
                "awk", "'/^S/{print \">$2; print $3}'", f"{params.prefix}.bp.p_ctg.gfa", ">", output.fasta
            ]
            cmd3 = ["rm", "-f", fq_gz]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd0) + "\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"hifiasm assembly failed for sample {wildcards.sample_id} with error: {e}\n")
            raise f"hifiasm assembly failed for sample {wildcards.sample_id} with error: {e}"


rule repeatmasker_init:
    """Build RepeatMasker library cache once to avoid concurrent makeblastdb conflicts."""
    input:
        dfam_h5 = config.get("Params", {}).get("RepeatMasker", {}).get("dfam_h5") or ""
    output:
        outdir + "/.repeatmasker_lib_init.done"
    log:
        logdir + "/all/RepeatMasker/repeatmasker_init.log"
    conda:
        "centromere.yaml"
    params:
        species = config.get("Params", {}).get("RepeatMasker", {}).get("species", "Mus musculus"),
        RepeatMasker = config.get("Procedure", {}).get("RepeatMasker") or "RepeatMasker"
    threads: 1
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="repeatmasker_init", log_file=log[0])
            logger.info("Building RepeatMasker library cache (one-time init)")
            tmpdir = tempfile.mkdtemp()
            dummy_fa = os.path.join(tmpdir, "dummy.fa")
            with open(dummy_fa, "w") as f:
                f.write(">dummy\nATCGATCGATCG\n")
            cmd1 = [
                params.RepeatMasker,
                "-species", f"'{params.species}'",
                "-pa", str(threads),
                "-dir", tmpdir,
                "-gff",
                dummy_fa
            ]
            cmd3 = ["touch", output[0]]
            cmd2 = ["rm", "-rf", tmpdir]
            cmd4 = ["echo", "RepeatMasker library cache initialized."]
            script = os.path.join(tmpdir, "init.sh")
            if input.dfam_h5 and os.path.exists(input.dfam_h5):
                logger.warning("RepeatMasker don't support specifying custom Dfam h5 directly, pipeline will try to find the corresponding library dir and link your dfam h5 to the library dir with a name that RepeatMasker can recognize. Please ensure the Dfam h5 you provided is compatible with your RepeatMasker version and includes relevant satellite DNA for target species.")
                RepeatMasker_dir = shell("which RepeatMasker", read=True)# Get RepeatMasker executable path, enviroment only play in shell
                if RepeatMasker_dir is None:
                    logger.error(f"RepeatMasker executable not found: {params.RepeatMasker}")
                    raise FileNotFoundError(f"RepeatMasker executable not found: {params.RepeatMasker}")
                RepeatMasker_lib_dir = os.path.join(os.path.dirname(RepeatMasker_dir), "../share/RepeatMasker/Libraries/famdb")
                if not os.path.exists(RepeatMasker_lib_dir):
                    logger.error(f"RepeatMasker library directory not found: {RepeatMasker_lib_dir}")
                    raise FileNotFoundError(f"RepeatMasker library directory not found: {RepeatMasker_lib_dir}")
                target_h5_path = os.path.join(RepeatMasker_lib_dir, os.path.basename(input.dfam_h5))
                if os.path.exists(target_h5_path):
                    logger.info(f"Custom Dfam h5 already exists in RepeatMasker library directory: {target_h5_path}, skipping linking.")
                else:
                    os.symlink(input.dfam_h5, target_h5_path)
                    logger.info(f"Linked custom Dfam h5 to RepeatMasker library directory: {target_h5_path}")
            else:
                logger.warning("No valid Dfam h5 provided for RepeatMasker, using default libraries. please ensure the default libraries include relevant satellite DNA for target species.")
            with open(script, "w") as f:
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                f.write(" ".join(cmd4) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"RepeatMasker init failed with error: {e}\n")
            raise f"RepeatMasker init failed with error: {e}"


def get_input_for_repeatmasker(wildcards):
    """Get input files for RepeatMasker."""
    in_dict = {}
    in_dict["fasta"] = os.path.join(outdir, f"{wildcards.sample_id}/assembly/asm.bp.p_ctg.fa")
    dfam_h5 = config.get("Params", {}).get("RepeatMasker", {}).get("dfam_h5")
    if dfam_h5 and os.path.exists(dfam_h5):
        in_dict["lib_init"] = os.path.join(outdir, ".repeatmasker_lib_init.done")
    return in_dict

rule repeatmasker_run:
    """Run RepeatMasker on assembled contigs to annotate satellite DNA."""
    input:
        unpack(get_input_for_repeatmasker)
    output:
        out = outdir + "/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.out",
        tbl = outdir + "/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.tbl"
    log:
        logdir + "/{sample_id}/repeatmasker.log"
    threads: 12
    conda:
        "centromere.yaml"
    params:
        species = config.get("Params", {}).get("RepeatMasker", {}).get("species", "Mus musculus"),
        rm_dir = outdir + "/{sample_id}/RepeatMasker",
        RepeatMasker = config.get("Procedure", {}).get("RepeatMasker") or "RepeatMasker",
        
    run:
        # setup_logger is now available from common.smk
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="repeatmasker_run", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start RepeatMasker for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/repeatmasker_{current_time}.sh")
            cmd = [
                params.RepeatMasker,
                "-species", f"'{params.species}'",
                "-pa", str(threads),
                "-dir", params.rm_dir,
                "-gff",
                input.fasta
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"RepeatMasker failed for sample {wildcards.sample_id} with error: {e}\n")
            raise f"RepeatMasker failed for sample {wildcards.sample_id} with error: {e}"


rule centromere_extract:
    """Extract satellite DNA statistics from RepeatMasker output."""
    input:
        out = outdir + "/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.out"
    output:
        stats = outdir + "/{sample_id}/{sample_id}.centromere_stats.txt"
    log:
        logdir + "/{sample_id}/centromere_extract.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/centromere/bin/extract_centromere_stats.py"),
        python = config.get("Procedure", {}).get("python") or "python"
    conda:
        "centromere.yaml"
    threads: 1
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="centromere_extract", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start centromere stats extraction for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/centromere_extract_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--rm_out", input.out,
                "--output", output.stats
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Centromere stats extraction failed for sample {wildcards.sample_id} with error: {e}\n")
            raise f"Centromere stats extraction failed for sample {wildcards.sample_id} with error: {e}"


rule centromere_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        stats = outdir + "/{sample_id}/{sample_id}.centromere_stats.txt"
