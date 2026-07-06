"""Telomere analysis module.

This module provides multiple approaches for telomere length measurement:

1. telogator2_run: Per-chromosome-arm telomere length (seed-extend matching)
2. assembly_telomere_scan: Scan assembly contig ends for telomeric repeats (Approach A)
3. read_density_telomere: Genome-wide average from read-level k-mer density (Approach B)
4. tidk_scan: Assembly-based telomere scan using tidk community tool (Approach C)

Note: For mouse telomeres (30-150kb), Approach A is recommended as HiFi reads
(~15-25kb) cannot span the full telomere for telogator2.
"""

include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""
species = config.get("Params", {}).get("telogator2", {}).get("species", "human")
# Directory containing hifiasm assembly (from centromere module)
assembly_dir = config.get("assembly_dir", "")


def get_telogator2_ref():
    """Get species-specific subtelomere reference path."""
    if species == "human":
        return ""  # use default
    # Find telogator2 installation and use non-human ref
    try:
        result = subprocess.run(["which", "telogator2"], capture_output=True, text=True)
        telogator2_bin = result.stdout.strip()
        # Resources are in site-packages/source/resources/
        resource_dir = os.path.join(os.path.dirname(telogator2_bin), "..",
            "lib", f"python{'.'.join(map(str, __import__('sys').version_info[:2]))}",
            "site-packages", "source", "resources")
        ref_file = os.path.join(resource_dir, "non-human", f"telogator-ref-{species}.fa.gz")
        if os.path.exists(ref_file):
            return ref_file
        logger.warning(f"No telogator2 reference found for species '{species}', using default (human)")
    except Exception as e:
        logger.warning(f"Error finding telogator2 ref: {e}")
    return ""


def get_input_for_telogator2(wildcards):
    logger.info(f"telogator2_run called with {wildcards}")
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bai")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bai")
    return in_dict


def get_input_for_read_density(wildcards):
    """Get BAM input for read-density telomere analysis."""
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bam")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
    return in_dict


def get_assembly_fasta(wildcards):
    """Get assembly FASTA path from centromere module output."""
    if assembly_dir:
        return os.path.join(assembly_dir, f"{wildcards.sample_id}/assembly/asm.bp.p_ctg.fa")
    return os.path.join(outdir, f"../centromere/{wildcards.sample_id}/assembly/asm.bp.p_ctg.fa")


def get_tidk_scan_input(wildcards):
    """Get input dict for tidk_scan.

    If tidk database already exists (~/.local/share/tidk/tidk_database.csv),
    skip tidk_init dependency. Otherwise include it.
    """
    in_dict = {"fasta": get_assembly_fasta(wildcards)}
    tidk_db = os.path.expanduser("~/.local/share/tidk/tidk_database.csv")
    if not os.path.exists(tidk_db):
        in_dict["db_init"] = outdir + "/.tidk_build.done"
    return in_dict


# ============================================================
# Rule 1: Telogator2 (existing, per-chromosome-arm)
# ============================================================
rule telogator2_run:
    """Run Telogator2 telomere length analysis on a single sample."""
    input:
        unpack(get_input_for_telogator2)
    output:
        allele_tsv = os.path.join(outdir, "{sample_id}/telogator2/tlens_by_allele.tsv"),
        allele_plot = os.path.join(outdir, "{sample_id}/telogator2/all_final_alleles.png"),
        atl_plot = os.path.join(outdir, "{sample_id}/telogator2/violin_atl.png")
    log:
        logdir + "/{sample_id}/telogator2.log"
    params:
        telogator2 = config.get("Procedure", {}).get("telogator2") or "telogator2",
        ref = get_telogator2_ref()
    threads: 16
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="telogator2_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start telogator2 for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/telogator2_{current_time}.sh")
            output_dir = os.path.join(outdir, f"{wildcards.sample_id}/telogator2")
            cmd = [
                params.telogator2, "-i", input.bam,
                "-o", output_dir,
                "-r", "hifi",
                "-p", str(threads)
            ]
            if params.ref:
                cmd.extend(["-t", params.ref])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"telogator2 failed for sample {wildcards.sample_id} with error: {e}\n")
            raise


# ============================================================
# Rule 2: Assembly-based telomere scan (Approach A)
# ============================================================
rule assembly_telomere_scan:
    """Scan hifiasm assembly contig ends for telomeric repeats.

    This is the recommended approach for mouse telomeres (30-150kb).
    Requires assembly from centromere module.
    """
    input:
        fasta = get_assembly_fasta
    output:
        tsv = os.path.join(outdir, "{sample_id}/assembly_scan/{sample_id}_assembly_telomere.tsv"),
        stats = os.path.join(outdir, "{sample_id}/assembly_scan/{sample_id}_assembly_telomere_stats.txt")
    log:
        logdir + "/{sample_id}/assembly_telomere_scan.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/bin/scan_assembly_telomere.py"),
        output_dir = os.path.join(outdir, "{sample_id}/assembly_scan"),
        scan_length = config.get("Params", {}).get("telogator2", {}).get("assembly_scan_length", 50000)
    threads: 1
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="assembly_telomere_scan", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start assembly telomere scan for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/assembly_telomere_scan_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--fasta", input.fasta,
                "--sample_name", wildcards.sample_id,
                "--output_dir", params.output_dir,
                "--scan_length", str(params.scan_length)
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Assembly telomere scan failed for sample {wildcards.sample_id} with error: {e}\n")
            raise


# ============================================================
# Rule 3: Read-level k-mer density (Approach B)
# ============================================================
rule read_density_telomere:
    """Estimate genome-wide average telomere length from read-level k-mer density.

    Counts telomeric k-mers across all HiFi reads and divides by chromosome arms.
    Gives a genome-wide average, not per-chromosome measurements.
    """
    input:
        unpack(get_input_for_read_density)
    output:
        tsv = os.path.join(outdir, "{sample_id}/read_density/{sample_id}_read_telomere.tsv"),
        stats = os.path.join(outdir, "{sample_id}/read_density/{sample_id}_read_telomere_stats.txt")
    log:
        logdir + "/{sample_id}/read_density_telomere.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/bin/read_density_telomere.py"),
        output_dir = os.path.join(outdir, "{sample_id}/read_density"),
        n_chrom_arms = config.get("Params", {}).get("telogator2", {}).get("n_chrom_arms", 40)
    threads: 1
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="read_density_telomere", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start read density telomere for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/read_density_telomere_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--bam", input.bam,
                "--sample_name", wildcards.sample_id,
                "--output_dir", params.output_dir,
                "--n_chrom_arms", str(params.n_chrom_arms)
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Read density telomere failed for sample {wildcards.sample_id} with error: {e}\n")
            raise


# ============================================================
# Rule 4: tidk scan (Approach C)
# ============================================================
rule tidk_init:
    """Build tidk reference database (one-time init)."""
    output:
        outdir + "/.tidk_build.done"
    log:
        logdir + "/all/tidk_init.log"
    params:
        tidk = config.get("Procedure", {}).get("tidk") or "tidk"
    threads: 1
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="tidk_init", log_file=log_path)
            logger.info("Building tidk reference database (one-time init)")
            script = os.path.join(outdir, "tidk_init.sh")
            cmd = [params.tidk, "build"]
            cmd_touch = ["touch", output[0]]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(" ".join(cmd_touch) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"tidk build failed with error: {e}\n")
            raise


rule tidk_scan:
    """Scan assembly for telomeric repeats using tidk community tool.

    tidk is a standardized tool for telomere identification from assemblies.
    """
    input:
        unpack(get_tidk_scan_input)
    output:
        tsv = os.path.join(outdir, "{sample_id}/tidk/{sample_id}_tidk_telomeres.tsv")
    log:
        logdir + "/{sample_id}/tidk_scan.log"
    params:
        output_dir = os.path.join(outdir, "{sample_id}/tidk"),
        tidk = config.get("Procedure", {}).get("tidk") or "tidk",
        motif = "TTAGGG"
    threads: 1
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="tidk_scan", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start tidk scan for sample {wildcards.sample_id} at {current_time}")

            script = os.path.join(outdir, f"{wildcards.sample_id}/tidk_scan_{current_time}.sh")
            # tidk search: --string, --output (prefix), --dir, positional FASTA
            cmd_search = [
                params.tidk, "search",
                "--string", params.motif,
                "--output", f"{wildcards.sample_id}_tidk_search",
                "--dir", params.output_dir,
                input.fasta
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("mkdir -p " + params.output_dir + "\n")
                f.write(" ".join(cmd_search) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            # tidk search outputs <output>_telomeric_repeat.tsv in --dir
            search_tsv = os.path.join(params.output_dir, f"{wildcards.sample_id}_tidk_search_telomeric_repeat.tsv")
            if os.path.exists(search_tsv):
                shutil.copy2(search_tsv, output.tsv)
                logger.info(f"Copied tidk search output to {output.tsv}")
            else:
                logger.warning(f"tidk search output not found: {search_tsv}")
                with open(output.tsv, "w") as f:
                    f.write("# No telomeric repeats found by tidk\n")

        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"tidk scan failed for sample {wildcards.sample_id} with error: {e}\n")
            raise


# ============================================================
# Result rules (aggregation points for subworkflow)
# ============================================================
rule telogator2_result:
    """Result aggregation for telogator2."""
    input:
        tsv = outdir + "/{sample_id}/telogator2/tlens_by_allele.tsv"


rule assembly_telomere_result:
    """Result aggregation for assembly-based telomere scan."""
    input:
        stats = outdir + "/{sample_id}/assembly_scan/{sample_id}_assembly_telomere_stats.txt"


rule read_density_telomere_result:
    """Result aggregation for read-density telomere estimation."""
    input:
        stats = outdir + "/{sample_id}/read_density/{sample_id}_read_telomere_stats.txt"


rule tidk_result:
    """Result aggregation for tidk scan."""
    input:
        tsv = outdir + "/{sample_id}/tidk/{sample_id}_tidk_telomeres.tsv"
