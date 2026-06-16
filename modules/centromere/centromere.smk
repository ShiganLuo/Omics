from snakemake.logging import logger
import os
import time
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

rm_species = config.get("Params", {}).get("RepeatMasker", {}).get("species", "mouse")
bam_substring = config.get("bam_substring") or ""

def get_input_for_hifiasm(wildcards):
    logger.info(f"hifiasm_assemble called with {wildcards}")
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
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start hifiasm assembly for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"{wildcards.sample_id}/hifiasm_{current_time}.sh")
        fq_gz = os.path.join(outdir, f"{wildcards.sample_id}/assembly/{wildcards.sample_id}.hifi.fq.gz")
        cmd0 = [
            "samtools", "fastq", "-@", str(threads), input.bam, "|", "gzip", "-c", ">", fq_gz
        ]
        cmd1 = [
            "hifiasm", "-o", params.prefix, "-t", str(threads), fq_gz
        ]
        cmd2 = [
            "awk", "'/^S/{print \">\"$2; print $3}'", f"{params.prefix}.bp.p_ctg.gfa", ">", output.fasta
        ]
        cmd3 = ["rm", "-f", fq_gz]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd0) + "\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(" ".join(cmd3) + "\n")
        shell(f"bash {script} > {log} 2>&1")



rule repeatmasker_run:
    """Run RepeatMasker on assembled contigs to annotate satellite DNA."""
    input:
        fasta = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.fa"
    output:
        out = outdir + "/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.out",
        tbl = outdir + "/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.tbl"
    log:
        logdir + "/{sample_id}/repeatmasker.log"
    threads: 12
    conda:
        "centromere.yaml"
    params:
        species = rm_species,
        rm_dir = outdir + "/{sample_id}/RepeatMasker",
        RepeatMasker = config.get("Procedure", {}).get("RepeatMasker") or "RepeatMasker"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start RepeatMasker for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"{wildcards.sample_id}/repeatmasker_{current_time}.sh")
        cmd = [
            params.RepeatMasker,
            "-species", params.species,
            "-pa", str(threads),
            "-dir", params.rm_dir,
            "-gff",
            input.fasta
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"mkdir -p {params.rm_dir}\n")
            f.write(" ".join(cmd) + "\n")
        shell(f"bash {script} > {log} 2>&1")


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
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start centromere stats extraction for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"{wildcards.sample_id}/centromere_extract_{current_time}.sh")
        cmd = [
            "python", params.script,
            "--rm_out", input.out,
            "--output", output.stats
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell(f"bash {script} > {log} 2>&1")



rule centromere_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        stats = outdir + "/{sample_id}/{sample_id}.centromere_stats.txt"
