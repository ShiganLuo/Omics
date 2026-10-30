include: "../common/common.smk"

import shutil

from snakemake.logging import logger
import os
import time

outdir = config.get("outdir", "results/censatbuild")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
genome_name = config.get("genome_name") or "genome"
p = config.get("Params", {}).get("censatbuild", {})

ARRAYS_OUT = outdir + f"/{genome_name}.sat_arrays.bed"
CENSAT_OUT = outdir + f"/{genome_name}.censat.bed"
RM_OUT = outdir + f"/{genome_name}/RepeatMasker/asm.bp.p_ctg.fa.out"


rule censatbuild_run:
    """Aggregate target: censat region-class BED + stats."""
    input:
        bed = CENSAT_OUT,
        stats = outdir + f"/{genome_name}.censat.stats.tsv",


# ---- RepeatMasker reuses modules/centromere (Dfam species library + FamDB) ----
# The centromere module consumes the hifiasm layout {outdir}/{sample_id}/
# assembly/asm.bp.p_ctg.fa; censatbuild_link_genome places the target genome
# there so repeatmasker_run is reused unmodified. FamDB (Dfam h5 partitions,
# e.g. dfam39_full.0.h5 root + dfam39_full.7.h5 Mammalia) is provided via a
# famdb.conf bind mount at runtime (see README "FamDB").
centromere_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": env,
    "indir": outdir,
    "outdir": outdir,
    "logdir": logdir,
    "samples": [genome_name],
    "Procedure": {},
    "Params": {"RepeatMasker": {"species": p.get("species", "Mus musculus")}},
}

module centromere:
    snakefile: "../centromere/centromere.smk"
    config: centromere_config

use rule repeatmasker_run from centromere as censatbuild_repeatmasker_run


rule censatbuild_link_genome:
    """Place the target genome at the layout expected by repeatmasker_run."""
    input:
        fasta = lambda wildcards: _require(genome_cfg.get("fasta"), "genome.fasta", "censatbuild"),
    output:
        link = outdir + f"/{genome_name}/assembly/asm.bp.p_ctg.fa",
    log:
        logdir + "/censatbuild_link_genome.log"
    threads: 1
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("censatbuild_link_genome", log_file=log_path)
        link_path = str(output.link)
        os.makedirs(os.path.dirname(link_path), exist_ok=True)
        target = os.path.abspath(str(input.fasta))
        if os.path.islink(link_path) or os.path.exists(link_path):
            os.remove(link_path)
        os.symlink(target, link_path)
        rule_logger.info(f"Linked genome {target} -> {link_path}")


rule censatbuild_arrays:
    """Cluster monomer hits into satellite arrays and classify them."""
    input:
        rm_out = RM_OUT,
    output:
        arrays = ARRAYS_OUT,
        hits = outdir + f"/{genome_name}.sat_hits.bed",
    log:
        logdir + "/censatbuild_arrays.log"
    threads: 1
    conda:
        "censatbuild.yaml"
    container:
        sif("censatbuild.yaml")
    params:
        script = ROOT_DIR + "/modules/censatbuild/bin/censatbuild_rm2bed.py",
        gap = p.get("gap", 10000),
        major_pattern = p.get("major_pattern", "major"),
        minor_pattern = p.get("minor_pattern", "minor"),
        class_pattern = p.get("class_pattern", "Satellite"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("censatbuild_arrays", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start satellite clustering for {genome_name} at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script = os.path.join(cmd_dir(logdir, "group"), f"censatbuild_arrays_{genome_name}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"censatbuild_arrays_{genome_name}.sh")
            cmd = [
                "python", params.script,
                "--rm-out", str(input.rm_out),
                "--arrays-out", str(output.arrays),
                "--hits-out", str(output.hits),
                "--gap", str(params.gap),
                "--major-pattern", str(params.major_pattern),
                "--minor-pattern", str(params.minor_pattern),
                "--class-pattern", str(params.class_pattern),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "censatbuild_arrays for {genome_name} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during censatbuild_arrays for {genome_name}: {e}\n")
            raise RuntimeError(f"Error occurred during censatbuild_arrays for {genome_name}: {e}\n")


rule censatbuild_flank:
    """Add flanking regions and write the final censat region-class BED."""
    input:
        arrays = ARRAYS_OUT,
        fai = lambda wildcards: _require(genome_cfg.get("fai"), "genome.fai", "censatbuild"),
    output:
        bed = CENSAT_OUT,
        stats = outdir + f"/{genome_name}.censat.stats.tsv",
    log:
        logdir + "/censatbuild_flank.log"
    threads: 1
    conda:
        "censatbuild.yaml"
    container:
        sif("censatbuild.yaml")
    params:
        flank = p.get("flank", 10000),
        flank_name = f"flank_{p.get('flank', 10000) // 1000}kb",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("censatbuild_flank", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start flank/merge for {genome_name} at {current_time}")
            script = os.path.join(cmd_dir(logdir, "group"), f"censatbuild_flank_{genome_name}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"censatbuild_flank_{genome_name}.sh")
            flank_bed = outdir + f"/{genome_name}.sat_flank.bed"
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(f"bedtools flank -i {input.arrays} -g {input.fai} -b {params.flank} | "
                        f"awk 'BEGIN{{OFS=\"\\t\"}}{{$4=\"{params.flank_name}\"; print}}' > {flank_bed}\n")
                f.write(f"cat {input.arrays} {flank_bed} | sort -k1,1 -k2,2n > {output.bed}\n")
                f.write(f"awk 'BEGIN{{OFS=\"\\t\"}} !/^#/ {{n[$4]++; s[$4]+=$3-$2}} "
                        f"END{{print \"class\",\"count\",\"total_bp\"; for (c in n) print c, n[c], s[c]}}' "
                        f"{output.bed} > {output.stats}\n")
                f.write(f'echo "censatbuild_flank for {genome_name} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during censatbuild_flank for {genome_name}: {e}\n")
            raise RuntimeError(f"Error occurred during censatbuild_flank for {genome_name}: {e}\n")


rule censatbuild_mcg_pileup:
    """Optional: CpG methylation pileup (bedGraph) from an aligned BAM.

    Consumed by modules/fiberseq/censat as genome.mcpg_pileup. Needs
    Params.censatbuild.bam (e.g. a fiberseq aligned BAM).
    """
    input:
        bam = lambda wildcards: _require(p.get("bam"), "Params.censatbuild.bam", "censatbuild"),
    output:
        bedgraph = outdir + f"/{genome_name}.mcpg_pileup.bedGraph",
    log:
        logdir + "/censatbuild_mcg_pileup.log"
    threads: p.get("threads", 32)
    conda:
        "censatbuild.yaml"
    container:
        sif("censatbuild.yaml")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("censatbuild_mcg_pileup", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start modkit CpG pileup at {current_time}")
            script = os.path.join(cmd_dir(logdir, "group"), f"censatbuild_mcg_pileup_{genome_name}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"censatbuild_mcg_pileup_{genome_name}.sh")
            tmp_bed = outdir + f"/{genome_name}.mcpg_pileup.bedMethyl"
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(f"modkit pileup --cpg --combine-strands {input.bam} -o {tmp_bed}\n")
                f.write(f"awk 'BEGIN{{OFS=\"\\t\"}} !/^#/ {{print $1,$2,$3,$11/100}}' "
                        f"{tmp_bed} > {output.bedgraph}\n")
                f.write(f'echo "censatbuild_mcg_pileup at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during censatbuild_mcg_pileup: {e}\n")
            raise RuntimeError(f"Error occurred during censatbuild_mcg_pileup: {e}\n")
