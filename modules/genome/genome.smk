include: "../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
smallrna_types = config.get("Params", {}).get("smallrna_types",
    ["miRNA", "snRNA", "snoRNA", "rRNA", "misc_RNA", "scRNA", "scaRNA", "vaultRNA"])
flank = config.get("Params", {}).get("smallrna_flank", 50)
ROOT_DIR = config.get("ROOT_DIR", os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

def get_input_for_chromosome_sizes(wildcards):
    fasta = config.get("genome", {}).get("fasta")
    if not fasta or not os.path.exists(fasta):
        raise ValueError(f"FASTA file for chromosome_sizes not found. Please check config['genome']['fasta']: {fasta}")
    return fasta
rule chromosome_sizes:
    input:
        fasta = get_input_for_chromosome_sizes
    output:
        chrom_sizes = outdir + "/genome/chrom.sizes"
    log:
        logdir + "/genome/chromosome_sizes.log"
    threads: 1
    conda:
        "genome.yaml"
    container:
        sif("genome.yaml")
    params:
        samtools = config.get('Procedure', {}).get('samtools') or 'samtools'
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("chromosome_sizes", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start chromosome_sizes at {current_time}")
            script = os.path.join(os.path.dirname(output.chrom_sizes), f"chromosome_sizes_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.samtools} faidx {input.fasta} 2>> {log_path}\n")
                f.write(f"cut -f1,2 {input.fasta}.fai > {output.chrom_sizes}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"chromosome_sizes completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"chromosome_sizes failed: {e}\n")
            raise e

# ── Extract smallRNA BED + FASTA from GENCODE GTF ───────────────────────
def get_input_for_extract_smallrna(wildcards):
    gtf = config.get("genome", {}).get("gtf")
    fasta = config.get("genome", {}).get("fasta")
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"GTF file for extract_smallrna not found. Please check config['genome']['gtf']: {gtf}")
    if not fasta or not os.path.exists(fasta):
        raise ValueError(f"FASTA file for extract_smallrna not found. Please check config['genome']['fasta']: {fasta}")
    in_dict = {
        "gtf": gtf,
        "fasta": fasta,
        "chrom_sizes": outdir + "/genome/chrom.sizes"
    }
    return gtf, fasta
rule extract_smallrna:
    input:
        unpack(get_input_for_extract_smallrna)
    output:
        bed = outdir + "/genome/smallrna/smallrna_genes.bed",
        fasta = outdir + "/genome/smallrna/smallrna_genes_flank.fa",
    log:
        logdir + "/genome/extract_smallrna.log"
    threads: 1
    conda:
        "genome.yaml"
    container:
        sif("genome.yaml")
    params:
        script = os.path.join(ROOT_DIR, "modules", "genome", "bin", "extract_smallrna.py"),
        outdir = outdir + "/genome/smallrna",
        flank = flank,
        types = " ".join(smallrna_types),
        bedtools = config.get("Procedure", {}).get("bedtools") or "bedtools",
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("extract_smallrna", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start extract_smallrna at {current_time}")
            script = os.path.join(params.outdir, f"extract_smallrna_{current_time}.sh")
            os.makedirs(params.outdir, exist_ok=True)
            with open(script, "w") as f:
                f.write(f"python {params.script} \\\n")
                f.write(f"    --gtf {input.gtf} \\\n")
                f.write(f"    --fasta {input.fasta} \\\n")
                f.write(f"    --chrom-sizes {input.chrom_sizes} \\\n")
                f.write(f"    --outdir {params.outdir} \\\n")
                f.write(f"    --flank {params.flank} \\\n")
                f.write(f"    --types {params.types} \\\n")
                f.write(f"    --bedtools {params.bedtools} \\\n")
                f.write(f"    > {log_path} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"extract_smallrna completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"extract_smallrna failed: {e}\n")
            raise e

def get_input_for_genome_index(wildcards):

    fasta = config.get("genome", {}).get("fasta")
    if not fasta or not os.path.exists(fasta):
        raise ValueError(f"FASTA file for genome_index not found. Please check config['genome']['fasta']: {fasta}")
    return fasta

rule genome_index:
    input:
        fasta = get_input_for_genome_index
    output:
        fai = config.get("genome", {}).get("fasta") + ".fai"
    log:
        logdir + "/genome/genome_index.log"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools",
    threads: 1
    run:
        log_path = str(log)
        try:
            open(log_path,"w").close
            rule_logger = setup_logger("genome_index", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start genome_index at {current_time}")
            script = os.path.join(params.outdir, f"genome_index{current_time}.sh")
            cmd = [
                params.samtools, "faidx", input.fasta
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'genone_index for {input.fasta} was successfuly build at {current_time}' ")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"genome_index failed: {e}\n")
            raise e