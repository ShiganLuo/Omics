include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"
BEDTOOLS = config.get("Procedure", {}).get("bedtools") or "bedtools"
smallrna_bed = config.get("genome", {}).get("smallrna_bed")
smallrna_fasta = config.get("genome", {}).get("smallrna_fasta")
pass1_outdir = config.get("pass1_outdir", "")

p3g = config.get("Params", {}).get("star_3pass_gene", {})
hard_clip_length = p3g.get("hard_clip_length", 10)

# ── Extract smallRNA reads from pass1 BAM and split into per-gene FASTQ ──────
# Uses bedtools intersect to find reads overlapping each smallRNA gene,
# then samtools fastq to produce per-gene FASTQ files.
rule star_3pg_extract_per_gene:
    input:
        bam = pass1_outdir + "/{sample_id}/{sample_id}.bam",
        bai = pass1_outdir + "/{sample_id}/{sample_id}.bam.bai",
        bed = smallrna_bed,
        fasta = smallrna_fasta,
    output:
        done = outdir + "/per_gene_fq/{sample_id}/.done",
    log:
        logdir + "/star3pg/{sample_id}/extract_per_gene.log"
    threads: 4
    conda:
        "star_3pass.yaml"
    run:
        import gzip
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("star_3pg_extract_per_gene", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3pg_extract_per_gene for sample {wildcards.sample_id} at {current_time}")

            outdir_fq = os.path.dirname(str(output.done))
            os.makedirs(outdir_fq, exist_ok=True)

            # Read BED to get gene list
            genes = []
            with open(input.bed) as f:
                for line in f:
                    cols = line.strip().split("\t")
                    if len(cols) >= 7:
                        genes.append({
                            "chrom": cols[0],
                            "start": cols[1],
                            "end": cols[2],
                            "gene_id": cols[3],
                            "strand": cols[5],
                            "gene_name": cols[6],
                        })

            rule_logger.info(f"Found {len(genes)} smallRNA genes in BED")

            script = os.path.join(outdir_fq, f"extract_per_gene_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                for g in genes:
                    gene_id = g["gene_id"]
                    gene_name = g["gene_name"]
                    # Safe filename: gene_id may contain dots
                    safe_id = gene_id.replace(".", "_")
                    gene_fq = os.path.join(outdir_fq, f"{wildcards.sample_id}_{safe_id}.fq.gz")

                    # Extract reads overlapping this gene, then convert to FASTQ
                    f.write(f"echo 'Processing gene: {gene_name} ({gene_id})'\n")
                    f.write(
                        f"{BEDTOOLS} intersect -abam {input.bam} -b <(printf '%s\\t%s\\t%s\\t%s\\t0\\t%s\\n' "
                        f"'{g['chrom']}' '{g['start']}' '{g['end']}' '{g['gene_id']}' '{g['strand']}') -u "
                        f"| {SAMTOOLS} fastq -0 /dev/null -s /dev/null - "
                        f"| gzip -c > {gene_fq} 2>/dev/null || true\n"
                    )

                f.write(f"touch {output.done}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"star_3pg_extract_per_gene completed for sample {wildcards.sample_id}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during star_3pg_extract_per_gene for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during star_3pg_extract_per_gene for sample {wildcards.sample_id}: {e}")
            raise e


# ── Per-gene local alignment: re-align each gene's FASTQ to its single-gene ──
# genomic sequence using STAR local alignment (three-pass style).
# Requires per-gene STAR indexes built from smallrna_fasta.
rule star_3pg_align_per_gene:
    input:
        done = outdir + "/per_gene_fq/{sample_id}/.done",
        index_dir = config.get("genome", {}).get("smallrna_star_index_dir"),
    output:
        done = outdir + "/per_gene_bam/{sample_id}/.done",
    log:
        logdir + "/star3pg/{sample_id}/align_per_gene.log"
    threads: 4
    conda:
        "star_3pass.yaml"
    params:
        STAR = config.get("Procedure", {}).get("STAR") or "STAR",
        smallrna_fasta = smallrna_fasta,
        per_gene_fq_dir = outdir + "/per_gene_fq/{sample_id}",
        per_gene_bam_dir = outdir + "/per_gene_bam/{sample_id}",
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("star_3pg_align_per_gene", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3pg_align_per_gene for sample {wildcards.sample_id} at {current_time}")

            bam_dir = params.per_gene_bam_dir.format(sample_id=wildcards.sample_id)
            fq_dir = params.per_gene_fq_dir.format(sample_id=wildcards.sample_id)
            os.makedirs(bam_dir, exist_ok=True)

            # Collect per-gene FASTQ files
            fq_files = [f for f in os.listdir(fq_dir) if f.endswith(".fq.gz") and f.startswith(f"{wildcards.sample_id}_")]
            rule_logger.info(f"Found {len(fq_files)} per-gene FASTQ files")

            script = os.path.join(bam_dir, f"align_per_gene_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                for fq_name in fq_files:
                    # gene_id is the part after {sample_id}_ and before .fq.gz
                    gene_key = fq_name[len(wildcards.sample_id) + 1:].replace(".fq.gz", "")
                    # Reverse the safe_id transformation to get original gene_id
                    # The FASTA headers in smallrna_fasta use gene_id::gene_name format
                    # STAR index was built from this FASTA
                    fq_path = os.path.join(fq_dir, fq_name)
                    out_prefix = os.path.join(bam_dir, f"{gene_key}.")
                    bam_out = os.path.join(bam_dir, f"{gene_key}.bam")

                    f.write(f"echo 'Aligning gene: {gene_key}'\n")
                    f.write(
                        f"if [ -s {fq_path} ]; then\n"
                        f"  {params.STAR} --runThreadN {threads} "
                        f"--genomeDir {input.index_dir} "
                        f"--readFilesIn {fq_path} "
                        f"--readFilesCommand zcat "
                        f"--alignEndsType Local "
                        f"--outSAMtype BAM Unsorted "
                        f"--outFileNamePrefix {out_prefix} "
                        f"--outReadsUnmapped None "
                        f"--outStd Log\n"
                        f"  if [ -f {out_prefix}Aligned.out.bam ]; then\n"
                        f"    mv {out_prefix}Aligned.out.bam {bam_out}\n"
                        f"    {SAMTOOLS} index {bam_out}\n"
                        f"  fi\n"
                        f"fi\n"
                    )
                f.write(f"touch {output.done}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"star_3pg_align_per_gene completed for sample {wildcards.sample_id}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during star_3pg_align_per_gene for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during star_3pg_align_per_gene for sample {wildcards.sample_id}: {e}")
            raise e


# ── Merge per-gene BAMs into a single sample BAM ─────────────────────────────
rule star_3pg_merge:
    input:
        done = outdir + "/per_gene_bam/{sample_id}/.done",
    output:
        bam = outdir + "/{sample_id}/{sample_id}.bam",
        bai = outdir + "/{sample_id}/{sample_id}.bam.bai",
    log:
        logdir + "/star3pg/{sample_id}/merge.log"
    threads: 4
    conda:
        "star_3pass.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("star_3pg_merge", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3pg_merge for sample {wildcards.sample_id} at {current_time}")

            bam_dir = os.path.dirname(str(output.bam))
            os.makedirs(bam_dir, exist_ok=True)

            per_gene_bam_dir = os.path.join(os.path.dirname(bam_dir), "per_gene_bam", wildcards.sample_id)
            gene_bams = []
            if os.path.isdir(per_gene_bam_dir):
                for f in sorted(os.listdir(per_gene_bam_dir)):
                    if f.endswith(".bam"):
                        gene_bams.append(os.path.join(per_gene_bam_dir, f))

            rule_logger.info(f"Found {len(gene_bams)} per-gene BAMs to merge")

            script = os.path.join(bam_dir, f"merge_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                if gene_bams:
                    # samtools merge with -f (force overwrite) and -b (input list)
                    bam_list = os.path.join(bam_dir, f"bam_list_{current_time}.txt")
                    with open(bam_list, "w") as bl:
                        for gb in gene_bams:
                            bl.write(f"{gb}\n")
                    f.write(f"{SAMTOOLS} merge -@ {threads} -f -b {bam_list} {output.bam}\n")
                    f.write(f"{SAMTOOLS} index -@ {threads} {output.bam}\n")
                else:
                    # No per-gene BAMs: create empty BAM
                    f.write(f"{SAMTOOLS} view -bS /dev/null > {output.bam}\n")
                    f.write(f"{SAMTOOLS} index {output.bam}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"star_3pg_merge completed for sample {wildcards.sample_id}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during star_3pg_merge for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during star_3pg_merge for sample {wildcards.sample_id}: {e}")
            raise e
