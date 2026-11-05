from snakemake.logging import logger
include: "../../common/common.smk"

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
fai_index = config.get("genome", {}).get("fai_index")
dict_index = config.get("genome", {}).get("dict_index")
fasta = config.get("genome", {}).get("fasta")
targets = config.get("genome", {}).get("targets")
known_sites = config.get("genome", {}).get("known_sites")
interval = config.get("genome", {}).get("interval")
# Germline suppression resources for Mutect2/FilterMutectCalls (both optional):
#   genome.germline_resource - population allele-frequency VCF (e.g. gnomAD),
#                       assembly-tied so it lives in genome.references.<name>
#   Params.gatk.Mutect2.panel_of_normals - sites observed in a SELF-BUILT
#                       normal panel (cohort/lab property, not a genome
#                       property); recurrent artifacts get hard-filtered
germline_resource = config.get("genome", {}).get("germline_resource")
panel_of_normals = config.get("Params", {}).get("gatk", {}).get("Mutect2", {}).get("panel_of_normals")
# Strain-aware germline filtering: full-genotype source VCF (e.g. MGP REL-2112
# with one sample column per strain) plus the per-pair strain name resolved by
# node.py from the meta ``strain`` column. When both are set, each pair gets a
# strain-specific germline resource (variants carried by THAT strain, AF from
# its genotype) built by rule germline_resource_strain; without them the
# all-strain germline_resource is used as-is.
strain_resource_src = config.get("genome", {}).get("strain_resource")
pair_strains = config.get("pair_strains") or {}  # {experimental_sample_id: strain}

def get_input_somaticMutect2(wildcards):
    """
    This function determines the appropriate input BAM file for the HaplotypeCaller step based on the presence of known_sites and interval parameters.
    If both known_sites and interval are provided, it assumes that BQSR has been performed and uses the BQSR-corrected BAM file. Otherwise, it falls back to using the MarkDuplicates BAM file. 
    It also checks for the presence of the corresponding index files and includes them as input if they exist.
    """
    in_dict = {}

    if known_sites and interval:
        logger.debug(f"Using known_sites: {known_sites} and interval: {interval}")
        include: "../gatk_bqsr/gatk_bqsr.smk"
        in_dict["normal_bam"] = f"{indir}/{wildcards.normal_sample_id}/{wildcards.normal_sample_id}.sorted_markdup.bqsr.bam"
        in_dict["experimental_bam"] = f"{indir}/{wildcards.experimental_sample_id}/{wildcards.experimental_sample_id}.sorted_markdup.bqsr.bam"
    else:
        logger.debug("No known_sites or interval specified in config, proceeding without them.")
        in_dict["normal_bam"] = f"{indir}/{wildcards.normal_sample_id}/{wildcards.normal_sample_id}.sorted_markdup.bam"
        in_dict["experimental_bam"] = f"{indir}/{wildcards.experimental_sample_id}/{wildcards.experimental_sample_id}.sorted_markdup.bam"

    if fai_index and dict_index:
        logger.debug(f"Using provided fai_index: {fai_index} and dict_index: {dict_index}")
        in_dict["fai"] = fai_index
        in_dict["dict"] = dict_index
        in_dict["fasta"] = fasta
    else:
        logger.debug("No fai_index or dict_index specified in config, using rule to generate them.")
        in_dict["fai"] = f"{indir}/index/genome.fa.fai"
        in_dict["dict"] = f"{indir}/index/genome.dict"
        in_dict["fasta"] = f"{indir}/index/genome.fa"
    return in_dict

rule somaticMutect2:
    input:
       unpack(get_input_somaticMutect2),
       targets = lambda wildcards: [targets] if targets else [],
       germline_resource = lambda wildcards: [germline_resource] if germline_resource else [],
       panel_of_normals = lambda wildcards: [panel_of_normals] if panel_of_normals else [],
    output:
        vcf = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/{normal_sample_id}_vs_{experimental_sample_id}.vcf.gz"
    log:
        logdir + "/{normal_sample_id}_{experimental_sample_id}_mutect2.log"
    conda:
        "../gatk.yaml"
    container:
        sif("../gatk.yaml", "somaticMutect2")
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = config.get("tmp_dir") or None
    threads: rule_res("somaticMutect2").get("threads", 10)
    resources:
        mem_mb = rule_res("somaticMutect2").get("mem_mb", 32000),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="somaticMutect2", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"somaticMutect2_{wildcards.normal_sample_id}_{current_time}.sh")
            sample_outdir = os.path.dirname(output.vcf)
            plain_script = os.path.join(sample_outdir, "somaticMutect2.sh")
            cmd = [
                params.gatk, "Mutect2",
                "--java-options", params.javaOptions,
                "-R", input.fasta,
                "-I", input.normal_bam,
                "-I", input.experimental_bam,
                "-normal", wildcards.normal_sample_id,
                "-O", output.vcf,
                "--native-pair-hmm-threads", str(threads)
            ]
            if input.targets:
                cmd.extend(["-L", str(input.targets[0])])
            # Germline suppression: population AF resource + normal panel.
            # Both flow through to FilterMutectCalls via the raw VCF header.
            if input.germline_resource:
                cmd.extend(["--germline-resource", str(input.germline_resource[0])])
            if input.panel_of_normals:
                cmd.extend(["--pon", str(input.panel_of_normals[0])])
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'somaticMutect2 completed successfully for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during somaticMutect2 execution: {str(e)}")
            raise RuntimeError(f"Error occurred while running somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id}: {e}, you can check the log file {log_path} for more details.")



def get_strain_for_pair(wildcards):
    """Resolve the strain name for a pair's tumor (experimental) sample."""
    strain = pair_strains.get(wildcards.experimental_sample_id)
    if not strain:
        raise ValueError(
            f"gatk_somatic: no strain in pair_strains for sample "
            f"'{wildcards.experimental_sample_id}'; fill the meta strain column "
            "or unset genome.strain_resource to use the all-strain resource."
        )
    return strain


rule germline_resource_strain:
    """Build a strain-specific germline resource (variants carried by one strain).

    Filters the full-genotype source VCF (one sample column per strain, e.g.
    MGP REL-2112) down to sites where the given strain carries an ALT allele;
    AF becomes that strain's allele dosage (0.5 het / 1.0 hom). Sites-only
    output, bgzip+tabix, ready as --germline-resource.
    """
    input:
        src = lambda wildcards: _require(strain_resource_src, "genome.strain_resource", "gatk_somatic"),
    output:
        vcf = outdir + "/germline_strain/{strain}.vcf.gz",
    log:
        logdir + "/gatk/germline_strain_{strain}.log",
    threads: rule_res("germline_resource_strain").get("threads", 2)
    resources:
        mem_mb = rule_res("germline_resource_strain").get("mem_mb", 16000),
    conda:
        "../gatk.yaml"
    container:
        sif("../gatk.yaml", "germline_resource_strain")
    params:
        strain = lambda wildcards: wildcards.strain,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("germline_resource_strain", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start germline resource build for strain {wildcards.strain} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"germline_resource_strain_{wildcards.strain}_{current_time}.sh")
            sample_outdir = os.path.dirname(str(output.vcf))
            os.makedirs(sample_outdir, exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"germline_resource_strain_{wildcards.strain}.sh")
            # bcftools view -s <strain>: keep only the strain's sample column;
            # -i 'GT[0]="alt"' keeps sites where THAT strain carries an ALT
            # allele (GT[0] = the single remaining sample after -s; the bare
            # GT="alt" form silently no-ops on some bcftools versions);
            # annotate -x FORMAT strips genotypes -> sites-only af resource.
            # The filter expression contains quotes/brackets -> MUST shlex.quote
            # it (bare " ".join would let bash eat the quotes: skill pitfall).
            import shlex
            strain_expr = shlex.quote('GT[0]="alt"')
            # -G in the final view drops ALL genotype columns (sites-only);
            # annotate -x FORMAT alone leaves the sample column behind (verified).
            cmd = [
                "bcftools", "view", "-s", params.strain,
                "-i", strain_expr, "-Ou", str(input.src), "|",
                "bcftools", "annotate", "-x", "FORMAT", "-Ou", "|",
                "bcftools", "view", "-G", "-Oz", "-o", str(output.vcf),
                "&&", "bcftools", "index", "-t", "-f", str(output.vcf),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "germline resource for strain {wildcards.strain} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during germline resource build for strain {wildcards.strain}: {e}")
            raise RuntimeError(f"Error occurred during germline resource build for strain {wildcards.strain}: {e}, you can check the log file {log_path} for more details.")


rule filterMutectCalls:
    """Apply GATK FilterMutectCalls to a raw Mutect2 VCF (per pair).

    Turns the unfiltered somatic VCF into the filtered one: germline (via the
    germline resource recorded in the raw header), panel of normals, orientation
    bias and statistical filters land in the FILTER column. When strain-aware
    filtering is configured (genome.strain_resource + pair strain), the
    strain-specific resource is re-attached with --germline-resource so the
    germline prior matches the strain actually being analyzed.
    """
    input:
        raw_vcf = lambda wildcards: f"{outdir}/{{normal_sample_id}}_vs_{{experimental_sample_id}}/{{normal_sample_id}}_vs_{{experimental_sample_id}}.vcf.gz".replace("{{normal_sample_id}}", wildcards.normal_sample_id).replace("{{experimental_sample_id}}", wildcards.experimental_sample_id),
        fai = lambda wildcards: fai_index if fai_index else f"{indir}/index/genome.fa.fai",
        dict = lambda wildcards: dict_index if dict_index else f"{indir}/index/genome.dict",
        fasta = lambda wildcards: fasta if fasta else f"{indir}/index/genome.fa",
        germline_resource = lambda wildcards: (
            [outdir + "/germline_strain/" + pair_strains[wildcards.experimental_sample_id] + ".vcf.gz"]
            if (strain_resource_src and pair_strains.get(wildcards.experimental_sample_id)) else
            ([germline_resource] if germline_resource else [])
        ),
    output:
        vcf = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/{normal_sample_id}_vs_{experimental_sample_id}.filtered.vcf.gz",
    log:
        logdir + "/{normal_sample_id}_{experimental_sample_id}_filterMutectCalls.log",
    threads: rule_res("filterMutectCalls").get("threads", 2)
    resources:
        mem_mb = rule_res("filterMutectCalls").get("mem_mb", 8000),
    conda:
        "../gatk.yaml"
    container:
        sif("../gatk.yaml", "filterMutectCalls")
    params:
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        javaOptions = config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx8g",
        tmp_dir = config.get("tmp_dir") or None,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("filterMutectCalls", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start FilterMutectCalls for {wildcards.normal_sample_id} vs {wildcards.experimental_sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"filterMutectCalls_{wildcards.experimental_sample_id}_{current_time}.sh")
            sample_outdir = os.path.dirname(str(output.vcf))
            plain_script = os.path.join(sample_outdir, f"filterMutectCalls_{wildcards.experimental_sample_id}.sh")
            cmd = [
                params.gatk, "FilterMutectCalls",
                "--java-options", params.javaOptions,
                "-R", str(input.fasta),
                "-V", str(input.raw_vcf),
                "-O", str(output.vcf),
            ]
            if input.germline_resource:
                cmd += ["--germline-resource", str(input.germline_resource[0])]
            if params.tmp_dir:
                cmd += ["--tmp-dir", params.tmp_dir]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "filterMutectCalls completed successfully for {wildcards.normal_sample_id} vs {wildcards.experimental_sample_id}"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during filterMutectCalls for {wildcards.experimental_sample_id}: {e}")
            raise RuntimeError(f"Error occurred during filterMutectCalls for {wildcards.experimental_sample_id}: {e}, you can check the log file {log_path} for more details.")


_design_pairs = config.get("design_pairs") or []
rule gatk_somatic_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        filtered = [
            f"{outdir}/{p['control']}_vs_{p['experimental']}/{p['control']}_vs_{p['experimental']}.filtered.vcf.gz"
            for p in _design_pairs
        ],
