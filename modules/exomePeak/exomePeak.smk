import logging
logger = logging.getLogger(__name__)
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
treated_ip_samples = config.get("treated_ip_samples", [])
treated_input_samples = config.get("treated_input_samples", [])

def get_input_for_diff_exomePeak(wildcards):
    logger.info(f"Getting input for diff_exomePeak with wildcards: {wildcards}")
    ip_bams = []
    input_bams = []
    treated_ip_bams = []
    treated_input_bams = []
    if wildcards.sample_id in ip_samples:
        ip_bams.append(indir + f"/{wildcards.sample_id}.bam")
    if wildcards.sample_id in input_samples:
        input_bams.append(indir + f"/{wildcards.sample_id}.bam")
    if wildcards.sample_id in treated_ip_samples:
        treated_ip_bams.append(indir + f"/{wildcards.sample_id}.bam")
    if wildcards.sample_id in treated_input_samples:
        treated_input_bams.append(indir + f"/{wildcards.sample_id}.bam")
    if len(ip_bams) == 0 or len(input_bams) == 0 or len(treated_ip_bams) == 0 or len(treated_input_bams) == 0:
        logger.error(f"diff exomePeak is missing BAM files for sample_id: {wildcards.sample_id},must have at least one BAM file for each category.")
        logger.error(f"IP BAMs: {ip_bams}")
        logger.error(f"Input BAMs: {input_bams}")
        logger.error(f"Treated IP BAMs: {treated_ip_bams}")
        logger.error(f"Treated Input BAMs: {treated_input_bams}")
        raise ValueError(f"diff exomePeak is missing BAM files for sample_id: {wildcards.sample_id},must have at least one BAM file for each category.")
    logger.info(f"IP BAMs: {ip_bams}")
    logger.info(f"Input BAMs: {input_bams}")
    logger.info(f"Treated IP BAMs: {treated_ip_bams}")
    logger.info(f"Treated Input BAMs: {treated_input_bams}")
    return ip_bams, input_bams, treated_ip_bams, treated_input_bams

rule diff_exomePeak:
    input:
        get_input_for_exomePeak
    output:
        diff_peak = outdir + "/exomePeak_diff_peaks.tsv",
        sig_siff_peak = outdir + "/exomePeak_sig_siff_peaks.tsv",
        con_sig_diff_peak = outdir + "/exomePeak_con_sig_diff_peaks.tsv"
    log:
        logdir + "/endpoint/exomePeak.log"
    conda:
        "exomePeak.yaml"
    threads: 12
    params:
        exomePeak = config.get('Procedure',{}).get('exomePeak') or 'exomePeak'
    shell:
        """
        Rscript bin/exomePeak.r \
        --ip_bams {input[0]} \
        --input_bams {input[1]} \
        --treated_ip_bams {input[2]} \
        --treated_input_bams {input[3]} \
        --outprefix {outdir}/exomePeak \
        > {log} 2>&1
        """

def get_input_for_call_exomePeak(wildcards):
    ip_bams = []
    input_bams = []
    if wildcards.sample_id in ip_samples:
        ip_bams.append(indir + f"/{wildcards.sample_id}.bam")
    if wildcards.sample_id in input_samples:
        input_bams.append(indir + f"/{wildcards.sample_id}.bam")
    if len(ip_bams) == 0 or len(input_bams) == 0:
        logger.error(f"call exomePeak is missing BAM files for sample_id: {wildcards.sample_id},must have at least one BAM file for IP and one BAM file for Input.")
        logger.error(f"IP BAMs: {ip_bams}")
        logger.error(f"Input BAMs: {input_bams}")
        raise ValueError(f"call exomePeak is missing BAM files for sample_id: {wildcards.sample_id},must have at least one BAM file for IP and one BAM file for Input.")
    logger.info(f"IP BAMs: {ip_bams}")
    logger.info(f"Input BAMs: {input_bams}")
    return ip_bams, input_bams
rule call_exomePeak:
    input:
        get_input_for_call_exomePeak
    output:
        consistent_peak = outdir + "/exomePeak_consistent_peaks.tsv"
        all_peak = outdir + "/exomePeak_all_peaks.tsv"
    log:
        logdir + "/endpoint/call_exomePeak.log"
    conda:
        "exomePeak.yaml"
    threads: 12
    params:
        exomePeak = config.get('Procedure',{}).get('exomePeak') or 'exomePeak'
    shell:
        """
        Rscript bin/exomePeak.r \
        --ip_bams {input[0]} \
        --input_bams {input[1]} \
        --outprefix {outdir}/exomePeak \
        > {log} 2>&1
        """
rule exomePeak_result:
    input:
        diff_peak = outdir + "/exomePeak_diff_peaks.tsv",
        sig_siff_peak = outdir + "/exomePeak_sig_siff_peaks.tsv",
        con_sig_diff_peak = outdir + "/exomePeak_con_sig_diff_peaks.tsv",
        consistent_peak = outdir + "/exomePeak_consistent_peaks.tsv",
        all_peak = outdir + "/exomePeak_all_peaks.tsv"

