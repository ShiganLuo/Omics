import logging
logger = logging.getLogger(__name__)
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
gtf = config.get("gtf","")
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
treated_ip_samples = config.get("treated_ip_samples", [])
treated_input_samples = config.get("treated_input_samples", [])

def get_input_for_diff_exomePeak():
    ip_bams = [indir + f"/{sample_id}.dedup.bam" for sample_id in ip_samples]
    input_bams = [indir + f"/{sample_id}.dedup.bam" for sample_id in input_samples]
    treated_ip_bams = [indir + f"/{sample_id}.dedup.bam" for sample_id in treated_ip_samples]
    treated_input_bams = [indir + f"/{sample_id}.dedup.bam" for sample_id in treated_input_samples]
    bam_dict = {
        "ip_bams": ip_bams,
        "input_bams": input_bams,
        "treated_ip_bams": treated_ip_bams,
        "treated_input_bams": treated_input_bams
    }
    logger.info(f"Input BAM files for diff_exomePeak: {bam_dict}")
    return bam_dict

rule diff_exomePeak:
    input:
        **get_input_for_diff_exomePeak()
    output:
        diff_peak_bed = outdir + "/diff_peaks_gene_names.bed",
        diff_peak_xls = outdir + "/diff_peaks_gene_names.xls",
        sig_siff_bed = outdir + "/sig_diff_peak_gene_names.bed",
        sig_siff_xls = outdir + "/sig_diff_peak_gene_names.xls",
        con_sig_diff_bed = outdir + "/con_sig_diff_peak_gene_names.bed"
        con_sig_diff_xls = outdir + "/con_sig_diff_peak_gene_names.xls"
    log:
        logdir + "/endpoint/exomePeak.log"
    conda:
        "exomePeak.yaml"
    threads: 1
    params:
        exomePeak = config.get('Procedure',{}).get('exomePeak') or 'exomePeak',
    shell:
        """
        Rscript bin/exomePeak.r \
            --gtf {gtf} \
            --ip_bams {input.ip_bams} \
            --input_bams {input.input_bams} \
            --treated_ip_bams {input.treated_ip_bams} \
            --treated_input_bams {input.treated_input_bams} \
            --outprefix {outdir} \
            > {log} 2>&1
        
        Rscript bin/geneId2name.r \
            --infile {outdir}/diff_peaks.bed \
            --gtf {gtf} \
            --outfile {outdir}/diff_peaks_gene_names.bed
        Rscript bin/geneId2name.r \
            --infile {outdir}/diff_peaks.xls \
            --gtf {gtf} \
            --outfile {outdir}/diff_peaks_gene_names.xls

        Rscript bin/geneId2name.r \
            --infile {outdir}/con_sig_diff_peak.bed \
            --gtf {gtf} \
            --outfile {outdir}/con_sig_diff_peak_gene_names.bed
        Rscript bin/geneId2name.r \
            --infile {outdir}/con_sig_diff_peak.xls \
            --gtf {gtf} \
            --outfile {outdir}/con_sig_diff_peak_gene_names.xls
        
        Rscript bin/geneId2name.r \
            --infile {outdir}/sig_diff_peak.bed \
            --gtf {gtf} \
            --outfile {outdir}/sig_diff_peak_gene_names.bed
        Rscript bin/geneId2name.r \
            --infile {outdir}/sig_diff_peak.xls \
            --gtf {gtf} \
            --outfile {outdir}/sig_diff_peak_gene_names.xls
        
        rm -f {outdir}/diff_peaks.bed {outdir}/diff_peaks.xls {outdir}/con_sig_diff_peak.bed {outdir}/con_sig_diff_peak.xls {outdir}/sig_diff_peak.bed {outdir}/sig_diff_peak.xls

        """

def get_input_for_call_exomePeak():
    ip_bams = [indir + f"/{sample_id}.bam" for sample_id in ip_samples]
    input_bams = [indir + f"/{sample_id}.bam" for sample_id in input_samples]
    bam_dict = {
        "ip_bams": ip_bams,
        "input_bams": input_bams
    }
    return bam_dict
rule call_exomePeak:
    input:
        get_input_for_call_exomePeak()
    output:
        all_peak_bed = outdir + "/all_peaks_gene_names.bed",
        all_peak_xls = outdir + "/all_peaks_gene_names.xls",
        con_peaks_bed = outdir + "/con_peaks_gene_names.bed",
        con_peaks_xls = outdir + "/con_peaks_gene_names.xls"
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
        --gtf {gtf} \
        --ip_bams {input.ip_bams} \
        --input_bams {input.input_bams} \
        --outprefix {outdir} \
        > {log} 2>&1

        Rscript bin/geneId2name.r \
            --infile {outdir}/all_peaks.bed \
            --gtf {gtf} \
            --outfile {outdir}/all_peaks_gene_names.bed
        Rscript bin/geneId2name.r \
            --infile {outdir}/all_peaks.xls \
            --gtf {gtf} \
            --outfile {outdir}/all_peaks_gene_names.xls
        
        Rscript bin/geneId2name.r \
            --infile {outdir}/con_peaks.bed \
            --gtf {gtf} \
            --outfile {outdir}/con_peaks_gene_names.bed
        Rscript bin/geneId2name.r \
            --infile {outdir}/con_peaks.xls \
            --gtf {gtf} \
            --outfile {outdir}/con_peaks_gene_names.xls
        rm -f {outdir}/all_peaks.bed {outdir}/all_peaks.xls {outdir}/con_peaks.bed {outdir}/con_peaks.xls
        """

rule exomePeak_result:
    input:
        call_exomePeak = rules.call_exomePeak.output,
        diff_exomePeak = rules.diff_exomePeak.output

