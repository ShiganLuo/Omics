include: "../common/common.smk"
from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
treated_ip_samples = config.get("treated_ip_samples", [])
treated_input_samples = config.get("treated_input_samples", [])

def get_input_for_diff_exomePeak():
    ip_bams = [indir + f"/{sample_id}/{sample_id}.dedup.bam" for sample_id in ip_samples]
    input_bams = [indir + f"/{sample_id}/{sample_id}.dedup.bam" for sample_id in input_samples]
    treated_ip_bams = [indir + f"/{sample_id}/{sample_id}.dedup.bam" for sample_id in treated_ip_samples]
    treated_input_bams = [indir + f"/{sample_id}/{sample_id}.dedup.bam" for sample_id in treated_input_samples]
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
        con_sig_diff_bed = outdir + "/con_sig_diff_peak_gene_names.bed",
        con_sig_diff_xls = outdir + "/con_sig_diff_peak_gene_names.xls"
    log:
        logdir + "/endpoint/exomePeak.log"
    conda:
        "exomePeak.yaml"
    threads: 1
    params:
        gtf = config.get("gtf",""),
        exomePeak_script = ROOT_DIR + "/modules/exomePeak/bin/exomePeak.r",
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("diff_exomePeak", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start diff_exomePeak at {current_time}")
            sample_outdir = outdir
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"diff_exomePeak_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"Rscript {params.exomePeak_script} \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --ip_bams {input.ip_bams} \\\n")
                f.write(f"    --input_bams {input.input_bams} \\\n")
                f.write(f"    --treated_ip_bams {input.treated_ip_bams} \\\n")
                f.write(f"    --treated_input_bams {input.treated_input_bams} \\\n")
                f.write(f"    --outprefix {outdir} \\\n")
                f.write(f"    > {log} 2>&1\n")
                f.write(f"\n")
                f.write(f"Rscript bin/geneId2name.r \\\n")
                f.write(f"    --infile {outdir}/diff_peaks.bed \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/diff_peaks_gene_names.bed\n")
                f.write(f"Rscript bin/geneId2name.r \\\n")
                f.write(f"    --infile {outdir}/diff_peaks.xls \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/diff_peaks_gene_names.xls\n")
                f.write(f"\n")
                f.write(f"Rscript bin/geneId2name.r \\\n")
                f.write(f"    --infile {outdir}/con_sig_diff_peak.bed \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/con_sig_diff_peak_gene_names.bed\n")
                f.write(f"Rscript bin/geneId2name.r \\\n")
                f.write(f"    --infile {outdir}/con_sig_diff_peak.xls \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/con_sig_diff_peak_gene_names.xls\n")
                f.write(f"\n")
                f.write(f"Rscript bin/geneId2name.r \\\n")
                f.write(f"    --infile {outdir}/sig_diff_peak.bed \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/sig_diff_peak_gene_names.bed\n")
                f.write(f"Rscript bin/geneId2name.r \\\n")
                f.write(f"    --infile {outdir}/sig_diff_peak.xls \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/sig_diff_peak_gene_names.xls\n")
                f.write(f"\n")
                f.write(f"rm -f {outdir}/diff_peaks.bed {outdir}/diff_peaks.xls {outdir}/con_sig_diff_peak.bed {outdir}/con_sig_diff_peak.xls {outdir}/sig_diff_peak.bed {outdir}/sig_diff_peak.xls\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during diff_exomePeak: {e}\n")
            logger.error(f"Error occurred during diff_exomePeak: {e}")
            raise e

def get_input_for_call_exomePeak():
    ip_bams = [indir + f"/{sample_id}/{sample_id}.bam" for sample_id in ip_samples]
    input_bams = [indir + f"/{sample_id}/{sample_id}.bam" for sample_id in input_samples]
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
        exomePeak_script = ROOT_DIR + "/modules/exomePeak/bin/exomePeak.r",
        geneId2name_script = ROOT_DIR + "/modules/exomePeak/bin/geneId2name.r",
        gtf = config.get("gtf","")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("call_exomePeak", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start call_exomePeak at {current_time}")
            sample_outdir = outdir
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"call_exomePeak_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"Rscript {params.exomePeak_script} \\\n")
                f.write(f"--gtf {params.gtf} \\\n")
                f.write(f"--ip_bams {input.ip_bams} \\\n")
                f.write(f"--input_bams {input.input_bams} \\\n")
                f.write(f"--outprefix {outdir} \\\n")
                f.write(f"> {log} 2>&1\n")
                f.write(f"\n")
                f.write(f"Rscript {params.geneId2name_script} \\\n")
                f.write(f"    --infile {outdir}/all_peaks.bed \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/all_peaks_gene_names.bed\n")
                f.write(f"Rscript {params.geneId2name_script} \\\n")
                f.write(f"    --infile {outdir}/all_peaks.xls \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/all_peaks_gene_names.xls\n")
                f.write(f"\n")
                f.write(f"Rscript {params.geneId2name_script} \\\n")
                f.write(f"    --infile {outdir}/con_peaks.bed \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/con_peaks_gene_names.bed\n")
                f.write(f"Rscript {params.geneId2name_script} \\\n")
                f.write(f"    --infile {outdir}/con_peaks.xls \\\n")
                f.write(f"    --gtf {params.gtf} \\\n")
                f.write(f"    --outfile {outdir}/con_peaks_gene_names.xls\n")
                f.write(f"rm -f {outdir}/all_peaks.bed {outdir}/all_peaks.xls {outdir}/con_peaks.bed {outdir}/con_peaks.xls\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during call_exomePeak: {e}\n")
            logger.error(f"Error occurred during call_exomePeak: {e}")
            raise e

rule exomePeak_result:
    input:
        call_exomePeak = rules.call_exomePeak.output,
        diff_exomePeak = rules.diff_exomePeak.output

