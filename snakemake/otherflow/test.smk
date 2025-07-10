shell.prefix("set -x; set -e;")
configfile: "config/RNASNP.yaml"
indir = config.get('indir', '../data')
outdir = config.get('outdir', '../output')
# 动态获取所有样本（包含单端和双端）
paired_samples = glob_wildcards(indir + "/fq/{sample_id}_1.fq.gz").sample_id
Allsamples = glob_wildcards(indir + "/fq/{sample_id}.fq.gz").sample_id
single_samples = [sample for sample in Allsamples if re.match(r"^SRR\d+$", sample)]
all_samples = paired_samples + single_samples
genomes=['mouse']
print(all_samples)

# 动态请求最终输出
def all_output_request():
    output = dict()
    output['RG'] = [expand(outdir + "/RG/{genome}/{sample_id}.bam",genome=genomes,sample_id=all_samples)]
    return list(output.values())

rule all:
    input: all_output_request()

rule trimming_Paired:
    input:
        fastq1=indir+"/fq/{sample_id}_1.fq.gz",
        fastq2=indir+"/fq/{sample_id}_2.fq.gz"
    output:
        fastq1=outdir+"/cutadapt/{sample_id}_1.fq.gz",
        fastq2=outdir+"/cutadapt/{sample_id}_2.fq.gz",
        report1=outdir+"/log/{sample_id}/trimming_statistics_1.txt",
        report2=outdir+"/log/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir=outdir+"/cutadapt",
        quality=30,
        trim_galore="/opt/TrimGalore-0.6.10/trim_galore"
    threads: 6
    log:
        log=outdir+"/log/{sample_id}/trimming.txt"
    shell:
        """
        # trim_galore can automatically judge the fq quality scoring system,it's no need to add such as --phred33 --phred64
        {params.trim_galore} --paired  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fq.gz_trimming_report.txt {output.report2}
        """
rule trimming_Single:
    input:
        fastq = indir+"/fq/{sample_id}.fq.gz"
    output:
        fastq = outdir+"/cutadapt/{sample_id}Single.fq.gz",
        report = outdir+"/log/{sample_id}/trimming_statistics.txt"
    params:
        outdir=outdir+"/cutadapt",
        quality=30,
        trim_galore="/opt/TrimGalore-0.6.10/trim_galore"
    threads: 6
    log:
        log=outdir+"/log/{sample_id}/trimming.txt"
    shell:
        """
        {params.trim_galore} --phred33  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_trimmed.fq.gz {output.fastq}
        mv {params.outdir}/{wildcards.sample_id}.fq.gz_trimming_report.txt {output.report}
        """

def get_alignment_input(wildcards):
    """动态判断输入文件类型"""
    # 构造可能的输入路径
    paired_r1 = f"{outdir}/cutadapt/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{outdir}/cutadapt/{wildcards.sample_id}_2.fq.gz"
    single = f"{outdir}/cutadapt/{wildcards.sample_id}Single.fq.gz"
    
    # 检查文件实际存在情况
    if wildcards.sample_id in paired_samples:
        print([paired_r1, paired_r2])
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        return [single]
    else:
        raise FileNotFoundError(
            f"Missing input files for sample {wildcards.sample_id}\n"
            f"Checked paths:\n- {paired_r1}\n- {paired_r2}\n- {single}"
        )

rule star_align:
    input:
        fastq = get_alignment_input,
        genome_index = lambda wildcards: config['STAR'][wildcards.genome]['genome_index']
    output:
        outfile = outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        log = outdir + "/log/{sample_id}/{genome}/star_align.log"
    threads: 25
    params:
        outPrefix = outdir + "/2pass/{sample_id}/{genome}/{sample_id}",
        STAR = config["STAR"]["procedure"],
        # 动态判断输入参数
        input_params = lambda wildcards, input: \
            f"{input[0]} {input[1]}" if len(input) == 2 else f"{input[0]}"
    shell:
        """
        mkdir -p $(dirname {params.outPrefix})
        {params.STAR} --runThreadN {threads} \
            --genomeDir {input.genome_index} \
            --twopassMode Basic \
            --readFilesCommand zcat \
            --readFilesIn {params.input_params} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NM \
            --outFileNamePrefix {params.outPrefix} > {log.log} 2>&1
        """

rule addReadsGroup:
    input:
        bam = rules.star_align.output.outfile
    output:
        bam = outdir + "/RG/{genome}/{sample_id}.bam",
        bai = outdir + "/RG/{genome}/{sample_id}.bam.bai"
    log:
        log = outdir + "/log/{genome}/{sample_id}/addReadsGroup.log"
    threads: 16
    conda:
        config['conda']['cfDNA_base']
    params:
        id = "{sample_id}",
        java = "--java-options -Xmx15G",
        RGLB = config["addReadsGroup"]["RGLB"],
        RGPL = config["addReadsGroup"]["RGPL"],
        RGPU = config["addReadsGroup"]["RGPU"]
    shell:
        """
        echo "Processing sample: {wildcards.sample_id}" > {log.log}
        gatk AddOrReplaceReadGroups {params.java} \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            -SO coordinate \
            --RGLB {params.RGLB} \
            --RGPL {params.RGPL} \
            --RGPU {params.RGPU} \
            --RGSM {params.id} >> {log.log} 2>&1
        
        echo "Indexing BAM file..." >> {log.log}
        samtools index -@ {threads} {output.bam} >> {log.log} 2>&1
        """
