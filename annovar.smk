indir = config.get('indir', '../data')  # 如果没有传递 indir，则使用 'data' 作为默认值
outdir = config.get('outdir', '../output')  # 如果没有传递 outdir，则使用 'output' 作为默认值
samples = glob_wildcards(outdir + "/vcf/{sample_id}.bed").sample_id
def request():
    output = dict()
    output['annovar'] = expand(outdir + "/annovar/{sample_id}/{sample_id}.GRCm39_multianno.csv",sample_id=samples)

    return list(output.values())


rule all:
    input:request()
rule commonExpression:
    input:
        infile = outdir + "/counts/{sample_id}/{genome}/TEcount.cntTable"
    output:
        outfile = outdir + "/counts/{sample_id}/{genome}/TEcountCommon.cntTable"
    log:
        log = outdir + "/log/{genome}/commonExpression.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = "scripts/SNP/commonExpression.py"
    shell:
    """
        python {params.script} --input {input.infile} --output {output.outfile} > {log.log} 2>&1
    """
rule getBed:
    input:
        infile = outdir + "/counts/{sample_id}/{genome}/TEcountCommon.cntTable"
    output:
        outfile = outdir + "/counts/{sample_id}/{genome}/TEcountCommon.bed"
    log:
        log = outdir + "/log/{genome}/getBed.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = "scripts/SNP/getBed.py",
        exon_gtf = config['getBed']['exon_gtf'],
        TE_gtf = config['getBed']['TE_gtf']
    shell:
    """
        python {params.script} \
            --input {input.infile} \
            --output {output.outfile} \
            --exonGtf {params.exon_gtf} \
            --TE_gtf {params.TE_gtf} > {log.log} 2>&1
    """
rule vcfIntersectBed:
    input:
        vcf = outdir + "/filter/vcf/{genome}/{sample_id}.vcf.gz",
        bed = outdir + "/counts/{sample_id}/{genome}/TEcountCommon.bed"
    output:
        outfile = outdir + "/filter/vcf/{genome}/{sample_id}Common.vcf",
    log:
        log = outdir + "/log/{genome}/{sample_id}/vcfIntersectBed.log"
    conda:
        config['conda']['RNA-SNP']
    shell:
    """
        bedtools intersect -a {input.vcf} -b {input.bed} -wa -wb > {output.outfile} 2>{log.log}
    """
rule annovar_convert:
    input:
        vcf = outdir + "/filter/vcf/{genome}/{sample_id}Common.vcf"
    output:
        avinput = outdir + "/annovar/{genome}/{sample_id}/{sample_id}.avinput"
    log:
        log = outdir + "/log/{genome}/{sample_id}/annovar_convert.log"
    params:
        convert = "/opt/annovar/convert2annovar.pl"
    shell:
        """
        /usr/bin/perl {params.convert} \
        -format vcf4 \
        -withfreq {input.vcf} > {output.avinput} 2>{log.log}
        """
rule annovar_table:
    input:
        avinput = outdir + "/annovar/{genome}/{sample_id}/{sample_id}.avinput"
    output:
        outfile = outdir + "/annovar/{genome}/{sample_id}/{sample_id}.GRCm39_multianno.csv"
    log:
        log = outdir + "/log/{genome}/{sample_id}/annovar_table.log"
    params:
        db = "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/annovar/GRCm39/",
        buildver = "GRCm39",
        # annotate = "/opt/annovar/annotate_variation.pl",
        table = "/opt/annovar/table_annovar.pl",
        out = outdir + "/annovar/{sample_id}/{sample_id}"
    shell:
        """
        /usr/bin/perl {params.table} \
        {input.avinput} {params.db} \
        -buildver {params.buildver} \
        -out {params.out} \
        -remove -protocol refGene \
        -operation g \
        -nastring . \
        -csvout > {log.log} 2>&1
        """