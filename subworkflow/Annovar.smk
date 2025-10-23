rule TEcoutCPM:
    input:
       infile = outdir + "/counts/humanTEcount.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCPM.cntTable"
    log:
        log = outdir + "/log/human/TEcoutCPM.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = "scripts/SNP/run-NormCountMat.R"
    shell:
        """
        /usr/bin/Rscript {params.script} -i {input.infile} -o {output.outfile} > {log.log} 2>&1
        """

rule commonExpression:
    input:
        infile = outdir + "/counts/humanTEcountCPM.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCommon.cntTable"
    log:
        log = outdir + "/log/human/commonExpression.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = "scripts/SNP/commonExpression.py"
    shell:
        """
            python {params.script} --input {input.infile} --output {output.outfile} --threshold 5 > {log.log} 2>&1
        """

rule getBed:
    input:
        infile = outdir + "/counts/humanTEcountCommon.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCommon.bed"
    log:
        log = outdir + "/log/human/getBed.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = "scripts/SNP/getBed.py",
        gtf = config['getBed']['human']['gtf'],
        TE_gtf = config['getBed']['human']['TE_gtf']
    shell:
        """
            python {params.script} \
                --input {input.infile} \
                --output {output.outfile} \
                --Gtf {params.gtf} \
                --TEGtf {params.TE_gtf} > {log.log} 2>&1
        """
rule vcfIntersectBed:
    input:
        vcf = outdir + "/filter/vcf/human/{sample_id}.vcf.gz",
        bed = outdir + "/counts/humanTEcountCommon.bed"
    output:
        outfile = outdir + "/filter/vcf/human/{sample_id}Common.vcf"
    log:
        log = outdir + "/log/human/{sample_id}/vcfIntersectBed.log"
    conda:
        config['conda']['RNA-SNP']
    threads:4 #防止同时执行太多，爆内存
    shell:
        """
            bedtools intersect -a {input.vcf} -b {input.bed} -wa -wb > {output.outfile} 2>{log.log}
        """
rule annovar_convert:
    input:
        vcf = outdir + "/filter/vcf/human/{sample_id}Common.vcf"
    output:
        avinput = outdir + "/annovar/human/{sample_id}/{sample_id}.avinput"
    log:
        log = outdir + "/log/human/{sample_id}/annovar_convert.log"
    params:
        convert = "/opt/annovar/convert2annovar.pl"
    threads:4 #防止同时执行太多，爆内存
    shell:
        """
        /usr/bin/perl {params.convert} \
        -format vcf4 \
        -withfreq {input.vcf} > {output.avinput} 2>{log.log}
        """
rule annovar_table:
    input:
        avinput = outdir + "/annovar/human/{sample_id}/{sample_id}.avinput"
    output:
        outfile = outdir + "/annovar/human/{sample_id}/{sample_id}.GRCh38_multianno.csv"
    log:
        log = outdir + "/log/human/{sample_id}/annovar_table.log"
    params:
        db = config['annovar']['human']['db'],
        buildver = config['annovar']['human']['buildver'],
        # annotate = "/opt/annovar/annotate_variation.pl",
        table = "/opt/annovar/table_annovar.pl",
        out = outdir + "/annovar/human/{sample_id}/{sample_id}"
    threads:4 #防止同时执行太多，爆内存
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