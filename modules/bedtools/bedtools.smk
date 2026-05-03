from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule iCLIP_bedtools:
    input:
        bam = indir + "/{sample_id}.dedup.bam",
        bai = indir + "/{sample_id}.dedup.bam.bai",
        chromosome_sizes = config.get('genome',{}).get('chrom_sizes')
    output:
        bed = outdir + "/{sample_id}/{sample_id}.bed",
        plus_bedgraph = outdir + "/{sample_id}/{sample_id}.plus.bedgraph",
        minus_bedgraph = outdir + "/{sample_id}/{sample_id}.minus.bedgraph"
    log:
        log = logdir + "/{sample_id}/bedtools.log"
    threads: 4
    conda:
        "bedtools.yaml"
    params:
        bedtools = config.get('Procedure',{}).get('bedtools') or 'bedtools'
    shell:
        """
        sample_name=$(basename {output.bed} .bed)
        outdir=$(dirname {output.bed})
        echo "Processing sample: ${{sample_name}}; output directory: ${{outdir}}" > {log}
        {params.bedtools} bamtobed -i {input.bam} > {output.bed} 2>> {log}

        {params.bedtools} shift -m 1 -p -1 -i {output.bed} -g {input.chromosome_sizes} > ${{outdir}}/${{sample_name}}.shifted.bed 2>> {log}
        {params.bedtools} genomecov -bg -strand + -5 -scale 1000000 -i ${{outdir}}/${{sample_name}}.shifted.bed -g {input.chromosome_sizes} > ${{outdir}}/${{sample_name}}.plus.bedgraph 2>> {log}
        {params.bedtools} genomecov -bg -strand - -5 -scale 1000000 -i ${{outdir}}/${{sample_name}}.shifted.bed -g {input.chromosome_sizes} > ${{outdir}}/${{sample_name}}.minus.bedgraph 2>> {log}
        """
