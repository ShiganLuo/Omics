
import os

def get_hisat2_index(wildcards):
    logger.info(f"[get_hisat2_index] called with wildcards: {wildcards}")
    config_index_prefix = config.get('genome',{}).get(wildcards.genome,{}).get('hisat2_index_prefx') or None
    if config_index_prefix:
        first_file = f"{config_index_prefix}.1.ht2"
        if os.path.exists(first_file):
            return [f"{config_index_prefix}.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]
    return expand(
        outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}.{{idx}}.ht2",
        idx = [1, 2, 3, 4, 5, 6, 7, 8]
    )

rule TEtranscript_prepare_hisat2:
    input:
        fastq = get_alignment_input,
        index = get_hisat2_index 
    output:
        outfile = temp(outdir + "/TEtranscripts/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam")
    log:
        outdir + "/log/TEtranscripts/{sample_id}/{genome}/hisat2_align.log"
    threads: 12
    conda:
        config['conda']['run']
    params:
        HISAT2 = config.get('Procedure',{}).get('hisat2') or 'hisat2',
        SAMTOOLS = config.get('Procedure',{}).get('samtools') or 'samtools',
        # 这里需要更健壮地提取前缀，兼容不同的后缀情况
        index_prefix = lambda wildcards, input: input.index[0].rsplit('.', 2)[0],
        input_params = lambda wildcards, input: \
            f"-1 {input.fastq[0]} -2 {input.fastq[1]}" if len(input.fastq) == 2 else f"-U {input.fastq[0]}"
    shell:
        """
        mkdir -p $(dirname {output.outfile})
        {params.HISAT2} -x {params.index_prefix} \
            {params.input_params} \
            --no-mixed \
            --no-discordant \
            -k 100 \
            --score-min L,0,-0.2 \
            -p {threads} 2> {log} | \
        {params.SAMTOOLS} sort -@ {threads} -o {output.outfile}
        """