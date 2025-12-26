rule bowtie2_index_small:
    """
    fasta < 4GB
    """
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta']
    output:
        bw2_1   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.1.bt2",
        bw2_2   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.2.bt2",
        bw2_3   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.3.bt2",
        bw2_4   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.4.bt2",
        bw2_rv1 = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.rev.1.bt2",
        bw2_rv2 = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.rev.2.bt2",
    log:
        outdir + "/log/Align/bowtie2/{genome}/bowtie2_index.log"
    params:
        bowtie2_build_cmd = config.get('tools', {}).get('bowtie2-build') or 'bowtie2-build',
        index_prefix = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome"
    conda:
        config['conda']['run']
    threads:
        8
    shell:
        """
        {params.bowtie2_build_cmd} --threads {threads} {input.fasta}  {params.index_prefix} > {log} 2>&1
        """

rule bowtie2_align:
    input:
        fastqs = get_alignment_input,
        bw2_1   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.1.bt2",
        bw2_2   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.2.bt2",
        bw2_3   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.3.bt2",
        bw2_4   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.4.bt2",
        bw2_rv1 = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.rev.1.bt2",
        bw2_rv2 = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.rev.2.bt2"
    output:
        bam = outdir + "/Align/bam/{genome}/{sample_id}.bam"
    log:
        outdir + "/log/Align/bowtie2/{genome}/{sample_id}/bowtie2_align.log"
    params:
        bowtie2_cmd = config.get('tools',{}).get('bowtie2') or 'bowtie2',
        samtools_cmd = config.get('tools',{}).get('samtools') or 'samtools',
        index_prefix = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome"
    conda:
        config['conda']['run']
    threads:
        10
    shell:
        """
        # 将 input.fastqs 转换为 bash 数组
        fq_array=({input.fastqs})
        fq_count=${{#fq_array[@]}}

        if [ "$fq_count" -eq 2 ]; then
            # 双端逻辑
            {params.bowtie2_cmd} -x {params.index_prefix} \
                -1 "${{fq_array[0]}}" -2 "${{fq_array[1]}}" \
                -N 1 -L 30 \
                --threads {threads} 2> {log} \
                | {params.samtools_cmd} sort -@ {threads} -o {output.bam} - >> {log} 2>&1

        elif [ "$fq_count" -eq 1 ]; then
            # 单端逻辑
            {params.bowtie2_cmd} -x {params.index_prefix} \
                -U "${{fq_array[0]}}" \
                -N 1 -L 30 \
                --threads {threads} 2> {log} \
                | {params.samtools_cmd} sort -@ {threads} -o {output.bam} - >> {log} 2>&1
        else
            echo "Error: Expected 1 or 2 fastq files, got $fq_count" > {log}
            exit 1
        fi
        """