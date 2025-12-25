rule hisat2_index:
    input:
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output:
        # 1. 移除 lambda。
        # 2. 使用 {genome} 作为通配符。
        # 3. 使用 {{genome}} 让 expand 忽略它，只展开 {idx}。
        # 4. 如果文件名部分也要动态，可以把前缀写在 config 里。
        index = expand(
            outdir + "/genome/{{genome}}/index/hista2/{{genome}}.{idx}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
    threads: 8
    params:
        # 这里的 prefix 必须和上面 output 的路径结构完全一致（去掉 .idx.ht2）
        # 这样 hisat2 才能正确找到输出文件的位置
        prefix = lambda wildcards: outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}",
        HISAT2_BUILD = config["Procedure"]["hisat2-build"]
    log:
        outdir + "/log/genome/{genome}/hisat2_build.log"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        {params.HISAT2_BUILD} -p {threads} {input.fasta} {params.prefix} > {log} 2>&1
        """

rule hisat2_align:
    input:
        fastq = get_alignment_input,
        # 这里在 lambda 内部使用 expand 是允许的，因为 input 支持函数
        index = lambda wildcards: expand(
            outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}.{{idx}}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
    output:
        outfile = outdir + "/Align/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        outdir + "/log/Align/{sample_id}/{genome}/hisat2_align.log"
    threads: 12
    params:
        HISAT2 = config['Procedure']['hisat2'],
        SAMTOOLS = config['Procedure']['samtools'],
        # 关键修正：HISAT2 的 -x 需要的是前缀，而不是文件列表
        # 我们从 input.index 列表中取第一个，然后去掉后缀
        index_prefix = lambda wildcards, input: input.index[0].replace(".1.ht2", ""),
        # 修正 input_params 逻辑：input 现在是一个多索引文件列表 + fastq 列表
        # 建议通过 input.fastq 来指定
        input_params = lambda wildcards, input: \
            f"-1 {input.fastq[0]} -2 {input.fastq[1]}" if len(input.fastq) == 2 else f"-U {input.fastq[0]}"
    shell:
        """
        mkdir -p $(dirname {output.outfile})
        {params.HISAT2} -x {params.index_prefix} \
            {params.input_params} \
            -p {threads} 2> {log} | \
        {params.SAMTOOLS} sort -@ {threads} -@ {threads} -o {output.outfile}
        """
