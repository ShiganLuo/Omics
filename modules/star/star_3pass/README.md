# star_3pass / star_3pass_gene 设计说明

适用范围：小RNA3'末端转录组测序

`star_3pass` 和 `star_3pass_gene` 的三轮比对均封装在一个 Snakemake 执行规则中。规则内部通过生成带日志和时间戳的 shell script，按顺序执行 pass 1、pass 2、pass 3；subworkflow 只负责 module 配置、`use rule` 和输出编排。

## `star_3pass`

为减少 canonical small RNA reads 错误比对到 small RNA variant genes，采用三轮比对策略。`star_3p_align` 的单个 `run:` 规则执行：

### Pass 1: 全基因组 end-to-end 比对

将 reads 比对到人 hg38 全基因组，使用宽松参数以捕获尽可能多的 reads：

```
STAR --outFilterMultimapNmax 1000
     --alignIntronMin 9999999
     --outFilterMultimapScoreRange 1
     --outFilterMismatchNoverLmax 0.2
     --alignEndsType EndToEnd
```

随后用 bedtools + samtools 提取比对到 small RNA 基因区域的 reads。

### Pass 2: canonical small RNA 参考 end-to-end 比对

将 pass 1 提取的 small RNA reads 重新比对到 canonical small RNA 基因 FASTA 数据库（每个基因含上下游各 50 bp 侧翼序列）。此步骤使用严格参数并 clip read 两端，以去除非模板加尾序列的干扰：

```
STAR --outFilterMultimapNmax 1000
     --outFilterMultimapScoreRange 0
     --outFilterMismatchNoverLmax 0.2
     --outFilterMismatchNoverReadLmax 0.05
     --clip5pNbases 20 0
     --clip3pNbases 0 20
     --alignIntronMin 9999999
     --alignMatesGapMax 500
     --alignEndsType EndToEnd
     --outReadsUnmapped Fastx
```

比对后分为两组：
- **mapped reads**: 成功比对到 canonical small RNA 基因数据库的 reads
- **unmapped reads**: 未能比对到 canonical small RNA 基因数据库的 reads

### Pass 3: 全基因组 local 重比对

将 pass 2 的两组 reads 分别重新比对到人 hg38 全基因组，使用严格 mismatch 和 Local 模式：

```
STAR --outFilterMultimapNmax 1000
     --outFilterMultimapScoreRange 0
     --outFilterMismatchNoverLmax 0.025
     --alignIntronMin 9999999
     --alignMatesGapMax 500
     --alignEndsType Local
```

- **pass 3a**: pass 2 mapped reads 的全基因组 local 比对，再用 bedtools 提取再次落在 small RNA 基因区域的 reads。此步骤用于去除带有测序错误、可能错误比对到 small RNA variant genes 的 canonical small RNA reads。
- **pass 3b**: pass 2 unmapped reads 的全基因组 local 比对，捕获未能比对到 canonical 数据库但在全基因组上有 local 比对的 reads。

### 合并

合并 pass 3a 和 pass 3b 的 BAM，建立 index，得到最终 BAM。

### 三轮比对数据流

```
原始 FASTQ
    │
    ▼ pass1: genome E2E (mismatch 0.2)
    │
    ├── bedtools intersect → small RNA reads
    │
    ▼ pass2: smallRNA ref E2E + clip (mismatch 0.2, readLmax 0.05)
    │
    ├── mapped reads ──→ pass3a: genome Local (mismatch 0.025)
    │                      │
    │                      └── bedtools intersect → canonical small RNA reads
    │
    └── unmapped reads ─→ pass3b: genome Local (mismatch 0.025)
                           │
                           └── non-canonical reads
    │
    ▼
  merge pass3a + pass3b → 最终 BAM
```

## `star_3pass_gene`

`star_3pg_gene_specific` 从已有 `star_3pass` final BAM 开始，在一个 `run:` 规则中执行：

1. 一次性完成全 BAM 到 gene 的 read assignment，并生成 gene-specific FASTQ 和 genomic reference（每个基因含上下游 flank 序列）。
2. 对每个 gene reference 建 STAR index。
3. 对每个 gene 的 per-gene FASTQ 执行与 `star_3pass` 相同结构的三轮 local alignment：
   - pass1: per-gene E2E 比对全部 reads
   - pass2: per-gene E2E + clip 重比对全部 reads，产生 mapped / unmapped 两组
   - pass3a: per-gene Local 比对 pass2 mapped reads
   - pass3b: per-gene Local 比对 pass2 unmapped reads
   - 合并 pass3a + pass3b
4. 取每个 gene 最终 BAM 中 reads > 0 的 pass BAM 执行 Tailer。
5. 拼接所有 gene 的 tail CSV。

默认排除 ambiguous read。

注意：dry-run 或规则注册只能证明 DAG 和命令结构，不能证明实际 STAR 结果或文献等价性。

参考文献：
- https://doi.org/10.1101/2025.01.31.635978
- https://doi.org/10.1073/pnas.2315259121
