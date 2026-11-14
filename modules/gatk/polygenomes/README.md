# gatk/polygenomes — 多基因组 GATK 预处理子模块(gatk_prepare)

GATK 变异检测前置:建参考索引、加读段组(Read Group)、标记重复
(MarkDuplicates),产出坐标排序去重 BAM + QC 度量,供 germline/somatic/
RNA-seq/population 变异检测子模块消费。C 形态子模块,共享父目录
`../gatk.yaml` 环境。

同目录子模块(各有独立 .smk,本 README 不覆盖):
`gatk_germline/`(HaplotypeCaller、filterHaplotypeCallerVcf)、
`gatk_somatic/`(somaticMutect2)、`gatk_RNAseq/`(SplitNCigarReads、
VarientCalling、vcf_filter)、`gatk_population/`(gatk_population_haplotype_caller、
gatk_population_joint_genotyping)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `input_bam_substring` | 形态 A 拼接,多基因组布局:`{indir}/{genome}/{sample_id}/{sample_id}[.{input_bam_substring}].bam`(`input_bam_substring` 默认 `""` 即裸名) |
| `genome.references.<g>.fasta` | 每基因组参考 FASTA;缺失报错。注意 .smk 读 `genome.references`(复数),json 示例平铺 `genome.<g>.fasta`,以 .smk 为准 |

### 1.2 输出

```
{outdir}/index/{genome}/{genome}.fa            # 参考软链
{outdir}/index/{genome}/{genome}.fa.fai        # samtools faidx 索引
{outdir}/index/{genome}/{genome}.dict          # GATK 序列字典
{outdir}/{genome}/{sample_id}/{sample_id}.addReadsGroup.bam(.bai)   # 加 RG(temp)
{outdir}/{genome}/{sample_id}/{sample_id}.sorted_markdup.bam(.bai)  # 去重 BAM
{outdir}/{genome}/{sample_id}/{sample_id}.Markdup-metrics.txt       # 重复度量
```

日志:`{logdir_combine}/index/{genome}/gatk_index.log`、
`{logdir}/{sample_id}/{genome}/addReadsGroup.log`、`MarkDuplicates.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录 |
| `logdir_combine` | — | `log` | 索引日志目录 |
| `input_bam_substring` | 文件名插入段 | `""` | 见输入 |
| `Procedure.gatk` | `gatk CreateSequenceDictionary` / `AddOrReplaceReadGroups` / `MarkDuplicates` | `gatk` | GATK 可执行路径 |
| `Procedure.samtools` | `samtools faidx` / `index` | `samtools` | 索引 |
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | JVM 选项 |
| `Params.gatk.tmp-dir` | `--tmp-dir` | `None` | 非空才传 |
| `Params.gatk.addReadsGroup.RGLB` | `--RGLB` | `lib1` | 文库名 |
| `Params.gatk.addReadsGroup.RGPL` | `--RGPL` | `illumina` | 测序平台 |
| `Params.gatk.addReadsGroup.RGPU` | `--RGPU` | `unit1` | 测序单元 |
| 固定 | `--RGSM <sample_id>`、`-SO coordinate`、`--CREATE_INDEX true`、`--VALIDATION_STRINGENCY SILENT` | — | 样本名取 wildcard |

环境:`conda: "../gatk.yaml"`。

直接调用示例(等价引擎命令):

```bash
gatk CreateSequenceDictionary --java-options -Xmx30g -R genome.fa -O genome.dict
samtools faidx genome.fa -o genome.fa.fai
gatk AddOrReplaceReadGroups --java-options "-Xmx30g" \
  --INPUT s1.bam --OUTPUT s1.addReadsGroup.bam -SO coordinate \
  --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM s1
gatk MarkDuplicates --java-options "-Xmx30g" \
  --INPUT s1.addReadsGroup.bam --OUTPUT s1.sorted_markdup.bam \
  --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
  --METRICS_FILE s1.Markdup-metrics.txt
```

## 2. 算法原理与生物学意义

纯工具封装(GATK Best Practices 预处理三步):① `CreateSequenceDictionary` +
`faidx` 提供 GATK/ samtools 索引;② `AddOrReplaceReadGroups` 写入 RGLB/RGPL/
RGPU/RGSM 读段组——GATK 下游(HaplotypeCaller 等)硬性要求 RG,RGSM 即样本
ID;③ `MarkDuplicates` 按 5' 端坐标识别 PCR/optical 重复并打 FLAG 0x400,
度量文件给重复率(文库复杂度 QC)。去重 BAM 是 BQSR、变异调用、突变检测的
标准起点。结果怎么读:`Markdup-metrics.txt` PERCENT_DUPLICATION 高(如 >30%)
提示文库过度扩增或起始量不足。

## 3. 参考文献

1. GATK — https://gatk.broadinstitute.org/ ,
   https://github.com/broadinstitute/gatk
2. Van der Auwera GA, O'Connor BD. *Genomics in the Cloud*. O'Reilly,
   2020(GATK Best Practices);McKenna A et al. *The Genome Analysis
   Toolkit: a MapReduce framework for analyzing next-generation DNA
   sequencing data*. Genome Research 20(9):1297–1303, 2010.
   DOI: 10.1101/gr.107524.110
