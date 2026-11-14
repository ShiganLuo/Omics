# gatk/gatk_RNAseq — RNA-seq 变异检测模块(SplitNCigarReads + HaplotypeCaller)

GATK RNA-seq 变异检测标准三步:SplitNCigarReads 拆分跨内含子 N-CIGAR
比对 → HaplotypeCaller 调用 → VariantFiltration 硬过滤 + bgzip。子模块
共享父目录 conda 环境(`../gatk.yaml`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 比对目录;input 函数拼 `{indir}/{sample_id}/{sample_id}.sorted_markdup.bam`(标记重复后的排序 BAM) |
| `genome.fasta` | 参考基因组 FASTA(必填,input 函数校验存在) |
| `genome.dict_index` / `genome.fai_index` | GATK 序列字典 / FASTA 索引;缺失时回退 `{indir}/index/genome.dict`、`{indir}/index/genome.fai`,并把 fasta 软链到 `{indir}/index/genome.fa`(让 .dict/.fai 与 fasta 同名,GATK 要求) |
| `samples` | 样本 ID 列表(经 subworkflow 通配符) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.split.bam       # SplitNCigarReads 产物(temp)
{outdir}/{sample_id}/{sample_id}.raw.vcf.gz      # HaplotypeCaller 原始调用
{outdir}/{sample_id}/{sample_id}.filtered.vcf.gz # 过滤后 VCF(bgzip)
```

(无 `_result` 聚合规则;下游直接消费 `filtered.vcf.gz`。)

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | GATK JVM 参数 |
| `Params.gatk.tmp-dir` | `--tmp-dir` | 无 | GATK 临时目录(配置了才传) |
| `Procedure.gatk` | — | `gatk` | GATK 可执行路径 |
| `Procedure.bgzip` | — | `bgzip` | VCF 压缩工具 |
| — | `HaplotypeCaller --dont-use-soft-clipped-bases` | 固定 | RNA-seq 标准建议:忽略软剪切碱基 |
| — | `--standard-min-confidence-threshold-for-calling 20` | 固定 | 调用置信阈值 20 |
| — | `VariantFiltration --window 35 --cluster 3` | 固定 | 35 bp 窗口内 ≥3 个变异聚簇标记 |
| — | `--filter-name FS --filter "FS > 30.0"` | 固定 | 链偏倚过滤 |
| — | `--filter-name QD --filter "QD < 2.0"` | 固定 | 按深度归一化质量过滤 |
| threads / mem_mb | — | 8 / 16000(Split、Calling),8 / 8000(filter) | `rule_res` 配置 |

注:JSON 模板/ schema 声明的 `Params.gatk.VariantFiltration.*_threshold`
(FS/QD/DP/QUAL)与 `genome.known_sites`/`genome.interval` 本 .smk 未读取
(过滤阈值在规则内固定为 FS>30.0、QD<2.0);以 .smk 实际行为为准。

直接调用示例(脱离流程):

```bash
gatk SplitNCigarReads --java-options -Xmx30g \
  -R genome.fa -I S1.sorted_markdup.bam -O S1.split.bam

gatk HaplotypeCaller --java-options -Xmx30g \
  -R genome.fa -I S1.split.bam -O S1.raw.vcf.gz \
  --dont-use-soft-clipped-bases \
  --standard-min-confidence-threshold-for-calling 20

gatk VariantFiltration --java-options -Xmx30g \
  -R genome.fa -V S1.raw.vcf.gz \
  --window 35 --cluster 3 \
  --filter-name FS --filter "FS > 30.0" \
  --filter-name QD --filter "QD < 2.0" \
  -O S1.filtered.vcf
bgzip S1.filtered.vcf
```

## 2. 算法原理与生物学意义

- **SplitNCigarReads**:RNA-seq 比对中跨内含子的 read 带 N(CIGAR 的
  跳跃)操作,直接交给变异调用会被误判;该工具把 read 在 N 处拆成
  外显子内的连续片段,并重修比对质量,是 RNA-seq GATK 最佳实践的
  必经步骤。
- **HaplotypeCaller**:局部重组装单倍型 + PairHMM 计算基因型似然,
  RNA-seq 模式下 `--dont-use-soft-clipped-bases` 避免剪接位点附近的
  软剪切碱基引入假阳性。
- **VariantFiltration**:硬过滤打 FLAG(FS>30.0 链偏倚、QD<2.0 低质量
  密度,35 bp 窗口内 3 个以上变异聚簇);被打标的位点保留但 FILTER
  非 PASS,便于按需取舍。
- **生物学意义**:RNA 层面的变异可用于验证 DNA 突变的转录证据、发现
  RNA 编辑与表达相关多态;剪接感知的预处理是其可信度的前提。结果
  读法:取 `filtered.vcf.gz` 中 FILTER=PASS 位点;FS/QD 标签位点多为
  剪接边缘或低复杂度区,谨慎解读。

## 3. 参考文献

1. GATK — https://gatk.broadinstitute.org/(工具文档:SplitNCigarReads、
   HaplotypeCaller、VariantFiltration)。
2. McKenna A, Hanna M, Banks E, et al. *The Genome Analysis Toolkit: a
   MapReduce framework for analyzing next-generation DNA sequencing data*.
   Genome Research 20(9):1297-1303, 2010. DOI: 10.1101/gr.107524.110
3. Van der Auwera GA, O'Connor BD. *Genomics in the Cloud*. O'Reilly
   Media, 2020(GATK4 最佳实践,RNA-seq 变异检测流程)。
