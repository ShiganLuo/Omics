# gatk/gatk_bqsr — 碱基质量重校准(BQSR)子模块

GATK BaseRecalibrator + ApplyBQSR 两步:先用已知位点估计碱基质量协变量
重校准表,再把校正后的质量写回 BAM。conda 环境共享父目录
`../gatk.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config,形态 A) | `{indir}/{sample_id}/{sample_id}.sorted_markdup.bam` + `.sorted_markdup.bai`:markdup 后的坐标排序 BAM 及索引 |
| `reference.fasta` | 参考基因组 FASTA(注意 config 键是 `reference`,不是 `genome`) |
| `reference.interval` | 目标区间文件(BED/interval list);提供时以 `-L` 限定分析范围 |
| `reference.known_sites` | 已知位点 VCF 列表(如 dbSNP、Mills indel),逐个转为 `--known-sites`;可为空列表 |

### 1.2 输出

```
{outdir}/bqsr/{sample_id}.recal_data.table            # BaseRecalibrator 输出的重校准表
{outdir}/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam   # ApplyBQSR 输出的校正后 BAM
{outdir}/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam.bai  # 索引(gatk_bqsr_result 聚合声明)
{logdir}/{sample_id}/gatk_bqsr.log / ApplyBQSR.log    # 日志
```

注意:BaseRecalibrator 的输出表路径为 `bqsr/` 子目录,而 ApplyBQSR 的
`table` 输入声明为 `{outdir}/{sample_id}/{sample_id}.recal_data.table`——
两规则对表路径的声明不一致,消费/调试时须留意(以上均为源码原样记录)。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.gatk` | 可执行文件本身 | `gatk` | GATK 可执行路径 |
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | JVM 参数(两规则共用) |
| `Params.gatk.tmp-dir` | `--tmp-dir` | 不传 | GATK 临时目录 |
| `reference.fasta` | `-R` | — | 参考 FASTA |
| `reference.interval` | `-L` | 不传 | 目标区间 |
| `reference.known_sites` | `--known-sites`(可重复) | `[]` | 已知位点 VCF |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

资源默认:BaseRecalibrator threads=8、mem_mb=32000;ApplyBQSR
threads=8(硬编码)、mem_mb=16000。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
gatk BaseRecalibrator --java-options -Xmx30g \
  -R ref.fa -I s1.sorted_markdup.bam \
  --known-sites dbsnp.vcf.gz --known-sites Mills.vcf.gz \
  -L targets.bed -O bqsr/s1.recal_data.table

gatk ApplyBQSR --java-options -Xmx30g \
  -R ref.fa -I s1.sorted_markdup.bam \
  -bqsr s1.recal_data.table \
  -L targets.bed -O s1.sorted_markdup.bqsr.bam
```

## 2. 算法原理与生物学意义

BQSR(Base Quality Score Recalibration)分两步:

1. **BaseRecalibrator**:统计比对中观测到的错配与已知变异位点
   (`--known-sites`,真实变异不计入误差),按协变量(读内机器 cycle、
   碱基上下文 k-mer、原始比对质量等)分桶,经验贝叶斯估计每桶的质量值
   偏差,输出重校准表(`recal_data.table`,含校正前后质量分布);
2. **ApplyBQSR**:按表把校正后的碱基质量写入新 BAM,原 BAM 不变。

生物学意义:测序仪给出的碱基质量存在系统性偏倚(特定 cycle/序列上下文
的错误率被高估或低估),直接变异检测会放大假阳性/假阴性;重校准后质量值
更贴近真实错误率,是 GATK 最佳实践(GATK Best Practices)变异检测的标准
前置步骤。结果阅读:recal 表的 Empirical Quality vs Reported Quality 曲线
偏离对角线的程度反映校正幅度;校正后 BAM 用 `samtools calmd`/`gatk
AnalyzeCovariates` 可复核。

## 3. 参考文献

1. GATK 官方文档(BaseRecalibrator / ApplyBQSR)——
   https://gatk.broadinstitute.org/
2. McKenna A, et al. *The Genome Analysis Toolkit: a MapReduce framework
   for analyzing next-generation DNA sequencing data*. Genome Research 2010,
   20(9):1297–1303. DOI: 10.1101/gr.107524.110
