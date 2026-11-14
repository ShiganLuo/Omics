# gatk/gatk_population — 群体联合基因分型子模块(GATK HaplotypeCaller gVCF + GenotypeGVCFs)

population(群体)变异检测:每样本 HaplotypeCaller 出 gVCF,再
GenomicsDBImport + GenotypeGVCFs 联合基因分型出群体 VCF。输入采用
形态 B(`sample_bams` 样本 → 完整路径字典),与布局解耦。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `samples` (config) | 样本 ID 列表;与 `sample_bams`、参考文件同为必需,缺失时 parse 期 `ValueError` |
| `sample_bams` (config) | 形态 B 字典:样本 → BAM 完整路径(字典值须引用上游 config 的 `outdir` 拼接);索引取 `<bam>.bai` |
| `genome.fasta` (config) | 参考基因组 FASTA(必需) |
| `genome.fai_index` / `genome.fai` (config) | FASTA 索引(必需,两个键名任一) |
| `genome.dict_index` / `genome.dict` (config) | 序列字典(必需,两个键名任一) |
| `genome.interval` (config) | 可选目标区间(`-L`) |
| `outdir` / `logdir` (config) | 输出/日志目录 |

环境共享父目录:`conda: "../gatk.yaml"`、容器 `sif("../gatk.yaml", ...)`
指向 `modules/gatk/gatk.yaml`。

### 1.2 输出

```
{outdir}/variants/gvcf/{sample}.g.vcf.gz          # 每样本 gVCF(-ERC GVCF)
{outdir}/variants/gvcf/{sample}.g.vcf.gz.tbi      # IndexFeatureFile 索引
{outdir}/variants/joint/genomicsdb/               # GenomicsDBImport workspace
{outdir}/variants/joint/joint.vcf.gz              # 群体联合分型 VCF
{outdir}/variants/joint/joint.vcf.gz.tbi          # 索引
{outdir}/variants/joint/sample_map.tsv            # 样本名 → gVCF 映射
{outdir}/variants/gvcf/gatk_population_haplotype_caller_{sample}.sh   # 回执
{outdir}/variants/joint/gatk_population_joint_genotyping.sh           # 回执
```

日志:`{logdir}/variants/{sample}/gatk_population_haplotype_caller.log`、
`{logdir}/variants/gatk_population_joint_genotyping.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.fasta` | `-R` | — | 必需 |
| `genome.fai_index`(`fai`) | — | — | 必需 |
| `genome.dict_index`(`dict`) | — | — | 必需 |
| `genome.interval` | `-L` | null | 非空才附加 |
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | JVM 堆 |
| `Params.gatk.tmp-dir` | `--tmp-dir` | null | 非空才附加 |
| `Procedure.gatk` | 引擎可执行 | `gatk` | — |
| (threads) | `--native-pair-hmm-threads` / `--reader-threads` | 8 | `rule_res` 可覆盖 |
| (mem_mb) | — | 32000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
# 每样本 gVCF
gatk --java-options -Xmx30g HaplotypeCaller \
    -R genome.fa -I sample.bam -O sample.g.vcf.gz \
    -ERC GVCF --native-pair-hmm-threads 8
gatk IndexFeatureFile -I sample.g.vcf.gz

# 联合分型
rm -rf outdir/variants/joint/genomicsdb
gatk --java-options -Xmx30g GenomicsDBImport \
    --sample-name-map sample_map.tsv \
    --genomicsdb-workspace-path outdir/variants/joint/genomicsdb \
    --reader-threads 8
gatk --java-options -Xmx30g GenotypeGVCFs \
    -R genome.fa -V gendb://outdir/variants/joint/genomicsdb \
    -O outdir/variants/joint/joint.vcf.gz
gatk IndexFeatureFile -I outdir/variants/joint/joint.vcf.gz
```

## 2. 算法原理与生物学意义

GATK 群体变异检测最佳实践两段式:

- **HaplotypeCaller `-ERC GVCF`**:对每样本在活性区域做局部重组装
  (de Bruijn 图),输出 gVCF——除变异位点外还带 `<NON_REF>` 块记录
  同质参考区间的证据,使低深度样本在联合分型时能利用"无变异"信息;
- **GenomicsDBImport**:把全部 gVCF 导入 GenomicsDB workspace
  (样本 × 位点列式存储),按位点随机访问;
- **GenotypeGVCFs**:逐位点合并全体样本证据,给出群体基因型
  (GT/AD/PL)并重估等位基因频率。

生物学意义:联合分型比逐样本 calling 灵敏度更高(低频/低深度位点
借群体信息补强),且产出统一位点集的群体 VCF,直接支撑群体遗传
分析(admixture、LD、选择扫掠、burden)。`sample_map.tsv` 留档可复现
GenomicsDB 构建。

## 3. 参考文献

1. GATK — https://github.com/broadinstitute/gatk
   (文档 https://gatk.broadinstitute.org/)
2. McKenna A et al. *The Genome Analysis Toolkit: a MapReduce framework
   for analyzing next-generation DNA sequencing data*. Genome Research
   20(9):1297-1303, 2010. DOI: 10.1101/gr.107524.110.
