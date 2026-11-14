# gatk/polygenomes/gatk_population — 多基因组群体联合分型子模块

C 形态子模块(HaplotypeCaller GVCF → GenomicsDBImport + GenotypeGVCFs):
按 genome × sample 产出 gVCF,再按 genome 做群体联合分型。conda 环境
共享 `../../gatk.yaml`(见 `modules/gatk/`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `sample_bams`(config) | 形态 B:样本 → BAM 完整路径字典,键为纯 `sample`;索引取 `{bam}.bai` |
| `samples` | 样本 ID 列表(与 `sample_bams` 均必填,缺一抛 `ValueError`) |
| `genome.{genome}.fasta` | 每 genome 参考 FASTA |
| `genome.{genome}.fai_index`(或 `.fai`)、`dict_index`(或 `.dict`) | 参考索引(两个键名都兼容) |
| `genome.{genome}.interval` | 可选区间(`-L`),空则全基因组 |

### 1.2 输出

```
{outdir}/{genome}/variants/gvcf/{sample}.g.vcf.gz        # 单样本 gVCF(-ERC GVCF)
{outdir}/{genome}/variants/gvcf/{sample}.g.vcf.gz.tbi    # 索引(IndexFeatureFile)
{outdir}/{genome}/variants/joint/joint.vcf.gz            # 联合分型 VCF
{outdir}/{genome}/variants/joint/joint.vcf.gz.tbi        # 索引
{outdir}/{genome}/variants/joint/sample_map.tsv          # 样本-gVCF 映射(GenomicsDBImport 输入)
{outdir}/{genome}/variants/joint/genomicsdb/             # GenomicsDB 工作区(每次运行 rm -rf 重建)
logdir/{genome}/variants/{sample}/gatk_population_haplotype_caller.log
logdir/{genome}/variants/gatk_population_joint_genotyping.log
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.gatk` | 可执行文件 | `gatk` | GATK 路径 |
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | JVM 参数 |
| `Params.gatk.tmp-dir` | `--tmp-dir` | null | 临时目录 |
| `sample_bams` | — | {} | 样本 → BAM 路径字典(必填) |
| `samples` | — | [] | 样本列表(必填) |
| `genome.{genome}.fasta` / `fai_index` / `dict_index` / `interval` | `-R` / `-L` | — | 见输入 |
| `outdir` / `logdir` | — | `output` / `logs` | 输出按 `{genome}/variants/...` 分层 |

线程/内存经 `rule_res(...)` 配置:两规则 threads 默认 8、mem_mb 默认
32000(`--native-pair-hmm-threads` / `--reader-threads` 跟随 threads)。

引擎命令(回执脚本原文):

```bash
gatk --java-options "-Xmx30g" HaplotypeCaller -R genome.fa -I {bam} -O {sample}.g.vcf.gz -ERC GVCF --native-pair-hmm-threads 8
gatk IndexFeatureFile -I {sample}.g.vcf.gz
rm -rf {outdir}/{g}/variants/joint/genomicsdb
gatk --java-options "-Xmx30g" GenomicsDBImport --sample-name-map sample_map.tsv --genomicsdb-workspace-path .../genomicsdb --reader-threads 8
gatk --java-options "-Xmx30g" GenotypeGVCFs -R genome.fa -V gendb://.../genomicsdb -O joint.vcf.gz
gatk IndexFeatureFile -I joint.vcf.gz
```

## 2. 算法原理与生物学意义

GATK 最佳实践的群体分型两段式:先对每个样本以 `-ERC GVCF` 输出带
<NON_REF> 参考置信块的 gVCF(记录每个位点的基因型似然,而不仅是变异);
再 GenomicsDBImport 汇总全部 gVCF 到 GenomicsDB,GenotypeGVCFs 对全群体
每个位点联合重新定基因型。生物学意义:联合分型让低深度样本借助群体
信息补全基因型,显著提高跨样本基因型一致性,是群体遗传(PCA、选择清除、
关联分析)的标准输入;多 genome 布局下每个参考基因组独立成套输出,互不
混杂。

## 3. 参考文献

1. GATK — https://gatk.broadinstitute.org ;McKenna A et al. *The Genome
   Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA
   sequencing data*. Genome Res 20(9):1297–1303, 2010.
   DOI: 10.1101/gr.107524.110
2. GATK Best Practices(germline short variant discovery)—
   https://gatk.broadinstitute.org/hc/en-us/sections/360007226651
