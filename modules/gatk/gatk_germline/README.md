# gatk/gatk_germline — 单样本胚系变异检出与过滤子模块

C 形态子模块(HaplotypeCaller → VariantFiltration → result 锚点)。conda
环境共享父目录 `../gatk.yaml`;`known_sites`+`interval` 同时提供时会
`include: "../gatk_bqsr/gatk_bqsr.smk"` 并消费 BQSR 修正后的 BAM。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam` | 已 BQSR 的 BAM(仅当 `reference.known_sites` 与 `reference.interval` 均提供) |
| `{indir}/{sample_id}/{sample_id}.sorted_markdup.bam` | 否则:markdup BAM(回退路径) |
| `genome.fasta` / `genome.fai_index` / `genome.dict_index` | 参考三件套;缺省时回退 `{indir}/index/genome.fa(.fai)` 与 `genome.dict`(上游 prepare 规则产物) |
| `genome.targets` | 可选目标区间文件(`-L`) |
| `reference.interval` / `reference.known_sites` | BQSR 触发键(见上) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.raw.vcf.gz          # HaplotypeCaller 原始变异
{outdir}/{sample_id}/{sample_id}.filtered.vcf.gz     # VariantFiltration 过滤后变异
{outdir}/{sample_id}/{sample_id}.filtered.vcf.gz.csi # bcftools index 索引
{outdir}/{sample_id}/HaplotypeCaller.sh              # 回执(filterHaplotypeCallerVcf.sh 同目录)
logdir/{sample_id}/haplotypeCaller.log               # 运行日志(过滤日志为 haplotypeCaller-filtered.log)
```

`gatk_haplotypeCaller_result` 无 output,只作 DAG 终点锚点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.gatk` | 可执行文件 | null(回退 `gatk`) | GATK 路径 |
| `Procedure.bcftools` | 可执行文件 | null(回退 `bcftools`) | bcftools 路径 |
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | JVM 参数 |
| `Params.gatk.tmp-dir` | `--tmp-dir` | null | 临时目录 |
| `Params.gatk.FS_threshold` | `--filter "FS > X"` | 20.0 | Fisher 链偏阈值,超过打 FS 标 |
| `Params.gatk.QD_threshold` | `--filter "QD < X"` | 2.0 | 按深度归一质量阈值,低于打 QD 标 |
| `Params.gatk.DP_threshold` | `--filter "DP < X"` | 10.0 | 深度阈值,低于打 DP 标 |
| `Params.gatk.QUAL_threshold` | `--filter "QUAL < X"` | 20.0 | 位点质量阈值,低于打 QUAL 标 |
| `genome.fasta` / `fai_index` / `dict_index` / `targets` | `-R` / `-L` | — | 见输入 |
| `reference.interval` / `reference.known_sites` | — | — | 见输入(BQSR 触发) |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录基路径 |

注意:`gatk_germline.json` 把过滤阈值嵌套在
`Params.gatk.VariantFiltration.*`、把 `known_sites`/`interval` 放在
`genome.*` 下,与 .smk 读取的扁平 `Params.gatk.*` 和 `reference.*` 键不
一致,接线以 .smk 为准。过滤的 JEXL 表达式经 `shlex.quote` 引用(`>` `<`
是 shell 重定向字符)。

线程/内存经 `rule_res(...)` 配置:HaplotypeCaller threads 10 / mem 32000,
filterHaplotypeCallerVcf threads 10 / mem 16000。

引擎命令(回执脚本原文):

```bash
gatk --java-options "-Xmx30g" HaplotypeCaller -R genome.fa -I {s}.sorted_markdup.bam -O {s}.raw.vcf.gz [ -L targets.bed ]
gatk --java-options "-Xmx30g" VariantFiltration -R genome.fa -V {s}.raw.vcf.gz -O {s}.filtered.vcf.gz \
    -window 35 -cluster 3 \
    --filter-name FS --filter 'FS > 20.0' --filter-name QD --filter 'QD < 2.0' \
    --filter-name DP --filter 'DP < 10.0' --filter-name QUAL --filter 'QUAL < 20.0'
bcftools index {s}.filtered.vcf.gz
```

## 2. 算法原理与生物学意义

HaplotypeCaller 在候选活跃区做局部重组装(de Bruijn 图),对单倍型重新
打分并给出基因型概率,输出 SNV/InDel(未过滤的 raw VCF)。VariantFiltration
按硬阈值打标(FS/QD/DP/QUAL)并在 35 bp 窗口内把 ≥3 个失败位点聚簇标记
(`-window 35 -cluster 3`)——被打标位点保留在 VCF 中(FILTER 列),由下游
按需排除。生物学意义:硬过滤是单样本胚系变异的经典策略,阈值(如
FS>20 链偏、QD<2 低质量、DP<10 低深度)分别对应测序伪影、比对错误与
覆盖不足;`.filtered.vcf.gz` 可直接用于注释或进入群体联合分型。

## 3. 参考文献

1. GATK — https://gatk.broadinstitute.org ;McKenna A et al. *The Genome
   Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA
   sequencing data*. Genome Res 20(9):1297–1303, 2010.
   DOI: 10.1101/gr.107524.110
2. bcftools — https://www.htslib.org/ ;Danecek P et al. *Twelve years of
   SAMtools and BCFtools*. GigaScience 10(2):giab008, 2021.
   DOI: 10.1093/gigascience/giab008
