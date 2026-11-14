# gatk/polygenomes/gatk_RNAseq — 多基因组 RNA-seq 变异检测模块

`gatk/gatk_RNAseq` 的 polygenomes 变体:同一套 GATK RNA-seq 变异检测
三步(SplitNCigarReads → HaplotypeCaller → VariantFiltration + bgzip),
按基因组(`{genome}` 通配符)分别执行,基因组名进输入/输出路径段。
conda 环境共享祖父目录 `../../gatk.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 比对目录;input 函数拼 `{indir}/{genome}/{sample_id}/{sample_id}.sorted_markdup.bam`(per-genome 布局) |
| `genome.references.<g>.fasta` | 每个基因组的参考 FASTA(`<g>` = `{genome}` 通配符取值;必填) |
| `genome.references.<g>.dict_index` / `.fai_index` | 序列字典 / FASTA 索引;缺失时回退 `{indir}/index/{genome}/{genome}.dict`、`{genome}.fa.fai`,并把 fasta 软链到 `{indir}/index/{genome}/{genome}.fa` |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.split.bam       # SplitNCigarReads(temp)
{outdir}/{genome}/{sample_id}/{sample_id}.raw.vcf.gz      # HaplotypeCaller 原始调用
{outdir}/{genome}/{sample_id}/{sample_id}.filtered.vcf.gz # 过滤后 VCF(bgzip)
```

(无 `_result` 聚合规则;下游直接消费 `filtered.vcf.gz`。)

### 1.3 参数

与 `gatk/gatk_RNAseq` 完全一致:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | GATK JVM 参数 |
| `Params.gatk.tmp-dir` | `--tmp-dir` | 无 | GATK 临时目录(配置了才传) |
| `Procedure.gatk` / `Procedure.bgzip` | — | `gatk` / `bgzip` | 可执行路径 |
| — | `HaplotypeCaller --dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20` | 固定 | RNA-seq 调用参数 |
| — | `VariantFiltration --window 35 --cluster 3` + `FS > 30.0` / `QD < 2.0` | 固定 | 硬过滤(阈值在规则内固定) |
| threads / mem_mb | — | 8 / 16000(Split、Calling),8 / 8000(filter) | `rule_res` 配置 |

直接调用示例:同 `gatk/gatk_RNAseq` 的 GATK 命令,只是参考与输入输出
换成对应基因组的路径(如 `-R genome1.fa -I in/genome1/S1/S1.sorted_markdup.bam
-O out/genome1/S1/S1.split.bam`)。

## 2. 算法原理与生物学意义

同 `gatk/gatk_RNAseq`:SplitNCigarReads 拆分 N-CIGAR → HaplotypeCaller
局部重组装调用(忽略软剪切碱基)→ VariantFiltration 按 FS/QD 硬过滤。
多基因组场景(如宿主/污染源双基因组、多物种比对)下每个基因组独立
调用,互不混淆参考坐标。详见 `gatk/gatk_RNAseq/README.md`。

## 3. 参考文献

1. GATK — https://gatk.broadinstitute.org/
2. McKenna A et al. *The Genome Analysis Toolkit: a MapReduce framework for
   analyzing next-generation DNA sequencing data*. Genome Research
   20(9):1297-1303, 2010. DOI: 10.1101/gr.107524.110
3. Van der Auwera GA, O'Connor BD. *Genomics in the Cloud*. O'Reilly
   Media, 2020。
