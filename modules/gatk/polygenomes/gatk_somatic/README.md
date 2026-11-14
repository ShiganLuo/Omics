# gatk/polygenomes/gatk_somatic — 多基因组体细胞突变检测子模块(GATK Mutect2)

polygenomes(每基因组独立目录)布局下的 tumor-vs-normal 体细胞
SNV/INDEL 检测:每基因组、每配对样本跑一次 Mutect2,输出 raw somatic
VCF。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| wildcards `genome` / `normal_sample_id` / `experimental_sample_id` | 基因组 ID 与配对样本(normal vs experimental/tumor) |
| `indir` (config) | 样本 BAM 目录(上游 config 的 `outdir`):`{indir}/{sample}/{sample}.sorted_markdup.bam`;当该基因组 `known_sites` 与 `interval` 同时配置时改用 `{indir}/{sample}/{sample}.sorted_markdup.bqsr.bam`(上游 gatk_bqsr 产出) |
| `genome.<genome_id>.fasta` / `.fai_index` / `.dict_index` (config) | 每基因组参考三件套;缺省时回落到 `{indir}/index/{genome}/genome.fa`(`.fa.fai` / `genome.dict`) |
| `genome.<genome_id>.known_sites` / `.interval` (config) | 同时非空 → 切换 BQSR BAM(见上) |
| `mutect2_parameters` (config,顶层) | 附加 Mutect2 参数字符串(按空白拆分追加) |
| `outdir` / `logdir` (config) | 输出/日志目录 |

环境共享父目录:`conda: "../../gatk.yaml"`、容器 `sif("../../gatk.yaml", ...)`
指向 `modules/gatk/gatk.yaml`。同级非 polygenomes 版见
`modules/gatk/gatk_somatic/README.md`(其 germline_resource / PON 等
增强键**不在**本子模块的 .smk 中,勿混用)。

### 1.2 输出

```
{outdir}/{genome}/mutect2-vcf/{normal}_vs_{experimental}/{normal}_vs_{experimental}.vcf.gz
{outdir}/{genome}/mutect2-vcf/{normal}_vs_{experimental}/mutect2.sh     # 回执
```

日志:`{logdir}/{genome}/all/gatk/{normal}_{experimental}/mutect2.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.<g>.fasta` | `-R` | 回落 `{indir}/index/{g}/genome.fa` | 参考 FASTA |
| `genome.<g>.known_sites` + `.interval` | — | 空 | 两者同时非空 → 输入切 `.bqsr.bam` |
| `mutect2_parameters` | 原样附加(空白拆分) | "" | 透传任意 Mutect2 参数 |
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | JVM 堆 |
| `Params.gatk.tmp-dir` | `--tmp-dir` | null | 非空才附加 |
| `Procedure.gatk` | 引擎可执行 | `gatk` | — |
| (threads) | `--native-pair-hmm-threads` | 10 | `rule_res("somaticMutect2")` 可覆盖 |
| (mem_mb) | — | 32000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
gatk Mutect2 \
    --java-options -Xmx30g \
    -R genome.fa \
    -I normal.sorted_markdup.bam \
    -I tumor.sorted_markdup.bam \
    -normal normal_sample_id \
    -O normal_vs_tumor.vcf.gz \
    --native-pair-hmm-threads 10
```

## 2. 算法原理与生物学意义

Mutect2 在活性区域做局部重组装,对 tumor/normal 两份证据联合建模:
位点在 tumor 中显著富集而 normal 缺席/低频时判为体细胞候选,并输出
TLOD/NLOD、POPAF 等 INFO 证据标签。本子模块只产出 **raw VCF**
(未接 FilterMutectCalls,FILTER 列为 `.`);下游(spectrum 频谱、注释)
需按 INFO 字段自行阈值过滤或另行接入过滤步骤。多基因组布局下同一对
样本对每个基因组各出一份 VCF,可比较不同参考下的体细胞调用差异。
配对语义:normal 由 `-normal` 显式标注,同一 normal 可配多个
experimental。

## 3. 参考文献

1. GATK — https://github.com/broadinstitute/gatk
   (文档 https://gatk.broadinstitute.org/)
2. McKenna A et al. *The Genome Analysis Toolkit: a MapReduce framework
   for analyzing next-generation DNA sequencing data*. Genome Research
   20(9):1297-1303, 2010. DOI: 10.1101/gr.107524.110.
