# deepvariant/polygenomes — 多参考基因组 DeepVariant 变异检测子模块

DeepVariant(深度学习变异检测)的多参考基因组(polygenomes)变体:按
`genome` wildcard 对每样本 × 每基因组输出 VCF/gVCF。conda 环境共享父目录
`../deepvariant.yaml`(容器同源);config 可另指定 Docker 容器。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_substring`(config,形态 A) | BAM:`{indir}/{sample_id}/{sample_id}[.{bam_substring}].bam`;索引:`{indir}/{sample_id}/{sample_id}[.{bam_substring}].bai`(注意索引后缀是 `.bai`,不是 `.bam.bai`)。`bam_substring` 空时省略该段 |
| `genome.<G>.fasta` / `genome.<G>.fai` | 每基因组参考 FASTA 及其 .fai 索引;`<G>` 为 genome wildcard |
| `samples` | 样本 ID 列表(供 subworkflow 注册驱动 DAG) |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.vcf.gz        # 变异 VCF
{outdir}/{genome}/{sample_id}/{sample_id}.vcf.gz.csi    # 索引(bcftools index)
{outdir}/{genome}/{sample_id}/{sample_id}.g.vcf.gz      # gVCF
{logdir}/{genome}/{sample_id}/deepvariant.log           # 日志
```

收尾规则 `deepvariant_result` 聚合 VCF + 索引。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.deepvariant` | 可执行文件本身 | `run_deepvariant` | DeepVariant 可执行路径 |
| `Procedure.bcftools` | 可执行文件本身 | `bcftools` | 用于 `bcftools index` |
| `Params.deepvariant.model_type` | `--model_type` | `PACBIO` | 模型类型(PACBIO / WGS / WES / ONT_R104 等,按测序类型选) |
| `bam_substring` | — | `""` | 输入 BAM 文件名中 `.bam` 前的修饰段(如 `dedup`) |
| `container.deepvariant` | — | `docker://google/deepvariant:1.10.0` | Docker 容器 URI;规则中其后另有 `sif("../deepvariant.yaml", ...)` 声明,按声明顺序后者生效 |
| `genome.<G>.fasta` / `.fai` | `--ref` | — | 参考 FASTA 与索引 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

规则默认资源:`threads = 8`(作为 `--num_shards` 并行分片数)、
`mem_mb = 32000`(经 `rule_res("deepvariant_run")` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
run_deepvariant \
  --num_shards 8 --model_type PACBIO \
  --ref GRCh38.fa --reads s1.bam \
  --output_vcf  output/GRCh38/s1/s1.vcf.gz \
  --output_gvcf output/GRCh38/s1/s1.g.vcf.gz
bcftools index output/GRCh38/s1/s1.vcf.gz
```

## 2. 算法原理与生物学意义

DeepVariant 用卷积神经网络把变异检测转化为图像分类任务:对每个候选位点
从比对 reads 渲染 pileup 图像,训练好的网络分类为杂合/纯合变异/纯合参考,
`model_type` 决定测序平台专属的候选生成与模型权重(PACBIO 面向 PacBio
HiFi)。输出 VCF 给出逐位点基因型与质量,gVCF 同时保留参考置信块,便于
多样本联合基因型分析(joint genotyping)。polygenomes 场景下同一样本对多
参考基因组分别调用(如人/鼠跨物种样本),用于 reads 归属核验与污染排查。
结果阅读:VCF 的 `QUAL`/`GQ`/`PL` 评估基因型可信度;比对质量差的区域
gVCF 的 `<NON_REF>` 块提示覆盖不足;`bcftools stats` 做总量 QC。

## 3. 参考文献

1. DeepVariant —— https://github.com/google/deepvariant
2. Poplin R, et al. *A universal SNP and small-indel variant caller using
   deep neural networks*. Nature Biotechnology 2018, 36(3):246–253.
   DOI: 10.1038/nbt.4235
