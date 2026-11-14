# deepvariant — 深度学习变异检测模块

封装 Google DeepVariant `run_deepvariant`:对每样本 BAM 用深度神经网络
调用 SNP/小 InDel,输出 VCF/gVCF 并建 CSI 索引。单规则模块
`deepvariant_run` + 收尾 `deepvariant_result`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | BAM 目录;`{indir}/{sample_id}/{sample_id}.bam` + `{sample_id}.bai`(设 `bam_substring` 时为 `{sample_id}.{bam_substring}.bam` + `{sample_id}.{bam_substring}.bai`) |
| `samples` | 样本 ID 列表 |
| `genome.fasta` | 参考基因组 FASTA(必填;缺失时 input 函数报错并给出指引) |
| `genome.fai` | FASTA 索引;缺省回退 `fasta + ".fai"` |
| `bam_substring`(config) | 可选:BAM 文件名中缀,用于区分同一样本多种比对 |
| `tmp_dir`(config) | 可选:临时目录,默认 `{outdir}/tmp`;脚本内 `export TMPDIR`(容器内 GNU parallel 需要可 bind 的 tmp 路径) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.vcf.gz      # 变异 VCF
{outdir}/{sample_id}/{sample_id}.vcf.gz.csi  # CSI 索引
{outdir}/{sample_id}/{sample_id}.g.vcf.gz    # 全基因组 gVCF(可供 GenomicsDB/joint calling)
```

收尾规则 `deepvariant_result` 消费 vcf + 索引。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.deepvariant.model_type` | `--model_type` | `PACBIO` | 测序平台模型(WGS/WES/PACBIO/ONT_R104 等,取 DeepVariant 支持的模型名) |
| `Procedure.deepvariant` | — | `run_deepvariant` | DeepVariant 入口 |
| `Procedure.bcftools` | — | `bcftools` | 索引工具 |
| — | `--num_shards` | threads | 并行分片数(默认 8,`rule_res("deepvariant_run")` 配置) |
| threads / mem_mb | — | 8 / 32000 | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
export TMPDIR=out/tmp
run_deepvariant \
  --num_shards 8 --model_type PACBIO \
  --ref genome.fa --reads S1.bam \
  --output_vcf S1.vcf.gz --output_gvcf S1.g.vcf.gz
bcftools index S1.vcf.gz
```

## 2. 算法原理与生物学意义

DeepVariant 把每个候选位点周围 pileup 渲染成类图像张量,用 Inception
深度网络分类基因型(SNP/InDel),再汇总为 VCF;gVCF 模式额外给出非变异
位点的覆盖度证据,支持队列 joint genotyping。相比传统统计模型,深度
模型对测序错误模式、比对偏差的学习更充分,在 GIAB 基准上 SNP/InDel
精度长期领先。生物学意义:胚系/体细胞变异检出的第一步;`model_type`
须与数据类型匹配(PacBio HiFi 用 PACBIO,短读长 WGS 用 WGS),否则
精度显著下降。结果读法:VCF FILTER=PASS 位点为主要交付,gVCF 供下游
合并基因型。

## 3. 参考文献

1. DeepVariant — https://github.com/google/deepvariant;
   https://github.com/google/deepvariant/blob/main/docs/deepvariant-details.md
2. Poplin R, Chang PC, Alexander D, et al. *A universal SNP and small-indel
   variant caller using deep neural networks*. Nature Biotechnology
   36(3):246-248, 2018. DOI: 10.1038/nbt.4235
