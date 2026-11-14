# manta — 结构变异(SV)检测模块(germline + somatic)

封装 Illumina Manta 结构变异检测:germline 单样本与 somatic 肿瘤-正常
配对两种模式,各分 configure + run 两步,输出候选/体细胞 SV VCF 及
tabix 索引。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | BAM 目录;input 函数拼 `{indir}/{sample_id}/{sample_id}{bam_suffix}.bam` |
| `bam_suffix`(config) | 可选:BAM 文件名中缀(如 `.dedup`) |
| `genome.fasta` | 参考基因组 FASTA(`configManta.py --referenceFasta`) |
| `samples` | 样本 ID 列表;germline 按样本,somatic 按 `{normal_sample_id}` / `{experimental_sample_id}` 通配符配对 |

### 1.2 输出

```
# germline
{outdir}/{sample_id}/runWorkflow.py(.config.pickle)              # configure 产物
{outdir}/{sample_id}/results/variants/candidateSV.vcf.gz(.tbi)   # 候选 SV

# somatic
{outdir}/{normal_sample_id}_vs_{experimental_sample_id}/runWorkflow.py(.config.pickle)
{outdir}/{normal_sample_id}_vs_{experimental_sample_id}/results/variants/somaticSV.vcf.gz(.tbi)
```

收尾规则 `manta_result` 聚合 germline 的 `candidateSV.vcf.gz(.tbi)`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.manta.threads` | `runWorkflow.py -j` | 8 | Manta 工作流并行 jobs(直读 config,不经 rule_res) |
| `bam_suffix` | — | "" | BAM 文件名中缀 |
| — | `configManta.py --bam`(germline) / `--tumorBam --normalBam`(somatic) | 固定 | configure 阶段参数 |
| `Procedure.manta` / `Procedure.samtools` | — | (模板声明;规则内直接调 `configManta.py`/`python`,未读取) | — |
| mem_mb | — | config 步 8000,run 步 32000 | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
# germline
configManta.py --bam S1.bam --referenceFasta genome.fa --runDir out/S1
python out/S1/runWorkflow.py -j 8
# 产物: out/S1/results/variants/candidateSV.vcf.gz

# somatic
configManta.py --tumorBam T1.bam --normalBam N1.bam \
  --referenceFasta genome.fa --runDir out/N1_vs_T1
python out/N1_vs_T1/runWorkflow.py -j 8
# 产物: out/N1_vs_T1/results/variants/somaticSV.vcf.gz
```

## 2. 算法原理与生物学意义

Manta 从 split/soft-clip reads、异常插入片段大小与配对方向构建 SV
证据图,组装候选断点并用评分模型筛选,一次运行同时输出 DEL/DUP/INV/
BND(易位)与小 InDel;somatic 模式用配对样本模型区分体细胞与胚系
事件,输出 `somaticSV.vcf.gz`。生物学意义:SV 是肿瘤基因组重排、拷贝
数改变与遗传病的重要来源;短读长 SV 检出中断点精度与假阳性控制是
关键,Manta 是 GIAB 短读长 SV 基准的常用调用器之一。结果读法:
`candidateSV.vcf.gz` 过滤态较多,常与 manta 配置的 diploid/somatic
模型评分配合取 PASS;somaticSV 中取 SOMATIC 标签事件。

## 3. 参考文献

1. Manta — https://github.com/Illumina/manta
2. Chen X, Schulz-Trieglaff O, Shaw R, et al. *Manta: rapid detection of
   structural variants and indels for germline and cancer sequencing
   applications*. Bioinformatics 32(8):1220-1222, 2016.
   DOI: 10.1093/bioinformatics/btw033
