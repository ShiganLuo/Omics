# hiphase — 长读长单倍型定相模块(HiPhase)

HiPhase 用长读长的相位信息把杂合变异定相到两条单倍型,同步输出定相后的
BAM(HP/PS 标签)与 phased VCF。目录下另有 `polygenomes/` 子模块(多参考
基因组场景),不在本文档范围。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `bam_dir` / `vcf_dir`(config) | BAM/VCF 根目录(均默认 `indir`);值由 subworkflow 按上游 config `outdir` 拼接(形态 B 的目录键) |
| BAM + 索引 | `{bam_dir}/{sample_id}/{sample_id}[.{input_bam_substring}].bam`;索引为 `.bam.bai`(无 substring 时)或 `{s}.{sub}.bai` |
| VCF + 索引 | `{vcf_dir}/{sample_id}/{sample_id}[.{input_vcf_substring}].vcf.gz` + `.csi` |
| `genome.fasta` | 参考基因组 FASTA |
| `samples` | 样本 ID 列表(供 subworkflow 注册驱动 DAG) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}[.{output_substring}].phased.vcf.gz     # 定相 VCF
{outdir}/{sample_id}/{sample_id}[.{output_substring}].phased.vcf.gz.tbi # 定相 VCF 索引(hiphase 自动写出,规则声明并注册)
{outdir}/{sample_id}/{sample_id}[.{output_substring}].phased.bam        # 定相 BAM(HP/PS 标签)
{logdir}/{sample_id}/hiphase.log
```

收尾规则 `hiphase_result` 聚合全部输出。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.hiphase` | 可执行文件本身 | `hiphase` | HiPhase 可执行路径 |
| `bam_dir` | `--bam` 的根 | `indir` | BAM 目录 |
| `vcf_dir` | `--vcf` 的根 | `indir` | VCF 目录 |
| `input_bam_substring` | — | `""` | 输入 BAM 文件名 `.bam` 前的修饰段 |
| `input_vcf_substring` | — | `""` | 输入 VCF 文件名前的修饰段 |
| `output_substring` | — | `""` | 输出文件名 `.phased.` 前的修饰段 |
| `genome.fasta` | `--reference` | — | 参考 FASTA |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

资源默认:`threads = 8`、`mem_mb = 16000`(经 `rule_res("hiphase_phase")`
可覆盖);conda 环境 `hiphase.yaml`(容器同源)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
hiphase --threads 8 \
  --reference GRCh38.fa \
  --bam s1.bam  --output-bam s1.phased.bam \
  --vcf  s1.vcf.gz --output-vcf s1.phased.vcf.gz
```

## 2. 算法原理与生物学意义

HiPhase 利用长读长(HiFi/ONT)单条 read 跨越多个杂合位点的特性,把杂合
变异聚类分配到两条染色体单倍型,并把定相结果回写:VCF 中以 `|` 分隔基因
型(HP/PS 标签标注单倍型块),BAM 中每条 read 标 `HP`/`PS`,使等位基因
特异的表达、甲基化与结构变异分析成为可能。相比 Whatshap 等,HiPhase 面向
PacBio 数据优化并同时产出定相 BAM 与 VCF。结果阅读:phased VCF 的
`PS` 块大小反映定相连续性;`hiphase` 输出的统计(stdout/日志)给出定相
switch 错误率相关指标;HP-tagged BAM 可直接做单倍型分层分析。

## 3. 参考文献

1. HiPhase —— https://github.com/PacificBiosciences/hiphase(源码与文档)
