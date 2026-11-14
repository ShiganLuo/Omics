# trgt/polygenomes — 多参考基因组 TRGT 串联重复基因分型子模块

TRGT(Tandem Repeat Genotyping Tool)的多参考基因组(polygenomes)变体:
按 `genome` wildcard 对每样本 × 每基因组做长读长串联重复基因分型并出图。
conda 环境共享父目录 `../trgt.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_substring`(config,形态 A) | `{indir}/{sample_id}/{sample_id}[.{bam_substring}].bam` + 对应 `.bam.bai`(长读长比对 BAM 及索引;`bam_substring` 空时省略该段) |
| `genome.<G>.fasta` / `genome.<G>.fai` | 每基因组参考 FASTA 及索引 |
| `genome.<G>.repeat_bed` | 串联重复注释 BED(TRGT 格式:位点 + 重复单元定义) |
| `samples` | 样本 ID 列表(供 subworkflow 注册驱动 DAG) |

### 1.2 输出

```
{outdir}/{genome}/genotype/{sample_id}/{sample_id}.trgt.vcf.gz            # 重复基因分型 VCF
{outdir}/{genome}/genotype/{sample_id}/{sample_id}.trgt.spanning.sorted.bam  # 跨越重复区的 reads
{outdir}/{genome}/plot/{sample_id}/{sample_id}.trgt.repeat.png            # 重复长度可视化
{logdir}/{genome}/{sample_id}/trgt_genotype.log / trgt_plot.log
```

收尾规则 `trgt_result` 聚合 VCF + PNG。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.trgt` | 可执行文件本身 | `trgt` | TRGT 可执行路径 |
| `Params.trgt.karyotype` | `--karyotype` | `XX` | 性染色体核型(XX/XY,影响 X/Y 重复位点基因型判定) |
| `Params.trgt.repeat_id` | `--repeat-id` | `HTT` | trgt plot 出图的重复位点 ID(默认亨廷顿 HTT 位点) |
| `bam_substring` | — | `""` | 输入 BAM 文件名 `.bam` 前的修饰段 |
| `genome.<G>.fasta/.fai` | `--genome` | — | 每基因组参考 |
| `genome.<G>.repeat_bed` | `--repeats` | — | 重复位点注释 BED |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

资源默认:trgt_genotype threads=4、mem_mb=32000;trgt_plot
mem_mb=16000(均经 `rule_res()` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
trgt genotype --num-threads 4 \
  --genome GRCh38.fa --repeats GRCh38_repeats.bed \
  --karyotype XX --reads s1.bam \
  --output-prefix output/GRCh38/genotype/s1/s1.trgt

trgt plot \
  --genome GRCh38.fa --repeats GRCh38_repeats.bed \
  --vcf output/GRCh38/genotype/s1/s1.trgt.vcf.gz \
  --spanning-reads output/GRCh38/genotype/s1/s1.trgt.spanning.sorted.bam \
  --repeat-id HTT \
  --image output/GRCh38/plot/s1/s1.trgt.repeat.png
```

## 2. 算法原理与生物学意义

TRGT 针对长读长数据做串联重复(VNTR/STR)基因分型:对 `repeat_bed` 定义
的每个重复位点,抽取两侧锚定到唯一序列的 spanning reads,逐 read 解析重复
单元数与序列构成(可含基序内变异),聚合出每个等位基因的重复长度/序列,
以 VCF 输出(含每等位基因的重复单元计数);`trgt plot` 按位点渲染重复长
度分布与单 read 结构。生物学意义:动态突变类疾病(亨廷顿病 HTT、脆性 X
FMR1 等)与大量 VNTR 表型由重复扩增致病,长读长 + TRGT 是目前解析重复区
基因型的主要手段;`--karyotype` 保证性染色体位点按 XX/XY 正确判型。结果
阅读:VCF 中 ALT 的重复单元数即等位基因长度;PNG 图看两峰(二倍体两等位
基因)分布与嵌合拖尾。

## 3. 参考文献

1. TRGT —— https://github.com/PacificBiosciences/trgt(源码与文档)
