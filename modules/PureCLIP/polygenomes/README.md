# PureCLIP/polygenomes — 多参考基因组 CLIP-seq 位点调用子模块

PureCLIP(CLIP-seq 蛋白-RNA 交联位点/区域调用)的多参考基因组
(polygenomes)变体:按 `genome` wildcard 分基因组调用结合位点。conda 环境
共享父目录 `../PureCLIP.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config,默认 `output`)+ 平铺布局 | `{indir}/{sample_id}.dedup.bam` + `.dedup.bam.bai`:去重复 CLIP-seq BAM 及索引(注意平铺文件名,无样本子目录) |
| `genome.<G>.fasta` | 每基因组参考 FASTA;`<G>` 为 genome wildcard |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}.pureclip.sites.bed    # 单核苷酸交联位点 BED
{outdir}/{genome}/{sample_id}.pureclip.region.bed   # 结合区域 BED
{outdir}/{genome}/{sample_id}/pureclip.sh           # 回执(单条引擎命令)
{logdir}/{genome}/{sample_id}/pureclip_run.txt      # 日志
```

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.pureclip` | 可执行文件本身 | `pureclip` | PureCLIP 可执行路径 |
| `Params.pureclip.ld` | `--ld` | true | 为 true 时附加 `--ld` 选项,否则省略 |
| `genome.<G>.fasta` | `-g` | — | 每基因组参考 FASTA |
| `indir` / `outdir` / `logdir` | — | `output` / `output` / `log` | 输入/输出/日志根目录 |

资源默认:`threads = 8`(转 `-nt`)、`mem_mb = 16000`(经
`rule_res("pureclip")` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
pureclip --ld -nt 8 \
  -i s1.dedup.bam -bai s1.dedup.bam.bai \
  -g GRCh38.fa \
  -o output/GRCh38/s1.pureclip.sites.bed \
  -or output/GRCh38/s1.pureclip.region.bed
```

## 2. 算法原理与生物学意义

PureCLIP 用半监督机器学习(HMM + SVM 分类)从 CLIP-seq(或 iCLIP/eCLIP)
比对中分辨蛋白结合峰与背景噪声,输出两类结果:`sites.bed` 给出单核苷酸
分辨率的交联位点(每个峰内最可能的交联位置),`region.bed` 给出结合区域。
`--ld` 启用 linked dinucleotides 模式(跨交联二核苷酸信号),提高交联位点
定位精度(语义详见 PureCLIP 文档)。生物学意义:CLIP-seq 在体内捕获
RNA 结合蛋白(RBP)的直接靶点,单核苷酸位点是 motif、顺式元件与 RBP
调控网络分析的基础。结果阅读:sites.bed 的 score 列为位点置信度;
region.bed 可送 HOMER/motif 分析。polygenomes 场景下每基因组独立调用
(跨物种样本)。

## 3. 参考文献

1. PureCLIP —— https://github.com/skrakau/PureCLIP(源码与文档)
2. Krakau S, Richard H, Marsico A. *PureCLIP: capturing target-specific
   protein-RNA interaction footprints from single-nucleotide CLIP-seq data*.
   Genome Biology 2017, 18:240. DOI: 10.1186/s13059-017-1364-2
