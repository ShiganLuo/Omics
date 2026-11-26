# pbmm2 — PacBio 长读长比对模块

封装 pbmm2(PacBio 版 minimap2):把长读长比对(未比对)BAM 比对到
参考基因组,输出坐标排序 + 索引的 BAM。单规则 `pbmm2_align` + 收尾
`pbmm2_result`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 输入 BAM 目录;input 函数拼 `{indir}/{sample_id}/{sample_id}{bam_suffix}`(`bam_suffix` 默认 `.bam`) |
| `bam_suffix`(config) | 输入文件名后缀,默认 `.bam` |
| `output_bam_suffix`(config) | 输出 BAM 名中缀,默认 `.sorted` |
| `genome.fasta` | 参考基因组 FASTA |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}{output_bam_suffix}.bam   # 排序 BAM(默认 S1.sorted.bam)
{outdir}/{sample_id}/{sample_id}{output_bam_suffix}.bai   # 索引(默认 S1.sorted.bai)
{logdir}/{sample_id}/pbmm2_align.log                     # 日志(pbmm2 INFO 输出 + samtools 消息 + 起止/Executing 标记)
```

(中间 `{sample_id}.raw.bam` 转换成功后删除。)收尾规则 `pbmm2_result`
聚合 bam + bai。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.pbmm2` | — | `pbmm2` | pbmm2 可执行路径 |
| `Procedure.samtools` | — | `samtools` | samtools 可执行路径 |
| — | `pbmm2 align --sort --num-threads` | threads | 比对并排序(线程默认 16,`rule_res("pbmm2_align")` 配置) |
| — | `pbmm2 align --log-level` | `INFO` | 固定 `INFO`;pbmm2 默认 `WARN` 级别正常运行零输出,不开则日志只有 samtools 消息 |
| — | `samtools sort -@` / `samtools index -@` | threads | 二次排序 + 索引 |
| mem_mb | — | 32000 | `rule_res("pbmm2_align")` 配置 |

直接调用示例(脱离流程):

```bash
pbmm2 align --sort --num-threads 16 --log-level INFO genome.fa S1.bam S1.raw.bam
samtools sort -@ 16 S1.raw.bam -o S1.sorted.bam
samtools index -@ 16 S1.sorted.bam -o S1.sorted.bai
rm S1.raw.bam
```

## 2. 算法原理与生物学意义

pbmm2 是 minimap2 的 PacBio 封装(默认参数面向 HiFi/CLR 长读长,
`--sort` 直接输出坐标排序 BAM,并处理 PacBio BAM 的 RG/电影名元数据);
minimap2 用 minimizer 种子 + 链接 + 仿射缺口 Smith-Waterman 做长读长
快速比对。模块在其上再做 samtools sort/index,保证下游(如 DeepVariant/
pbsv/TRGT)拿到标准排序索引 BAM。生物学意义:HiFi 长读长比对是结构
变异、重复区、甲基化分析的基础;pbmm2 对 PacBio 元数据的兼容省去手工
修 RG 的麻烦。结果读法:比对统计可用 samtools flagstat/idxstats 检查。

## 3. 参考文献

1. pbmm2 — https://github.com/PacificBiosciences/pbmm2
2. Li H. *Minimap2: pairwise alignment for nucleotide sequences*.
   Bioinformatics 34(18):3094-3100, 2018. DOI: 10.1093/bioinformatics/bty191
3. samtools — https://www.htslib.org/
