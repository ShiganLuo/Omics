# bwa-mem2 — 短读长比对模块(index + mem)

封装 bwa-mem2:构建参考基因组索引,把双端/单端 FASTQ 比对到参考基因组,
经 samtools 转 BAM、排序、索引输出坐标排序 BAM。两规则独立,
`bwaMemm2_result` 聚合(规则名拼写即为 `bwaMemm2_result`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | FASTQ 目录;双端 `{indir}/{sample_id}/{sample_id}_1.fq.gz` + `_2.fq.gz`,单端 `{indir}/{sample_id}/{sample_id}single.fq.gz`(注意单端文件名无下划线,见 `get_alignment_input`) |
| `paired_samples` / `single_samples`(config) | 双端/单端样本 ID 列表,决定每个样本取哪种输入;两者皆无则 input 函数报错 |
| `genome.fasta` | 参考基因组 FASTA(`bwaMem2_index` 消费,缺文件建图期即失败) |
| `genome.index_prefix` | 可选:已有 bwa-mem2 索引前缀;存在 `<prefix>.0123` 时直接复用,否则回退模块自建 `{outdir}/index/genome.*` |

### 1.2 输出

```
{outdir}/index/genome.0123            # bwa-mem2 索引(5 文件:.0123 .amb .ann .bwt.2bit.64 .pac)
{outdir}/index/genome.amb
{outdir}/index/genome.ann
{outdir}/index/genome.bwt.2bit.64
{outdir}/index/genome.pac
{outdir}/{sample_id}/{sample_id}.bam  # 坐标排序 BAM(附 .bai 索引,未声明为输出)
```

收尾规则 `bwaMemm2_result` 消费 `{sample_id}.bam`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.bwaMem2` | — | `bwa-mem2` | bwa-mem2 可执行路径(schema/模板键名拼作 `Procedur`,以 .smk 读取的 `Procedure` 为准) |
| `Procedure.samtools` | — | `samtools` | samtools 可执行路径 |
| — | `bwa-mem2 mem -T 0` | 固定 | 最低输出比对分 0(保留全部比对) |
| — | `bwa-mem2 mem -t` | threads | 比对线程(默认 15,`rule_res("bwaMem2_alignment")` 配置) |
| — | `samtools sort -@` | threads | 排序线程 |
| — | `samtools index -@` | threads | 索引线程 |
| threads | — | index/alignment 均 15 | `rule_res` 配置 |
| mem_mb | — | index=64000, alignment=32000 | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
bwa-mem2 index -p index/genome genome.fa

bwa-mem2 mem -T 0 -t 15 index/genome S1_1.fq.gz S1_2.fq.gz 2>> align.log \
  | samtools view -b - \
  | samtools sort -@ 15 > S1.bam
samtools index -@ 15 S1.bam
```

## 2. 算法原理与生物学意义

bwa-mem2 是 BWA-MEM 算法的 SIMD 加速重实现(seed-and-extend:用
FM-index 找 maximal exact matches 作种子,再做带仿射缺口罚分的 Smith-
Waterman 延伸),输出与 bwa-mem 兼容而速度数倍;比对分数低于 `-T`
阈值的比对不输出(本模块 `-T 0` 全保留,过滤交给下游)。samtools
管道把 SAM 流转 BAM、按坐标排序、建索引,得到可随机访问的比对文件。

生物学意义:重测序/靶向测序分析的第一步;比对质量直接决定变异检出的
敏感度与假阳性率。单端/双端由样本列表自动选择;已有索引经
`genome.index_prefix` 复用,避免重复建索引(索引构建内存/时间开销大)。

## 3. 参考文献

1. bwa-mem2 — https://github.com/bwa-mem2/bwa-mem2
2. Li H. *Aligning sequence reads, clone sequences and assembly contigs with
   BWA-MEM*. arXiv:1303.3997, 2013(bwa-mem 算法原文,bwa-mem2 兼容其
   输出格式与参数约定)。
3. samtools — https://www.htslib.org/;Li H et al. *The Sequence
   Alignment/Map format and SAMtools*. Bioinformatics 25(16):2078-2079,
   2009. DOI: 10.1093/bioinformatics/btp352
