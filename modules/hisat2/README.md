# hisat2 — RNA-seq 剪接感知比对模块(HISAT2)

HISAT2 建索引与 RNA-seq 比对封装:双端/单端 FASTQ 比对到参考基因组,
输出坐标排序 BAM 与新剪接位点表。纯工具封装模块。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | FASTQ 目录(上游 config 的 `outdir`):双端 `{indir}/{sample_id}/{sample_id}_1.fq.gz` / `_2.fq.gz`,单端 `{indir}/{sample_id}/{sample_id}.single.fq.gz` |
| `paired_samples` / `single_samples` (config) | 双端/单端样本 ID 列表;样本必须在其中之一,否则 `ValueError`;按列表成员身份路由 PE/SE |
| `genome.fasta` (config) | 参考基因组 FASTA(hisat2_index 输入) |
| `genome.index_prefix` (config) | 可选已有索引前缀;`<prefix>.1.ht2` 存在则复用,不构建 |
| `outdir` / `logdir` (config) | 输出/日志目录 |

### 1.2 输出

```
# hisat2_index
{outdir}/index/genome.{1..8}.ht2        # 索引 8 件套
{outdir}/index/hisat2_index.sh          # 回执

# hisat2_align(每样本)
{outdir}/{sample_id}/{sample_id}.bam            # samtools sort 坐标排序 BAM
{outdir}/{sample_id}/{sample_id}.splice.txt     # --novel-splicesite-outfile 新剪接位点
{outdir}/{sample_id}/hisat2_align.sh            # 回执
```

日志:`{logdir}/index/hisat2_build.log`、`{logdir}/{sample_id}/hisat2_align.log`。
`hisat2_result` 为聚合规则(消费 BAM,供 node.py 注册收尾)。

### 1.3 参数

config 位于 `Params.hisat2`、`Procedure`、`genome` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.fasta` | `hisat2-build` 输入 | — | 建索引参考 FASTA |
| `genome.index_prefix` | `-x` | `{outdir}/index/genome` | 已有索引前缀,存在则复用 |
| `Params.hisat2.score_min` | `--score-min` | `L,0,-0.2` | 最小比对分函数(线性:0·len+(-0.2)·len?按 HISAT2 约定 `L,0,-0.2` 表示随长度线性下降的阈值) |
| `Params.hisat2.k` | `-k` | 5 | 每个读段最多报告 5 个比对位置 |
| `Params.hisat2.no-spliced-alignment` | `--no-spliced-alignment` | false | 开启后禁止跨内含子比对(如基因组 DNA 数据) |
| `Params.hisat2.flag_params` | 原样附加 | "" | 透传任意额外 hisat2 参数 |
| `Procedure.hisat2` / `samtools` / `hisat2-build` | 引擎可执行 | `hisat2` / `samtools` / `hisat2-build` | — |
| (threads) | `-p` / `-@` | 8(index)/ 12(align) | `rule_res` 可覆盖 |
| (mem_mb) | — | 8000(index)/ 32000(align) | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
hisat2-build -p 8 genome.fa outdir/index/genome

hisat2 -x outdir/index/genome \
    --score-min L,0,-0.2 \
    -k 5 \
    --novel-splicesite-outfile sample.splice.txt \
    --un-conc-gz outdir/sample/unmapped \
    -1 R1.fq.gz -2 R2.fq.gz \
    -p 12 \
    | samtools sort -@ 12 -o sample.bam
```

## 2. 算法原理与生物学意义

HISAT2 基于 FM-index 的分层图索引,把剪接位点建进图结构,能以极低
错误率做剪接感知比对(RNA-seq 跨内含子读段)。关键参数含义:

- `--score-min L,0,-0.2`:最低可接受比对分随读段长度线性变化,决定
  多低分的次优比对仍被接受;
- `-k 5`:保留至多 5 个比对位置(保留多重比对信息,供下游判断);
- `--novel-splicesite-outfile`:比对中发现的、注释未收录的新剪接位点
  (chr、位点、链),是可变剪接/新转录本发现的线索;
- `--no-spliced-alignment`:关闭跨内含子比对,用于基因组 DNA 数据。

BAM 经 `samtools sort` 坐标排序,供 StringTie 定量、变异检测等下游
使用;`.splice.txt` 与比对率(日志)共同反映 RNA 完整度与注释质量。

## 3. 参考文献

1. HISAT2 — https://github.com/DaehwanKimLab/hisat2
   (文档 https://daehwankimlab.github.io/hisat2/)
2. Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. *Graph-based genome
   alignment and genotyping with HISAT2 and HISAT-genotype*.
   Nature Biotechnology 37:907-915, 2019. DOI: 10.1038/nbt.4220.
