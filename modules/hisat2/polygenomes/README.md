# hisat2/polygenomes — 多基因组 RNA/短读长比对子模块(hisat2)

按基因组分层的 hisat2 比对模块:为每个 genome 建/取 HISAT2 索引,把每
(基因组, 样本) 的 FASTQ(单/双端自动判定)比对、排序为 BAM,并输出新剪接
位点。多基因组(polygenomes)形态:输出按 `{genome}/{sample_id}` 分层。
C 形态子模块,共享父目录 `hisat2.yaml` 环境。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config)/ 布局参数 | 上游 FASTQ 目录;`get_alignment_input` 按样本归属拼路径:双端 `{indir}/{sample_id}/{sample_id}_1.fq.gz`、`_2.fq.gz`,单端 `{indir}/{sample_id}/{sample_id}.single.fq.gz` |
| `genome_paired_samples` / `genome_single_samples` | genome → 样本列表字典,决定每样本单/双端;样本不在任一列表即报错 |
| `genome.references.<g>.fasta` | 建索引用参考 FASTA(`get_input_for_hisat2_index` 取用,缺失报错)。注意 .smk 读 `genome.references`(复数),json 示例为 `genome.<g>.fasta`,以 .smk 为准 |
| `genome.<g>.hisat2_index_prefix` | 可选:预建索引前缀;存在 `<prefix>.1.ht2` 即直接消费,否则用模块自建 `{outdir}/index/{g}/{g}`(json 示例键 `index_prefix` 不被本 .smk 消费) |

### 1.2 输出

```
{outdir}/index/{genome}/{genome}.{1..8}.ht2             # HISAT2 索引(8 个文件)
{outdir}/{genome}/{sample_id}/{sample_id}.bam           # 比对 + samtools sort 排序
{outdir}/{genome}/{sample_id}/{sample_id}.splice.txt    # 新剪接位点(--novel-splicesite-outfile)
```

未比对读段前缀:`{outdir}/{genome}/{sample_id}/unmapped`(`--un-conc-gz`)。
日志:`{logdir_combine}/index/{genome}/hisat2_build.log`
(`logdir_combine` 默认 `log/combine`)、`{logdir}/{sample_id}/{genome}/hisat2_align.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志目录 |
| `logdir_combine` | — | `log/combine` | 建索引日志目录 |
| `Procedure.hisat2` | `hisat2 -x ...` | `hisat2` | hisat2 可执行路径 |
| `Procedure.hisat2-build` | `hisat2-build -p ...` | `hisat2-build` | 建索引可执行路径 |
| `Procedure.samtools` | `samtools sort` | `samtools` | 排序可执行路径 |
| `Params.hisat2.score_min` | `--score-min` | `L,0,-0.2` | 最小比对得分函数 |
| `Params.hisat2.k` | `-k` | `5` | 报告至多 k 个比对位置 |
| `Params.hisat2.no-spliced-alignment` | `--no-spliced-alignment` | `false` | true 时附加(禁用剪接比对) |
| `Params.hisat2.flag_params` | 原样透传 | `""` | 追加任意 bowtie 风格 flag |
| 线程/内存 | `-p` / `-@` | index: `threads=8, mem_mb=8000`;align: `threads=12, mem_mb=32000` | `rule_res(...)` 可覆盖 |

注:`Params.hisat2.*` 四个字段仅由 .smk 消费(json/schema 未声明),缺省即用默认值。
环境:共享父目录 yaml —— `conda: "../hisat2.yaml"`(align 规则中 conda 指令
出现两次,`"../hisat2.yaml"` 与 `"hisat2.yaml"`,以实际 Snakemake 解析为准)。

直接调用示例(等价引擎命令):

```bash
hisat2-build -p 8 genome.fa output/index/GRCh38/GRCh38
hisat2 -x output/index/GRCh38/GRCh38 --score-min L,0,-0.2 -k 5 \
  --novel-splicesite-outfile s1.splice.txt --un-conc-gz s1_unmapped \
  -1 s1_1.fq.gz -2 s1_2.fq.gz | samtools sort -@ 12 -o s1.bam
```

## 2. 算法原理与生物学意义

纯工具封装:HISAT2 基于分层全基因组 FM 索引(Hierarchical Graph FM
index)的剪接感知比对器,能把跨内含子的 RNA-seq 读段拆段比对,并从数据中
学习新剪接位点。`--novel-splicesite-outfile` 记录样本新发现的供体/受体位点,
供 StringTie 等组装工具参考;`-k 5` 保留多重比对信息,利于转录本定量。
结果怎么读:BAM 为坐标排序;`splice.txt` 为剪接位点表(染色体、供体、受体、
链),可与注释比较找新转录事件。

## 3. 参考文献

1. HISAT2 — https://github.com/DaehwanKimLab/hisat2 ,
   https://daehwankimlab.github.io/hisat2/
2. Kim D, Paggi JM, Park C, et al. *Graph-based genome alignment and
   genotyping with HISAT2 and HISAT-genotype*. Nature Biotechnology
   37(8):907–915, 2019. DOI: 10.1038/s41587-019-0201-4
