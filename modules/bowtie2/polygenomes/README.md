# bowtie2/polygenomes — 多基因组短读长比对子模块(bowtie2)

按基因组分层的 bowtie2 比对模块:为每个 genome 建/取 bowtie2 索引,把每
(基因组, 样本) 的 FASTQ 比对为 BAM,并输出未比对读段与比对度量。多基因组
(polygenomes)形态:输出按 `{genome}/{sample_id}` 分层。C 形态子模块,
共享父目录 `bowtie2.yaml` 环境。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config)/ 布局参数 | 上游 FASTQ 目录,形态 A 拼接:双端 `{indir}/{sample_id}/{sample_id}_1.fq.gz`、`_2.fq.gz`;单端 `{indir}/{sample_id}/{sample_id}.single.fq.gz` |
| `genome.references.<g>.fasta` | 每基因组参考 FASTA(`get_input_for_bowtie2_index` 取用,缺失即报错) |
| `genome.references.<g>.index_prefix` | 可选:预建 bowtie2 索引前缀;存在 `<prefix>.1.bt2` 时直接消费,否则回退模块自建索引 `{outdir}/index/{g}/{g}` |
| `genome_paired_samples` / `genome_single_samples` | genome → 样本列表字典;`get_alignment_input` 据此判单/双端(当前该函数定义于模块但未被规则 input 引用;规则直接按静态路径取文件) |

### 1.2 输出

```
{outdir}/index/{genome}/{genome}.{1.bt2,2.bt2,3.bt2,4.bt2,rev.1.bt2,rev.2.bt2}   # bowtie2 索引
{outdir}/{genome}/{sample_id}/{sample_id}.bam                     # 比对结果(BAM)
{outdir}/{genome}/{sample_id}/{sample_id}_bowtie2_metrics.txt     # bowtie2 --met-file 度量
{outdir}/{genome}/{sample_id}/{sample_id}_unmapped_1.fq.gz        # 双端未比对 R1
{outdir}/{genome}/{sample_id}/{sample_id}_unmapped_2.fq.gz        # 双端未比对 R2
{outdir}/{genome}/{sample_id}/{sample_id}_unmapped.single.fq.gz   # 单端未比对
```

日志:`{logdir}/index/{genome}/bowtie2_index.log`、
`{logdir}/{sample_id}/{genome}/bowtie2_align.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志目录 |
| `logdir_combine` | — | `log/group` | 顶层读取,当前 .smk 未消费 |
| `Procedure.bowtie2` | `bowtie2 -x ...` | `bowtie2` | bowtie2 可执行路径 |
| `Procedure.bowtie2-build` | `bowtie2-build --threads ...` | `bowtie2-build` | 建索引可执行路径 |
| `Params.bowtie2.sam-append-comment` | `--sam-append-comment` | `false` | 为 true 时附加该 flag |
| 线程/内存(index) | `--threads` | `threads=12`,`mem_mb=8000` | `rule_res("bowtie2_index")` 可覆盖 |
| 线程/内存(align) | `--threads` / `samtools view -@` | `threads=8`,`mem_mb=16000` | `rule_res("bowtie2_align_*")` 可覆盖 |

环境:共享父目录 yaml —— `conda: "../bowtie2.yaml"`。

直接调用示例(等价引擎命令):

```bash
bowtie2-build --threads 12 GRCh38.fa output/index/GRCh38/GRCh38
bowtie2 -x output/index/GRCh38/GRCh38 -1 s1_1.fq.gz -2 s1_2.fq.gz \
  --threads 8 --met-file s1_bowtie2_metrics.txt \
  --un-conc-gz s1_unmapped \
  | samtools view -@ 8 -bS - > s1.bam
```

## 2. 算法原理与生物学意义

纯工具封装:bowtie2 是 Burrows-Wheeler 变换(BWT)索引上的短读长比对器,
以 FM-index 搜索 + 动态规划打分(Ma,Bowtie2 论文的并行分块搜索)找到与参考
基因组最优/次优匹配,输出 SAM(此处经 `samtools view -bS` 转 BAM)。多基因组
形态用于宿主-微生物/异种移植等混合样本:每基因组独立比对,后续可用
disambiguate 类工具按唯一比对归属拆分来源。结果怎么读:`_bowtie2_metrics.txt`
给出总体/单端比对率、多重比对率;未比对 FASTQ 可用于下游物种归属或从头分析。

## 3. 参考文献

1. Bowtie 2 — https://github.com/BenLangmead/bowtie2 ,
   https://bowtie-bio.sourceforge.net/bowtie2/
2. Langmead B, Salzberg SL. *Fast gapped-read alignment with Bowtie 2*.
   Nature Methods 9(4):357–359, 2012. DOI: 10.1038/nmeth.1923
