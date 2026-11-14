# bowtie2 — 短读长比对模块(Bowtie2 建库 + 比对)

Bowtie2 参考基因组建索引与短读长比对封装:双端/单端 FASTQ 比对到参考
基因组,输出 BAM、比对统计与未比对读段。纯工具封装模块。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | FASTQ 目录(上游 config 的 `outdir`):双端 `{indir}/{sample_id}/{sample_id}_1.fq.gz` / `_2.fq.gz`,单端 `{indir}/{sample_id}/{sample_id}.single.fq.gz` |
| `paired_samples` / `single_samples` (config) | 双端/单端样本 ID 列表,决定走 `bowtie2_align_paired` 还是 `bowtie2_align_single` |
| `genome.fasta` (config) | 参考基因组 FASTA;**必填**,parse 期缺失直接 `ValueError` |
| `genome.index_prefix` (config) | 可选:已有 Bowtie2 索引前缀;若 `<prefix>.1.bt2` 存在则直接使用,不再构建 |
| `outdir` / `logdir` / `logdir_combine` (config) | 输出/日志目录 |

注:`get_alignment_input` 中单端候选路径写作 `{indir}/{sample_id}/{sample_id}single.fq.gz`(无点号),与规则 input 的 `{sample_id}.single.fq.gz` 拼写不一致——样本路由以 `paired_samples`/`single_samples` 列表成员身份为准,不以文件存在性为准。已知问题:`bowtie2_align_single` 的 run 块中 `plain_script` 未定义(`shutil.copyfile(script, plain_script)` 处会 NameError),单端路径接线前需修复。

### 1.2 输出

```
# bowtie2_index
{outdir}/index/genome.{1,2,3,4,rev.1,rev.2}.bt2    # 索引 6 件套
{outdir}/index/bowtie2_index.sh                    # 回执

# bowtie2_align_paired / bowtie2_align_single(每样本)
{outdir}/{sample_id}/{sample_id}.bam                        # 比对结果(samtools view -bS)
{outdir}/{sample_id}/{sample_id}_bowtie2_metrics.txt        # --met-file 比对统计
{outdir}/{sample_id}/{sample_id}_unmapped_1.fq.gz           # 未比对 R1(PE)
{outdir}/{sample_id}/{sample_id}_unmapped_2.fq.gz           # 未比对 R2(PE)
{outdir}/{sample_id}/{sample_id}_unmapped.single.fq.gz      # 未比对读段(SE)
{outdir}/{sample_id}/bowtie2_align.sh                       # 回执
```

日志:`{logdir_combine}/bowtie2/bowtie2_index.log`(建库)、
`{logdir}/{sample_id}/bowtie2_align.log`(比对)。

### 1.3 参数

config 位于 `Params.bowtie2`、`Procedure`、`genome` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.fasta` | `bowtie2-build` 输入 | — | 必填;parse 期校验 |
| `genome.index_prefix` | `-x` | `{outdir}/index/genome` | 已有索引前缀;`<prefix>.1.bt2` 存在则复用,否则回落到本模块构建的索引 |
| `Params.bowtie2.sam-append-comment` | `--sam-append-comment` | false(示例 config 为 true) | 在 SAM comment 保留 FASTQ 头部注释 |
| `Procedure.bowtie2-build` | 建库引擎 | `bowtie2-build` | — |
| `Procedure.bowtie2` | 比对引擎 | `bowtie2` | — |
| `Procedure.samtools` | 转 BAM | `samtools` | — |
| (threads) | `--threads` / `-@` | 12(index)/ 8(align) | `rule_res` 可覆盖 |
| (mem_mb) | — | 8000(index)/ 16000(align) | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
# 建库
bowtie2-build --threads 12 genome.fa outdir/index/genome

# 双端比对(sam 转 bam)
bowtie2 -x outdir/index/genome \
    -1 R1.fq.gz -2 R2.fq.gz \
    --threads 8 --sam-append-comment \
    --met-file sample_bowtie2_metrics.txt \
    --un-conc-gz sample_unmapped \
    | samtools view -@ 8 -bS - > sample.bam
```

## 2. 算法原理与生物学意义

Bowtie2 将短读长用 FM-index(BWT)启发式搜索定位到参考基因组,支持
局部比对(soft-clip)与端到端比对,输出 SAM/BAM。本模块把 SAM 管道给
`samtools view -bS` 直接转 BAM,并用 `--un-conc-gz`/`--un-gz` 回收未
比对读段、`--met-file` 输出比对统计(总体比对率、多比对率等)。
`--sam-append-comment` 在 SAM 头部 comment 保留读段原始注释,便于追溯
UMI/样本信息。BAM 与未比对读段供下游变异检测、宿主/共生菌区分
(如 disambiguate 类流程)使用;`*_bowtie2_metrics.txt` 的比对率是文库
质量与参考选择是否合适的快检指标。

## 3. 参考文献

1. Bowtie2 — https://github.com/BenLangmead/bowtie2
   (文档 http://bowtie-bio.sourceforge.net/bowtie2/)
2. Langmead B, Salzberg SL. *Fast gapped-read alignment with Bowtie 2*.
   Nature Methods 9(4):357-359, 2012. DOI: 10.1038/nmeth.1923.
