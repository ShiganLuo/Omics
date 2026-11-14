# featureCounts — 基因水平读段计数模块(Subread/featureCounts)

用 featureCounts 把比对 BAM 按 GTF 注释汇总为基因 × 样本计数矩阵,
双端/单端分表输出,供下游差异表达分析使用。纯工具封装模块。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | BAM 目录(上游 config 的 `outdir`),形态 A 布局:`{indir}/{sample_id}/{sample_id}{.bam_substring}.bam` |
| `bam_substring` (config,顶层) | BAM 文件名中样本 ID 与 `.bam` 之间的后缀段(如 `sorted` → `sample.sorted.bam`);缺省为空 |
| `single_samples` / `paired_samples` (config) | 单端/双端样本列表,分别进入 `featureCounts_single_noMultiple` / `featureCounts_paired_noMultiple`;各自为空时报 `ValueError` |
| `genome.gtf` (config) | 基因注释 GTF(必需,直接下标访问 `config["genome"]["gtf"]`) |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 输出/日志/项目根目录 |

### 1.2 输出

```
{outdir}/all_single_featureCounts.tsv      # 全部单端样本合并计数表
{outdir}/all_paired_featureCounts.tsv      # 全部双端样本合并计数表
{outdir}/featureCounts.done                # 聚合标记(featureCounts_result)
```

回执:`{outdir}/featureCounts_single_noMultiple.sh`、
`featureCounts_paired_noMultiple.sh`;日志
`{logdir}/all/featureCounts/*.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.gtf` | `-a` | — | 注释 GTF,必需 |
| `bam_substring` | — | "" | 见输入;影响 BAM 路径拼接 |
| `Procedure.featureCounts` | 引擎可执行 | `featureCounts` | — |
| (threads) | `-T` | 10 | `rule_res` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

固定命令行(不在 config 暴露):

```bash
# 单端(featureCounts_single_noMultiple)
featureCounts -T 10 -t exon -g gene_id -a annotation.gtf \
    -o all_single_featureCounts.tsv sample1.bam ...

# 双端(featureCounts_paired_noMultiple)
featureCounts -T 10 -B -p --countReadPairs -t exon -g gene_id \
    -a annotation.gtf -o all_paired_featureCounts.tsv sample2.bam ...
```

即:`-t exon -g gene_id` 按 exon 特征计数、gene_id 聚合;双端额外
`-B -p --countReadPairs`(两端均比对、按片段计数)。两规则都**未**加
`-M`,多重比对读段不计入(noMultiple)。

## 2. 算法原理与生物学意义

featureCounts 统计比对到注释特征(本模块:exon)上的读段数,按 `gene_id`
聚合为基因水平计数矩阵,一行一个基因、一列一个样本,末列 `Status` 为
特征属性。要点:

- **noMultiple(不加 `-M`)**:多重比对读段不计数,牺牲一部分灵敏度换取
  计数特异性,避免重复基因家族/假基因读段污染表达量;
- **双端 `-B -p --countReadPairs`**:要求两条 mate 均比对成功,且按
  fragment(read pair)而非单条 read 计数,符合双端文库的分子数语义;
- 输出 TSV 前两列为 Geneid、Chr/Start/End/Strand/Length 等注释列,之后
  每样本一列整数计数,直接作为 DESeq2/edgeR 等差异表达工具的输入。

生物学意义:基因 × 样本计数矩阵是转录组定量的通用货币,驱动差异表达、
富集分析与共表达网络;计数策略(多重比对、片段计数)直接影响低表达基因
与基因家族的定量可靠性。

## 3. 参考文献

1. Subread/featureCounts — https://subread.sourceforge.net/
2. Liao Y, Smyth GK, Shi W. *featureCounts: an efficient general purpose
   program for assigning sequence reads to genomic features*.
   Bioinformatics 30(7):923-930, 2014. DOI: 10.1093/bioinformatics/btt656.
