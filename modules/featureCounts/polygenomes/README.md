# featureCounts/polygenomes — 多基因组基因表达定量子模块(featureCounts)

按基因组分层的 featureCounts 定量:每基因组把单端、双端样本分别汇总一次
计数,再合并为一个 count 矩阵。C 形态子模块,共享父目录
`../featureCounts.yaml` 环境。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 布局参数 | 形态 A 拼接,多基因组布局:`{indir}/{genome}/{sample_id}/{sample_id}.bam` |
| `genome_paired_samples` / `genome_single_samples` | genome → 样本列表字典;同基因组样本在一次 featureCounts 调用中合并计数(列 = 样本) |
| `genome.references.<g>.gtf` | 每基因组注释 GTF(`-a`);缺失报错。.smk 读 `genome.references`(复数),json 示例为 `genome.<g>.gtf` 平铺,以 .smk 为准 |
| `ROOT_DIR` | json/schema 声明,顶层读取,当前规则未实际使用 |

### 1.2 输出

```
{outdir}/{genome}/{genome}_single_featureCounts.tsv   # 单端样本合计数表
{outdir}/{genome}/{genome}_paired_featureCounts.tsv   # 双端样本合计数表
{outdir}/{genome}/{genome}_featureCounts.tsv          # 合并 count 矩阵
```

合并逻辑(`featureCounts_merge`,pandas):PE+SE 并存时按 `Geneid` 列 merge,
SE 侧丢弃 `Chr/Start/End/Strand/Length` 注释列、列名正则规整为样本 ID;仅
PE 或仅 SE 时直接复制。日志:`{logdir}/{genome}/featureCounts_*_noMultiple.log`、
`featureCounts_merge.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `output` / `output` / `log` | 输入/输出/日志目录 |
| `Procedure.featureCounts` | `featureCounts -T ... -t exon -g gene_id` | `featureCounts` | featureCounts 可执行路径 |
| `Params.featureCounts.M` | `-M` | `false` | 仅双端规则;计入多重比对读段 |
| `Params.featureCounts.O` | `-O` | `false` | 仅双端规则;允许读段归属多个特征 |
| 线程/内存 | `-T` | `threads=10`,`mem_mb=8000` | `rule_res(...)` 可覆盖 |

固定 flag:计数单位 `-t exon -g gene_id`;双端规则固定 `-B -p --countReadPairs`
(按读段对计数,仅计正确配对)。环境:`conda: "../featureCounts.yaml"`。

直接调用示例(等价引擎命令):

```bash
# 单多样本一次计数
featureCounts -T 10 -t exon -g gene_id -a annotation.gtf \
  -o GRCh38_single_featureCounts.tsv s1.bam s2.bam
# 双端
featureCounts -T 10 -B -p --countReadPairs -t exon -g gene_id -a annotation.gtf \
  -o GRCh38_paired_featureCounts.tsv s3.bam s4.bam
```

## 2. 算法原理与生物学意义

纯工具封装:featureCounts(SubRead 套件)把比对读段按 GTF 注释的
exon→gene_id 层级归属到基因,输出 raw count 矩阵,是 RNA-seq 差异表达
(DESeq2/edgeR)的标准定量步。`--countReadPairs -B -p` 使双端数据按片段
计数且只计 proper pair,避免双端读段重复计数;`-M/-O` 处理多重比对与
多重归属,更敏感但会引入模糊计数,默认关闭。结果怎么读:输出为
featureCounts 原生表,前 6 列注释(Geneid/Chr/Start/End/Strand/Length),
其后每样本一列整数计数;`*_featureCounts.tsv` 合并矩阵可直接进差异分析。

## 3. 参考文献

1. SubRead/featureCounts — https://subread.sourceforge.net/ ,
   https://github.com/ShiLab-Bioinformatics/subread
2. Liao Y, Smyth GK, Shi W. *featureCounts: an efficient general purpose
   program for assigning sequence reads to genomic features*. Bioinformatics
   30(7):923–930, 2014. DOI: 10.1093/bioinformatics/btt656
