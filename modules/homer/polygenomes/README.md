# homer/polygenomes — 多参考基因组 HOMER 峰注释子模块

HOMER `annotatePeaks.pl`:对 narrowPeak 峰做基因组特征注释(启动子、内含
子、基因间区等)与最近基因指派。按 `genome` wildcard 对每基因组注释。conda
环境共享父目录 `../homer.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 样本 ID 布局(形态 A) | `{indir}/{sample_id}/{sample_id}_peaks.narrowPeak`:MACS2 等产出的 narrowPeak 峰文件 |
| `genome.<G>.fasta` | 每基因组参考 FASTA |
| `genome.<G>.gtf` | 每基因组基因注释 GTF(`-gtf` 指定,决定最近基因与特征分类) |
| `samples` | 样本 ID 列表(供 subworkflow 注册驱动 DAG) |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}_peaks.annotatePeaks.txt   # 注释表
{logdir}/{genome}/{sample_id}/homer_annotatepeaks.log
```

收尾规则 `homer_annotatepeaks_result` 聚合注释表。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.annotatePeaks` | 可执行文件本身 | `annotatePeaks.pl` | HOMER 注释工具可执行路径 |
| `genome.<G>.fasta` | 第 2 位置参数(基因组) | — | 参考 FASTA |
| `genome.<G>.gtf` | `-gtf` | — | 基因注释 GTF |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

资源默认:`threads = 1`、`mem_mb = 16000`(经
`rule_res("homer_annotatepeaks")` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
annotatePeaks.pl s1_peaks.narrowPeak GRCh38.fa -gtf GRCh38.gtf \
  > s1_peaks.annotatePeaks.txt
```

## 2. 算法原理与生物学意义

纯工具封装:`annotatePeaks.pl` 对每个峰区间按注释分类基因组位置
(promoter/TSS、exon、intron、intergenic 等)、给出最近基因及距离、TSS
距离分布。生物学意义:ChIP-seq/ATAC-seq 峰的功能解读第一步——判断峰集
整体富集在启动子还是增强子区、把峰关联到靶基因,是 motif 分析与通路富集
的前置。结果阅读:注释表每行一个峰,`Annotation` 列为特征分类,`Closest
Promoter`/`Distance to TSS` 列关联基因;按 Annotation 计数分布可作峰集
QC。polygenomes 场景下每个基因组各自注释(跨物种样本峰来自不同基因组)。

## 3. 参考文献

1. HOMER(官方文档,annotatePeaks.pl)——
   http://homer.ucsd.edu/homer/ngs/annotation.html
2. Heinz S, et al. *Simple combinations of lineage-determining transcription
   factors prime cis-regulatory elements required for macrophage and B cell
   identities*. Molecular Cell 2010, 38(4):576–589.
   DOI: 10.1016/j.molcel.2010.05.004
