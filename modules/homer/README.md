# homer — Peak 功能注释模块(HOMER annotatePeaks)

封装 HOMER `annotatePeaks.pl`:对 narrowPeak 峰做基因组功能注释
(promoter/exon/intron/intergenic 等)并给出最近基因。单规则模块
`homer_annotatepeaks` + 收尾 `homer_annotatepeaks_result`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 峰目录;input 函数拼 `{indir}/{sample_id}/{sample_id}_peaks.narrowPeak`(MACS2/3 narrowPeak 格式) |
| `genome.fasta` | 参考基因组 FASTA(必填,缺失模块加载即报错) |
| `genome.gtf` | 基因注释 GTF/GFF(必填,缺失模块加载即报错;用于 peak 注释) |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_peaks.annotatePeaks.txt  # HOMER 注释表
```

收尾规则 `homer_annotatepeaks_result` 聚合注释表。

### 1.3 参数

无 `Params` 配置段;命令参数固定。

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.annotatePeaks` | — | `annotatePeaks.pl` | HOMER 注释脚本可执行路径 |
| — | `annotatePeaks.pl <peak> <fasta> -gtf <gtf>` | 固定 | peak 文件 + 基因组 + GTF 注释 |
| threads / mem_mb | — | 1 / 16000 | `rule_res("homer_annotatepeaks")` 配置 |

直接调用示例(脱离流程):

```bash
annotatePeaks.pl S1_peaks.narrowPeak genome.fa -gtf gencode.gtf \
  > S1_peaks.annotatePeaks.txt
```

## 2. 算法原理与生物学意义

HOMER `annotatePeaks.pl` 对每个峰区间按 GTF 注释判定功能类别
(promoter/TSS、exon、intron、5'/3' UTR、intergenic 等),给出最近/重叠
基因的 ID、距离 TSS 的距离、链向等(输出为以 `Annotation` 列为核心的
TSV)。生物学意义:ChIP-seq/ATAC-seq 峰的功能分布回答"结合/开放事件
集中在启动子还是远端调控区、关联哪些基因",是 peak 分析的标准注释步。
结果读法:`Annotation` 列看类别分布,`Nearest Refseq`(或同类列)与
`TSS Distance` 列做基因归属;后续可叠加 motif 分析(HOMER
`findMotifsGenome.pl`)。本模块为纯工具封装。

## 3. 参考文献

1. HOMER — http://homer.ucsd.edu/homer/(官方文档:annotatePeaks.pl)。
2. Heinz S, Benner C, Spann N, et al. *Simple combinations of
   lineage-determining transcription factors prime cis-regulatory elements
   required for macrophage and B cell identities*. Molecular Cell
   38(4):576-589, 2010. DOI: 10.1016/j.molcel.2010.05.004
