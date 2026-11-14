# bedtools/polygenomes — iCLIP 交联位点链特异 bigWig 轨道子模块

C 形态子模块(单规则 `iCLIP_bedtools`):把 dedup BAM 转成 iCLIP 交联位点
的链特异 5' 端覆盖 bigWig,按 genome × sample 布局输出。conda 环境共享父
目录 `../bedtools.yaml`(见 `modules/bedtools/`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}.dedup.bam` + `.bam.bai` | 去重 BAM/索引(布局固定后缀 `.dedup`,形态 A 平铺目录) |
| `genome.{genome}.chrom_sizes` | 每个 genome 的染色体长度表(bedtools shift/genomecov 与 bedGraphToBigWig 共用);`{genome}` 来自 wildcard |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.bed         # bamtobed 转换的 BED
{outdir}/{genome}/{sample_id}/{sample_id}.plus.bw     # 正链交联位点 bigWig
{outdir}/{genome}/{sample_id}/{sample_id}.minus.bw    # 负链交联位点 bigWig
{outdir}/{genome}/{sample_id}/iCLIP_bedtools.sh       # 回执(整条流水线脚本)
logdir/{genome}/{sample_id}/bedtools.log              # 运行日志
```

流水线中间文件(shifted.bed、两个 bedgraph、两个 sorted.bedgraph)在脚本
末尾 `rm` 清理。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.bedtools` | 可执行文件 | null(回退 `bedtools`) | bedtools 可执行路径 |
| `Procedure.bedGraphToBigWig` | 可执行文件 | null(回退 `bedGraphToBigWig`) | UCSC bedGraphToBigWig 路径 |
| `genome.{genome}.chrom_sizes` | `-g` / 位置参数 | — | 染色体长度表 |
| `indir` | — | `input` | 输入目录(平铺 `{sample_id}.dedup.bam`) |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

线程/内存经 `rule_res("iCLIP_bedtools")` 配置:threads 默认 4,mem_mb 默认 8000。

引擎命令(回执脚本原文,`{s}` = 样本名):

```bash
bedtools bamtobed -i {s}.dedup.bam > {s}.bed
bedtools shift -m 1 -p -1 -i {s}.bed -g chrom.sizes > {s}.shifted.bed
bedtools genomecov -bg -strand + -5 -scale 1000000 -i {s}.shifted.bed -g chrom.sizes > {s}.plus.bedgraph
bedtools genomecov -bg -strand - -5 -scale 1000000 -i {s}.shifted.bed -g chrom.sizes > {s}.minus.bedgraph
sort -k1,1 -k2,2n {s}.plus.bedgraph > {s}.plus.sorted.bedgraph   # LC_COLLATE=C
bedGraphToBigWig {s}.plus.sorted.bedgraph chrom.sizes {s}.plus.bw
```

## 2. 算法原理与生物学意义

iCLIP(CLIP-seq)中逆转录在 RNA-蛋白交联碱基处截断,cDNA 的 5' 端即
交联位点指示。本模块的判定链:`bamtobed` 转 BED 后 `shift -m 1 -p -1`
把正链特征左移 1 bp、负链右移 1 bp,对齐交联碱基;`genomecov -bg -5`
按链统计 5' 端覆盖并 `-scale 1000000` 归一到 CPM;排序后转 bigWig。
生物学意义:plus/minus 两个 bigWig 给出链特异的蛋白-RNA 交联密度轨道,
在 IGV 中可直接读出目标 RBP 的结合峰与结合基序位置,且多 genome 多样本
输出可横向比较。

## 3. 参考文献

1. bedtools — https://bedtools.readthedocs.io/ ;Quinlan AR, Hall IM.
   *BEDTools: a flexible suite of utilities for comparing genomic features*.
   Bioinformatics 26(6):841–842, 2010. DOI: 10.1093/bioinformatics/btq033
2. bedGraphToBigWig / BigWig — Kent WJ et al. *BigWig and BigBed: enabling
   browsing of large distributed datasets*. Bioinformatics 26(17):2204–2207,
   2010. DOI: 10.1093/bioinformatics/btq351
