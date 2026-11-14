# samtools/polygenomes/sort — 多基因组 BAM 排序/索引子模块

polygenomes 布局下的 samtools sort + index:按基因组(`{genome}` 通配符)
整理 BAM 输出布局,产出坐标排序 BAM 与 `.bam.bai` 索引。子模块共享
祖父目录 conda 环境(`../../samtools.yaml`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | BAM 目录;input 函数拼 `{indir}/{sample_id}/{sample_id}.bam` |
| `genome.<g>.fasta` | 每基因组参考 FASTA(在 params 中按 `{genome}` 解析;命令本身不使用) |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.bam       # 坐标排序 BAM
{outdir}/{genome}/{sample_id}/{sample_id}.bam.bai   # BAM 索引
```

收尾规则 `samtools_sort_index_result` 聚合 bam + bai。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.samtools` | — | `samtools` | samtools 可执行路径 |
| — | `samtools sort -@` | threads | 排序线程(默认 8,`rule_res("bam_sort")` 配置) |
| — | `samtools index -@` | threads | 索引线程 |
| mem_mb | — | 32000 | `rule_res("bam_sort")` 配置 |

直接调用示例(脱离流程):

```bash
samtools sort -@ 8 -o out/GRCh38/S1/S1.bam in/S1/S1.bam
samtools index -@ 8 out/GRCh38/S1/S1.bam out/GRCh38/S1/S1.bam.bai
```

## 2. 算法原理与生物学意义

`samtools sort` 做坐标的外部归并排序(内存上限内分块排序再归并),
`index` 建 BAI 随机访问索引;多基因组场景把每基因组的排序结果放到
`{outdir}/{genome}/...` 布局下,供各基因组的下游分析分别消费。纯工具
封装,生物学意义同 samtools 通用排序/索引(下游按区间随机访问的前提)。

## 3. 参考文献

1. samtools — https://www.htslib.org/;Li H et al. *The Sequence
   Alignment/Map format and SAMtools*. Bioinformatics 25(16):2078-2079,
   2009. DOI: 10.1093/bioinformatics/btp352
