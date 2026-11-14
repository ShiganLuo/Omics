# samtools/sort — BAM 排序与索引子模块

对 BAM 做坐标排序(`samtools sort`)并建索引(`samtools index`),
输出排序后的 BAM + BAI。目录下另有 `polygenomes/` 子模块(多参考基因组
并行场景),不在本文档范围。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config,默认 `input`)+ 样本 ID 布局(形态 A) | `{indir}/{sample_id}/{sample_id}.bam`:未排序(或需重排)的 BAM。路径由模块 input 函数按布局拼接,样本集合由 subworkflow 注册驱动 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.bam       # 坐标排序后的 BAM
{outdir}/{sample_id}/{sample_id}.bam.bai   # 索引
{logdir}/{sample_id}/samtools_sort.log     # 日志
```

收尾规则 `samtools_sort_index_result` 按同样路径聚合,保证两个产物齐备
后 DAG 才算完成。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.samtools` | 可执行文件本身 | `samtools` | 工具可执行路径 |
| `indir` | — | `input` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `log` | 日志目录 |
| `genome.fasta` | — | — | 配置中声明并被读取,但当前规则未消费(预留字段) |

规则默认 `threads = 8`、`mem_mb = 32000`(经 `rule_res("bam_sort")` 可被
流程资源配置覆盖);conda 环境共享父目录 `../samtools.yaml`(容器同源)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
samtools sort -@ 8 -o output/{sample_id}/{sample_id}.bam input/{sample_id}/{sample_id}.bam
samtools index -@ 8 output/{sample_id}/{sample_id}.bam output/{sample_id}/{sample_id}.bam.bai
```

## 2. 算法原理与生物学意义

纯工具封装:`samtools sort` 按参考序列坐标(默认)对 BAM 记录排序,比较器
依次比较参考序列 ID 与 POS,是建索引、IGV 浏览、区域查询与 GATK/变异检测
的前置要求;`samtools index` 生成 BAI 随机访问索引,使按区间抽取 reads 不需
要全文件扫描。生物学意义:坐标有序 + 索引是几乎所有下游分析(覆盖度、变异
calling、可视化)的通用数据形态。结果阅读:`samtools idxstats` / `samtools
flagstat` 可做健全性检查;日志与回执命令记录实际执行参数。

## 3. 参考文献

1. samtools(官方文档)—— https://www.htslib.org/doc/samtools-sort.html
   / https://www.htslib.org/doc/samtools-index.html
2. Li H, Handsaker B, Wysoker A, et al. *The Sequence Alignment/Map format
   and SAMtools*. Bioinformatics 2009, 25(16):2078–2079.
   DOI: 10.1093/bioinformatics/btp352
