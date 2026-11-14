# RmrRNA/polygenomes — 多参考基因组 rRNA 序列提取子模块

从基因组注释提取 rRNA 序列生成 `rRNA.fasta`,供下游 rRNA 去除(把 reads
比对到 rRNA 参考过滤核糖体 RNA 污染)使用。conda 环境共享父目录
`../RmrRNA.yaml`(容器同源);引擎脚本共享父目录
`modules/RmrRNA/bin/extract_rRNA.py`(经 `config["ROOT_DIR"]` 拼接)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `genomes.references.<G>.fasta` | 每基因组参考 FASTA;`<G>` 为 genome wildcard。缺失或文件不存在时抛 `ValueError` |
| `genomes.references.<G>.gtf` | 每基因组基因注释 GTF(rRNA 基因特征来源)。缺失或不存在同样抛 `ValueError` |

注意:config 键形态以 .smk 为准(`genomes.references.<G>`);目录内
`RmrRNA.json` 示例写的 `genome.<G>` 段与 .smk 实际读取的键不一致,接线时
按 .smk 形态提供。

### 1.2 输出

```
{outdir}/{genome}/rRNA.fasta                  # 提取的 rRNA 序列
{logdir}/RmrRNA/{genome}/extract_rRNA.log     # 日志
```

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genomes.references.<G>.fasta` | `--fasta` | — | 参考基因组 FASTA |
| `genomes.references.<G>.gtf` | `--gtf` | — | 基因注释 GTF |
| `ROOT_DIR` | — | `.` | 项目根目录,用于定位父目录 bin/ 脚本 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输出/日志根目录 |

资源默认:`threads = 2`(转 `--threads`)、`mem_mb = 8000`(经
`rule_res("extract_rRNA")` 可覆盖)。

直接调用示例(共享引擎脚本 `modules/RmrRNA/bin/extract_rRNA.py`):

```bash
python $ROOT_DIR/modules/RmrRNA/bin/extract_rRNA.py \
  --fasta GRCh38.fa --gtf GRCh38.gtf \
  --output output/GRCh38/rRNA.fasta --threads 2
```

## 2. 算法原理与生物学意义

纯工具封装:脚本按 GTF 注释筛选 rRNA 基因特征(18S/28S/5.8S/5S 等),
从基因组 FASTA 提取对应序列写成 rRNA 参考 FASTA。生物学意义:总 RNA/
转录组数据中核糖体 RNA 占绝大多数 reads,构建 rRNA 参考后可把 reads 先
比对过滤 rRNA,提高 mRNA/smallRNA 分析的有效深度;polygenomes 场景下每个
参考基因组各自生成自己的 rRNA 序列(跨物种样本 rRNA 序列不同)。结果
阅读:输出 FASTA 的序列数与注释 rRNA 基因数一致即可;日志记录提取过程。

## 3. 参考文献

1. 本模块引擎脚本:`modules/RmrRNA/bin/extract_rRNA.py`;
   序列提取基于 samtools faidx 生态 —— https://www.htslib.org/
