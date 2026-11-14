# population_metadata — 群体样本映射元数据生成模块

纯 Python 生成群体遗传分析所需的样本—群体映射表(单规则
`population_metadata`,不调用外部引擎)。

## 1. 模块输入输出与参数

### 1.1 输入

无文件输入;唯一输入是 config 中的 `populations` 字典(群体名 → 样本 ID
列表),群体名排序后逐个落盘。

### 1.2 输出

```
{outdir}/metadata/population.tsv                 # 样本—群体映射表
{outdir}/metadata/populations/{population}.txt   # 每群体一个样本清单
logdir/population_metadata.log                   # 运行日志
```

`population.tsv` 内容:首行 `sample<TAB>population`,随后每行
`{sample}<TAB>{population}`;`{population}.txt` 每行一个样本 ID。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `populations` | — | {} | 群体名 → 样本 ID 列表;为空时只生成表头 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

线程/内存:`rule_res("population_metadata")` 配置 mem_mb,默认 8000。

## 2. 算法原理与生物学意义

纯元数据生成,无算法可述。用途:群体遗传/群体基因组分析(如
`gatk_population` 的 GenotypeGVCFs 联合分型、PCA、选择清除扫描)需要
明确的样本—群体归属:`population.tsv` 供工具读取群体标签,
`populations/{population}.txt` 供按群体抽样/分层计算。

## 3. 参考文献

无(纯 Python 生成,不依赖外部软件)。
