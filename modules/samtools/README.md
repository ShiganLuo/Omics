# samtools — BAM 比对统计(flagstat)模块

对 `indir` 下每样本 BAM 运行 `samtools flagstat`,产出比对质控统计表。
同目录下 `sort/`(比对排序)、`polygenomes/`(多基因组版)为独立子模块,各有
自己的 `.smk` 与环境 yaml,本 README 只覆盖 `samtools.smk`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config)/ 布局参数 | 上游输出目录;模块按形态 A 自行拼 `{indir}/{sample_id}/{sample_id}.bam`。样本集合不由模块声明,由 subworkflow/node.py 注册的目标文件驱动 |

文件形态:BAM(`samtools flagstat` 读 BAM 头与比对记录,无需索引)。

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_flagstat.txt   # flagstat 统计表
{outdir}/{sample_id}/bam_flagstat.sh            # 回执(单条引擎命令)
```

日志:`{logdir}/{sample_id}/flagstat.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` | — | `input` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `log` | 日志目录 |
| `Procedure.samtools` | `samtools flagstat` | `samtools` | samtools 可执行路径 |
| `Params.samtools.onlykeep_properpair` | — | `false` | json/schema 中声明,当前 `samtools.smk` 未消费(仅保留 proper pair 的开关预留位) |
| 资源 | — | `mem_mb=8000` | `rule_res("bam_flagstat")` 可覆盖 |

环境:`conda: "samtools.yaml"` / `container: sif("samtools.yaml", ...)`(模块目录内 yaml)。

直接调用示例(等价引擎命令,脱离流程):

```bash
samtools flagstat input/sample/sample.bam > output/sample/sample_flagstat.txt
```

## 2. 算法原理与生物学意义

纯工具封装:`samtools flagstat` 逐条扫描 BAM 比对记录,按 flag 分类计数并输出
「QC-passed + QC-failed」双计数行:总读段数、secondary/supplementary、
duplicates、mapped(映射率)、paired in sequencing、read1/read2、properly
paired(双端正确配对率)、单端/双端 mate 未映射等。结果怎么读:重点看
`mapped`(比对率)与 `properly paired`(文库质量)两行的百分比;flagstat
是比对后第一道质控快检,异常低的映射率提示参考基因组/物种不匹配或文库问题。

## 3. 参考文献

1. samtools — https://www.htslib.org/ ,https://github.com/samtools/samtools
2. Li H et al. *The Sequence Alignment/Map format and SAMtools*.
   Bioinformatics 25(16):2078–2079, 2009. DOI: 10.1093/bioinformatics/btp352
