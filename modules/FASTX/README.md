# FASTX — fastq_quality_filter 单端 reads 质量过滤模块

把单端 FASTQ 按碱基质量阈值过滤(FASTX-Toolkit 的 `fastq_quality_filter`),
输出过滤后的 FASTQ,供下游比对/定量使用。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config,默认 `output/raw_fastq`)+ 样本 ID 布局(形态 A) | `{indir}/{sample_id}.fq.gz`:单端 gz 压缩 FASTQ,文件名 = 样本 ID + `.fq.gz`。样本集合由 subworkflow 侧注册(node.py outfiles 驱动 DAG),.smk 本身不读 `samples` 列表 |

注意:规则名为 `fastx_quality_filter_single`,仅处理单端数据;源码注释标有
`# not test`,该规则尚未实测。

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.fq.gz                  # 过滤后 FASTQ(temp(),下游消费后由 Snakemake 删除)
{outdir}/{sample_id}/fastx_quality_filter_single.sh     # 回执(单条引擎命令)
{logdir}/{sample_id}/fastx_quality_filter_single.log    # 日志
```

### 1.3 参数

config 位于顶层 + `Params.fastx_toolkit` 段:

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.fastx_toolkit` | 可执行文件本身 | `fastq_quality_filter` | 工具可执行路径 |
| `Params.fastx_toolkit.q` | `-q` | 10 | 保留 read 的最低碱基质量阈值 |
| `Params.fastx_toolkit.p` | `-p` | 100 | read 中质量 ≥ q 的碱基必须达到的百分比 |
| `Params.fastx_toolkit.Q` | `-Q` | 33 | FASTQ 质量值 ASCII 偏移(33 = Phred+33) |
| `indir` | — | `output/raw_fastq` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `log` | 日志目录 |

直接调用示例(模块无 bin/ 引擎脚本;下列为规则拼出的命令形态):

```bash
fastq_quality_filter -Q 33 -q 10 -p 100 \
  -i input/{sample_id}/{sample_id}.fq.gz \
  -o output/{sample_id}/{sample_id}.fq.gz
```

## 2. 算法原理与生物学意义

`fastq_quality_filter` 是 FASTX-Toolkit 的 reads 级质量过滤器,判定规则:
一条 read 仅当其质量值 ≥ `-q` 的碱基占比达到 `-p`% 时才整体保留,否则丢弃
(`-Q` 为质量值的 ASCII 偏移量,决定如何把 FASTQ 字符解码为 Phred 质量)。

生物学意义:测序下机数据含低质量 read,直接进入比对会引入错配与不可靠的
定量;按碱基质量占比过滤能在保留可用数据的同时剔除整体质量差的分子,提升
下游比对、变异与表达定量的可信度。结果阅读:输出即过滤后的 FASTQ,可与
原始 read 数对比回收率;模块为纯工具封装,无自定义判定规则,回执 `.sh`
保存实际执行命令,日志记录起止与报错信息。

## 3. 参考文献

1. FASTX-Toolkit(官方主页,`fastq_quality_filter` 文档)——
   http://hannonlab.cshl.edu/fastx_toolkit/
2. fastx_toolkit 源码 —— https://github.com/agordon/fastx_toolkit
