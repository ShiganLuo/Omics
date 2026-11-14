# UmiTools/extract — UMI/barcode 提取子模块

umi_tools `extract` 的流程封装:从 FASTQ 读段 5' 端按碱基模式字符串切出
UMI(单端)或 UMI + barcode(双端),写入读段名并输出去 UMI 的 FASTQ,
供下游 UMI 感知的比对/去重使用。本目录是 UmiTools 模块的 C 形态子模块。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config)/ 布局参数 | 上游 FASTQ 目录;模块按形态 A 拼路径。单端:`{indir}/{sample_id}.single.fq.gz`;双端:`{indir}/{sample_id}_1.fq.gz` + `{indir}/{sample_id}_2.fq.gz` |

文件形态:gz 压缩 FASTQ;单端/双端各有独立规则,目标由 subworkflow 注册。

### 1.2 输出

```
{outdir}/{sample_id}.umi.single.fq.gz   # 单端去 UMI(temp)
{outdir}/{sample_id}_1.umi.fq.gz        # 双端 R1 去 UMI(temp)
{outdir}/{sample_id}_2.umi.fq.gz        # 双端 R2 去 UMI(temp)
```

输出均为 `temp()`,即下游消费后可被 Snakemake 清理。日志:
`{logdir}/{sample_id}/umi_tools_extract_{single|paired}_run.txt`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` | — | `output` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `log` | 日志目录 |
| `Procedure.umi_tools` | `umi_tools extract` | `umi_tools` | umi_tools 可执行路径 |
| `Procedure.extract_method` | `--extract-method` | `string` | 提取方式;注意 .smk 从 `Procedure` 读取该字段,json 示例把它写在 `Params.umi_tools.extract_method`(该位置不被本 .smk 消费) |
| `Params.umi_tools.bc_pattern` | `--bc-pattern` | `NNNXXXXNN` | R1 碱基模式:N=UMI 碱基,X=barcode 碱基,其余为固定序列 |
| `Params.umi_tools.bc_pattern2` | `--bc-pattern2` | `NNNXXXXNN` | R2 碱基模式(仅双端规则) |
| 线程/内存 | — | `threads=2`,`mem_mb=8000` | `rule_res(...)` 可覆盖 |

环境:共享父目录 yaml —— `conda: "../UmiTools.yaml"`(本子模块目录内无独立
yaml)。

直接调用示例(等价引擎命令):

```bash
# 单端
umi_tools extract --extract-method=string --bc-pattern=NNNXXXXNN \
  -I sample.single.fq.gz -S sample.umi.single.fq.gz
# 双端
umi_tools extract --extract-method=string --bc-pattern=NNNXXXXNN \
  -S sample_1.umi.fq.gz --bc-pattern2=NNNXXXXNN \
  --read2-in=sample_2.fq.gz --read2-out=sample_2.umi.fq.gz
```

## 2. 算法原理与生物学意义

纯工具封装:umi_tools extract 按模式字符串(如 `NNNXXXXNN`:3 bp UMI +
4 bp barcode + 2 bp 校验/固定)从读段 5' 端切分 UMI,把 UMI 追加进读段名
(`@read:UMI`),再输出剥掉 UMI 区段的读段。UMI 是分子标签:同一原始分子的
PCR 拷贝共享同一 UMI,下游按 (位点, UMI) 聚合即可去除 PCR 重复并做纠错
计数,使 RNA/DNA 定量不受扩增偏倚影响。双端模式下 R2 的模式经
`--bc-pattern2` 独立指定。结果怎么读:输出 FASTQ 读段名含 UMI 字段,
统计行数应与输入一致。

## 3. 参考文献

1. UMI-tools — https://github.com/CGATOxford/UMI-tools ,
   https://umi-tools.readthedocs.io/
2. Smith T et al. *UMI-tools: modeling sequencing errors in Unique Molecular
   Identifiers to improve quantification accuracy*. Genome Research
   27(3):491–499, 2017. DOI: 10.1101/gr.209601.116
