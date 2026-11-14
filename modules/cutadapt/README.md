# cutadapt — 接头/低质量修剪模块

用 cutadapt 修剪 FASTQ 的接头、低质量末端与长度异常 read(双端/单端各一
规则);`mode: UMI` 时消费 UMI 提取后的 FASTQ。

## 1. 模块输入输出与参数

### 1.1 输入

形态 A:`indir` + 样本列表;`mode` 决定文件后缀:

| 输入 | 规则 | 说明 |
|---|---|---|
| `{indir}/{sample_id}/{sample_id}_1.fq.gz` + `_2.fq.gz` | trimming_Paired | 双端(gzip FASTQ);`mode: UMI` 时为 `_1.umi.fq.gz` / `_2.umi.fq.gz` |
| `{indir}/{sample_id}/{sample_id}.single.fq.gz` | trimming_Single | 单端;`mode: UMI` 时为 `{sample_id}.umi.single.fq.gz` |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_1.fq.gz                    # 修剪后 R1
{outdir}/{sample_id}/{sample_id}_2.fq.gz                    # 修剪后 R2
{outdir}/{sample_id}/{sample_id}.cutadapt_report.txt        # 双端修剪报告(stdout 重定向)
{outdir}/{sample_id}/{sample_id}.single.fq.gz               # 修剪后单端
{outdir}/{sample_id}/{sample_id}.single.cutadapt_report.txt # 单端修剪报告
{outdir}/{sample_id}/cutadapt_paired.sh                     # 回执(单端为 cutadapt_single.sh)
logdir/{sample_id}/cutadapt_paired.log                      # 运行日志(stderr 重定向)
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.cutadapt` | 可执行文件 | null(回退 `cutadapt`) | cutadapt 可执行路径 |
| `Params.cutadapt.quality` | `-q` | 20 | 3' 端质量修剪阈值(Phred) |
| `Params.cutadapt.adapter_r1` | `-a`(每项一个) | null | R1/单端接头序列列表 |
| `Params.cutadapt.adapter_r2` | `-A`(每项一个) | null | R2 接头序列列表(仅双端规则) |
| `Params.cutadapt.minimum_length` | `-m` | null | 修剪后最短长度,短于则丢弃 |
| `Params.cutadapt.maximum_length` | `-M` | null | 修剪后最长长度 |
| `Params.cutadapt.match_read_wildcards` | `--match-read-wildcards` | false | 接头中的 IUPAC 通配符匹配 read |
| `Params.cutadapt.cut` | `--cut` | null | 每条 read 两端硬切除碱基数 |
| `Params.cutadapt.trimmed_only` | `--trimmed-only` | false | 只保留发生了接头修剪的 read |
| `Params.cutadapt.report` | `--report` | `full` | 报告详细度(full/summary/…) |
| `mode` | — | null | 设为 `UMI` 时输入换 `.umi.` 后缀文件 |
| `indir` | — | `input` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |
| `paired_samples` / `single_samples` | — | [] | 读取但当前 input 函数不消费(见 .smk) |

线程/内存经 `rule_res(...)` 配置:threads 默认 6(`-j`),mem_mb 双端 8000 /
单端 4000。

引擎命令(回执脚本原文,双端;单端去 `-p/-A`,input 单条):

```bash
cutadapt -j 6 -q 20 --report full \
    -o {outdir}/{s}/{s}_1.fq.gz -p {outdir}/{s}/{s}_2.fq.gz \
    [ -a <R1接头> ] [ -A <R2接头> ] [ -m N ] [ -M N ] [ --cut N ] [ --trimmed-only ] \
    R1.fq.gz R2.fq.gz
```

## 2. 算法原理与生物学意义

纯工具封装:cutadapt 用局部比对在 read 中查找接头序列并切除,再按 `-q`
修剪低质量 3' 末端,按 `-m/-M` 过滤长度,`--cut` 做固定位置硬切除
(如 UMI/随机引物残留)。生物学意义:接头与低质量碱基会污染比对、变异
检出与定量;`--trimmed-only` 常用于 3' tag 类文库,只保留带目标接头的
特异 read。报告统计每类接头命中数与修剪碱基数,用于评估文库质量。

## 3. 参考文献

1. cutadapt — https://cutadapt.readthedocs.io/
2. Martin M. *Cutadapt removes adapter sequences from high-throughput
   sequencing reads*. EMBnet.journal 17(1):10–12, 2011.
   DOI: 10.14806/ej.17.1.200
