# trim-galore — 接头/低质量修剪模块(Trim Galore!)

用 Trim Galore!(cutadapt + FastQC 封装)修剪 FASTQ 的接头与低质量末端,
双端、单端各一条规则;`mode: UMI` 时消费 UMI 提取后的 FASTQ。

## 1. 模块输入输出与参数

### 1.1 输入

形态 A:`indir` + 样本列表;`mode` 决定文件后缀:

| 输入 | 规则 | 说明 |
|---|---|---|
| `{indir}/{sample_id}/{sample_id}_1.fq.gz` + `_2.fq.gz` | trimming_Paired | 双端(gzip FASTQ);`mode: UMI` 时为 `_1.umi.fq.gz` / `_2.umi.fq.gz` |
| `{indir}/{sample_id}/{sample_id}.single.fq.gz` | trimming_Single | 单端;`mode: UMI` 时为 `{sample_id}.umi.single.fq.gz` |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_1.fq.gz            # 修剪后 R1(由 {s}_val_1.fq.gz 改名)
{outdir}/{sample_id}/{sample_id}_2.fq.gz            # 修剪后 R2(由 {s}_val_2.fq.gz 改名)
{outdir}/{sample_id}/trimming_statistics_1.txt      # R1 修剪报告(由 {输入名}_trimming_report.txt 改名)
{outdir}/{sample_id}/trimming_statistics_2.txt      # R2 修剪报告
{outdir}/{sample_id}/{sample_id}.single.fq.gz       # 修剪后单端(由 {s}_trimmed.fq.gz 改名)
{outdir}/{sample_id}/trimming_statistics.txt        # 单端修剪报告
{outdir}/{sample_id}/trimming_Paired.sh             # 回执(单端为 trimming_Single.sh)
logdir/{sample_id}/trimming.txt                     # 运行日志
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.trim_galore` | 可执行文件 | null(回退 `trim_galore`) | trim_galore 可执行路径 |
| `Params.trim_galore.quality` | `--quality` | 25 | 末端质量修剪阈值(Phred) |
| `Params.trim_galore.adapters` | `--adapter` | null | 接头序列;null 时交给 Trim Galore 自动检测 |
| `mode` | — | null | 设为 `UMI` 时输入换 `.umi.` 后缀文件 |
| `indir` | — | `output/raw_fastq` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

线程/内存经 `rule_res(...)` 配置:threads 默认 6,mem_mb 双端 8000 / 单端 4000。
单端规则固定加 `--phred33`;双端规则不指定质量编码(Trim Galore 自动判断)。

引擎命令(回执脚本原文,双端):

```bash
trim_galore --paired --cores 6 --quality 25 -o {outdir}/{s} --basename {s} [ --adapter <seq> ] R1.fq.gz R2.fq.gz
mv {outdir}/{s}/{s}_val_1.fq.gz {outdir}/{s}/{s}_1.fq.gz
mv {outdir}/{s}/{s}_val_2.fq.gz {outdir}/{s}/{s}_2.fq.gz
mv {outdir}/{s}/{R1文件名}_trimming_report.txt {outdir}/{s}/trimming_statistics_1.txt
```

## 2. 算法原理与生物学意义

纯工具封装:Trim Galore 先调用 FastQC 推断接头类型(未指定 `--adapter` 时),
再以 cutadapt 做接头切除与 3' 端质量修剪(质量低于 `--quality` 的末端碱基
切除)。生物学意义:接头序列会污染比对与定量结果,低质量末端碱基抬高假
错配率;修剪后的 clean reads 是比对、组装、定量的合格输入。修剪统计
(trimming_statistics*.txt)给出每步切除的 reads 数与碱基数,可用于核对
接头污染程度。

## 3. 参考文献

1. Trim Galore — https://github.com/FelixKrueger/TrimGalore
   (Trim Galore 是 cutadapt + FastQC 的封装,无独立论文;引方法见 cutadapt)
2. Martin M. *Cutadapt removes adapter sequences from high-throughput
   sequencing reads*. EMBnet.journal 17(1):10–12, 2011.
   DOI: 10.14806/ej.17.1.200
