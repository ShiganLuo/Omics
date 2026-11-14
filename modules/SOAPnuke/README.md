# SOAPnuke — FASTQ 质控过滤模块(双端/单端)

用 SOAPnuke `filter` 对原始 FASTQ 做质量过滤(低质量 read、N 比例、长度),
双端、单端各一条规则,另有两条 result 规则作为 DAG 汇聚锚点。

## 1. 模块输入输出与参数

### 1.1 输入

形态 A:`indir` + 样本列表,input 函数按样本拼路径:

| 输入 | 规则 | 说明 |
|---|---|---|
| `{indir}/{sample_id}/{sample_id}_1.fq.gz` | soapnuke_filter_paired | 双端 R1,gzip 压缩 FASTQ |
| `{indir}/{sample_id}/{sample_id}_2.fq.gz` | soapnuke_filter_paired | 双端 R2 |
| `{indir}/{sample_id}/{sample_id}.single.fq.gz` | soapnuke_filter_single | 单端 reads |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_1.fq.gz          # 过滤后 R1(双端)
{outdir}/{sample_id}/{sample_id}_2.fq.gz          # 过滤后 R2(双端)
{outdir}/{sample_id}/{sample_id}.single.fq.gz     # 过滤后单端 reads
{outdir}/{sample_id}/soapnuke_filter_paired.sh    # 回执(双端;单端同名 _single.sh)
logdir/{sample_id}/soapnuke_filter_paired.log     # 运行日志(单端为 _single.log)
```

SOAPnuke 先以 `-o {outdir}/{sample_id}` 为工作目录、`-C/-D` 指定产物文件名,
再 `mv` 到声明输出路径。`soapnuke_filter_paired_result` /
`soapnuke_filter_single_result` 无 output,只把 clean FASTQ 作为 input,
供 subworkflow 注册终点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.SOAPnuke` | 可执行文件 | null(回退 `SOAPnuke`) | SOAPnuke 可执行路径 |
| `indir` | — | `output/raw_fastq` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

以下过滤参数**硬编码于规则 params**(不在 config 暴露):

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| min_len | `-l` | 15 | 过滤后 read 最小保留长度 |
| low_qual | `-q` | 0.2 | 低质量碱基占比阈值,超过则丢弃整条 read |
| max_n | `-n` | 0.05 | N 碱基占比上限,超过则丢弃整条 read |
| qual_type | `-Q` | 2 | 质量值编码类型(对照见 SOAPnuke 文档) |

线程/内存经 `rule_res(...)` 配置:threads 默认 8,mem_mb 默认 8000。

引擎命令(双端,回执脚本原文;单端去掉 `-2/-D`,input 走 `-1`):

```bash
SOAPnuke filter -l 15 -q 0.2 -n 0.05 \
    -1 {indir}/{s}/{s}_1.fq.gz -2 {indir}/{s}/{s}_2.fq.gz \
    -C {s}_1.fq.gz -D {s}_2.fq.gz -Q 2 \
    -o {outdir}/{s} -T 8
```

## 2. 算法原理与生物学意义

纯工具封装:SOAPnuke 是 BGI 的 FASTQ 质控/预处理工具(支持 MapReduce 加速)。
本模块用 `filter` 子命令按三条规则清洗原始 reads:低质量碱基占比超 `-q`、
N 碱基占比超 `-n`、或修剪后长度不足 `-l` 的 read 被丢弃。生物学意义:
低质量碱基与未知碱基(N)会在比对、变异检出中制造假错配与假变异,按统一
阈值过滤保证多样本数据质量一致,供下游分析使用。

## 3. 参考文献

1. SOAPnuke — https://github.com/BGI-flexbase/SOAPnuke
