# fastqc — FASTQ 质控报告模块

用 FastQC 为每个样本的 FASTQ 生成 HTML/ZIP 质控报告(单规则 `fastqc`)。

## 1. 模块输入输出与参数

### 1.1 输入

形态 A:`indir` + 样本列表;`get_fastqc_input` 按 `paired_samples` /
`single_samples` 名单选择文件形态:

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}_1.fq.gz` + `_2.fq.gz` | 双端样本(`sample_id ∈ paired_samples`) |
| `{indir}/{sample_id}/{sample_id}.single.fq.gz` | 单端样本(`sample_id ∈ single_samples`) |

样本不在任一名单时报 `ValueError` 拒绝运行;`sample_id` 受
`wildcard_constraints`(不含 `/`)约束。

### 1.2 输出

```
{outdir}/{sample_id}/                            # 样本输出目录(directory() 声明)
{outdir}/{sample_id}/{sample_id}_fastqc.html     # FastQC 报告(工具产物)
{outdir}/{sample_id}/{sample_id}_fastqc.zip      # FastQC 数据(工具产物)
{outdir}/{sample_id}/{sample_id}.fastqc.{log_suffix}  # 完成 flag(touch)
{outdir}/{sample_id}/fastqc.sh                   # 回执
logdir/{sample_id}/fastqc.{log_suffix}           # 运行日志
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.fastqc` | 可执行文件 | null(回退 `fastqc`) | fastqc 可执行路径 |
| `paired_samples` | — | [] | 双端样本 ID 列表 |
| `single_samples` | — | [] | 单端样本 ID 列表 |
| `log_suffix` | — | `txt` | 日志/flag 文件后缀 |
| `indir` | — | `input` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

线程/内存经 `rule_res("fastqc")` 配置:threads 默认 2,mem_mb 默认 8000。

引擎命令(回执脚本原文):

```bash
fastqc --threads 2 -o {outdir}/{sample_id} -t 2 {输入 FASTQ...}
touch {outdir}/{sample_id}/{sample_id}.fastqc.txt
```

## 2. 算法原理与生物学意义

纯工具封装:FastQC 对 FASTQ 逐 read 统计碱基质量分布、GC 含量、接头污染、
重复序列等指标,输出带 PASS/WARN/FAIL 标注的 HTML 报告。生物学意义:在
进入比对/定量前识别测序质量异常(接头残留、低质量末端、文库污染),是
原始数据交付与流程质控的标准环节。逐项指标解读见 FastQC 官方文档
("Modules" 章节)。

## 3. 参考文献

1. FastQC — https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
2. Andrews S. *FastQC: A Quality Control Tool for High Throughput Sequence
   Data*. (无期刊论文;引用官方工具即可)
