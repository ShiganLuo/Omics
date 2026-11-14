# hiphase/polygenomes — 多基因组单样本单倍型定相子模块

C 形态子模块(单规则 `hiphase_phase` + 锚点 `hiphase_result`):对每个
genome × sample 用 HiPhase 把杂合变异定相到两条单倍型,输出 phased VCF
与 phased BAM。conda 环境共享父目录 `../hiphase.yaml`(见
`modules/hiphase/`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{bam_dir}/{sample_id}/{sample_id}{<.input_bam_substring>}.bam` | 长读长比对 BAM;`bam_dir` 独立于 `indir`(形态 B 目录键) |
| 对应 `.bai` | 索引:`input_bam_substring` 非空时为 `{s}.{sub}.bai`,空时为 `{s}.bam.bai` |
| `{vcf_dir}/{sample_id}/{sample_id}{<.input_vcf_substring>}.vcf.gz` + `.csi` | 样本变异集(待定相) |
| `genome.{genome}.fasta` | 每个 genome 的参考基因组 FASTA(`{genome}` 来自 wildcard) |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.phased.vcf.gz   # 定相 VCF
{outdir}/{genome}/{sample_id}/{sample_id}.phased.bam      # 按单倍型标注的 BAM
{outdir}/{genome}/{sample_id}/hiphase.sh                  # 回执
logdir/{genome}/{sample_id}/hiphase.log                   # 运行日志
```

`output_substring` 非空时输出名变为 `{sample_id}.{output_substring}.phased.*`。
`hiphase_result` 无 output,只作 DAG 终点锚点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.hiphase` | 可执行文件 | `hiphase` | hiphase 可执行路径 |
| `bam_dir` | — | `input`(json 同值) | BAM 所在根目录(每样本一子目录) |
| `vcf_dir` | — | `input` | VCF 所在根目录(每样本一子目录) |
| `input_bam_substring` | — | "" | 输入 BAM 名后缀(如 `dedup` → `{s}.dedup.bam`) |
| `input_vcf_substring` | — | "" | 输入 VCF 名后缀(如 `sv` → `{s}.sv.vcf.gz`) |
| `output_substring` | — | "" | 输出文件名中 `phased` 前的附加后缀 |
| `genome.{genome}.fasta` | `--reference` | — | 每 genome 参考 FASTA |
| `samples` | — | [] | 样本 ID 列表 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录基路径(输出按 `{genome}/{sample_id}` 分层) |

线程/内存经 `rule_res("hiphase_phase")` 配置:threads 默认 8,mem_mb 默认 16000。

引擎命令(回执脚本原文):

```bash
hiphase --threads 8 --reference {genome.fa} \
    --bam {in.bam} --output-bam {outdir}/{g}/{s}/{s}.phased.bam \
    --vcf {in.vcf.gz} --output-vcf {outdir}/{g}/{s}/{s}.phased.vcf.gz
```

## 2. 算法原理与生物学意义

纯工具封装:HiPhase 综合长读长的相位信息(同一 read/HiFi 分子携带的相邻
杂合位点)与 VCF 中的杂合变异,把等位基因定相成两条单倍型,并把每条
read 归属到 hap1/hap2(HP 标签),输出定相 VCF 与 haplotype-resolved BAM。
生物学意义:单倍型分辨把复合杂合变异拆到两条染色体上,是结构变异定相、
等位基因特异表达/甲基化分析的基础;多 genome 布局使每个参考基因组独立
定相、互不串扰。

## 3. 参考文献

1. HiPhase(PacBio)— https://github.com/PacificBiosciences/HiPhase
