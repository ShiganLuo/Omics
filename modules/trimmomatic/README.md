# trimmomatic — Illumina 读段接头修剪与质量过滤模块

用 Trimmomatic(PE/SE 两种模式)对 FASTQ 做接头切除与质量裁剪,输出干净
读段、未配对读段与 `-summary` 报告,供下游比对使用。纯工具封装模块。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | 输入 FASTQ 目录(上游 config 的 `outdir`);双端为 `{indir}/{sample_id}_1.fq.gz` / `{indir}/{sample_id}_2.fq.gz`,单端为 `{indir}/{sample_id}.single.fq.gz`(均在 indir 根下,无样本子目录) |
| `outdir` / `logdir` (config) | 输出/日志目录 |
| `Params.trimmomatic.adapter_pe` | 双端模式接头 FASTA(ILLUMINACLIP 参考) |
| `Params.trimmomatic.adapter_se` | 单端模式接头 FASTA |
| `Procedure.trimmomatic` | Trimmomatic 可执行路径;以 `.jar` 结尾时自动包 `java -jar` |

### 1.2 输出

```
# trimmomatic_Paired(双端)
{outdir}/{sample_id}/{sample_id}_1.fq.gz                  # 干净 R1(temp)
{outdir}/{sample_id}/{sample_id}_2.fq.gz                  # 干净 R2(temp)
{outdir}/{sample_id}/{sample_id}_1.unpaired.fq.gz         # R1 未配对读段
{outdir}/{sample_id}/{sample_id}_2.unpaired.fq.gz         # R2 未配对读段
{outdir}/{sample_id}/trimmomatic_report.txt               # -summary 报告
{outdir}/{sample_id}/trimmomatic_Paired.sh                # 回执(单条引擎命令)

# trimmomatic_Single(单端)
{outdir}/{sample_id}/{sample_id}.single.fq.gz              # 干净读段(temp)
{outdir}/{sample_id}/trimmomatic_report.txt                # -summary 报告
{outdir}/{sample_id}/trimmomatic_Single.sh                 # 回执
```

日志:`{logdir}/{sample_id}/trimmomatic.txt`;每次运行的带时间戳脚本另存
`RECEIPT_DIR`。

### 1.3 参数

config 位于 `Params.trimmomatic` 与 `Procedure.trimmomatic` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.trimmomatic.adapter_pe` | `ILLUMINACLIP:<file>:2:30:10` | — | 双端接头 FASTA;启用双端模式时必填 |
| `Params.trimmomatic.adapter_se` | `ILLUMINACLIP:<file>:2:30:10` | — | 单端接头 FASTA;启用单端模式时必填 |
| `Procedure.trimmomatic` | 引擎可执行 | `trimmomatic` | `.jar` 结尾自动 `java -jar` 包裹 |
| (threads) | `-threads` | 6 | `rule_res("trimmomatic_*")` 可覆盖 |
| (mem_mb) | — | 8000(PE)/ 4000(SE) | 资源声明 |

环境:`conda: "trimmomatic.yaml"`(bioconda `trimmomatic=0.40` + openjdk,对应
SIF `/data/16T/Database/env/trimmomatic/trimmomatic.sif`)。`Procedure.trimmomatic`
也可指向独立 jar(以 `.jar` 结尾自动 `java -jar`);实测 0.40 与 0.41 jar 的
`-summary` 格式一致(逐行键值),qc_summary 的 `parse_trimmomatic` 通吃。

固定裁剪步骤(不在 config 暴露,写死于规则):

```
ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:80
```

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
trimmomatic PE -threads 6 \
    R1.fq.gz R2.fq.gz \
    -summary trimmomatic_report.txt \
    out_1.fq.gz out_1.unpaired.fq.gz \
    out_2.fq.gz out_2.unpaired.fq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:80
```

## 2. 算法原理与生物学意义

Trimmomatic 对 Illumina 读段做接头与质量修剪:

- **ILLUMINACLIP:adapters.fa:2:30:10**:以接头序列为参考切除接头污染;
  `2` = seed 最大错配数,`30` = palindrome 模式(双端接头互补)简单裁剪
  阈值,`10` = simple 模式(单端/单侧接头)匹配阈值,分值越高要求越严。
- **LEADING:3 / TRAILING:3**:从读段 5'/3' 端逐碱基裁掉质量 < 3 的尾部。
- **SLIDINGWINDOW:4:15**:4 bp 滑动窗口,窗口平均质量低于 15 时从该处
  截断。
- **MINLEN:80**:裁剪后短于 80 bp 的读段直接丢弃。

生物学意义:接头残留会引入假比对/假变异,低质量碱基会拉低比对与变异
调用可信度;修剪后读段更可靠。`trimmomatic_report.txt`(summary)给出
各步骤读段/碱基保留比例,可作数据健全性快检——若接头比例异常高提示
文库插入片段过短或接头选择错误。

## 3. 参考文献

1. Trimmomatic — https://github.com/usadellab/Trimmomatic
2. Bolger AM, Lohse M, Usadel B. *Trimmomatic: a flexible trimmer for
   Illumina sequence data*. Bioinformatics 30(15):2114-2120, 2014.
   DOI: 10.1093/bioinformatics/btu170.
