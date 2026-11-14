# demultiplexer — 随机分子标签修剪 + PCR 去重模块(jla-trim)

用 jla-trim 对双端 FASTQ 做 randommer 修剪与 PCR 去重(单规则
`demultiplex_trim_dedup`);引擎来自 pip 包 `jla-demultiplexer`。

## 1. 模块输入输出与参数

### 1.1 输入

形态 A:`indir` + 样本列表,仅双端:

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}_1.fq.gz` | R1,gzip 压缩 FASTQ |
| `{indir}/{sample_id}/{sample_id}_2.fq.gz` | R2 |

规则先用 Python `gzip` 解压为临时 `.tmp.fq`(jla-trim 不接受压缩输入)。

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_1.fq.gz   # 修剪/去重后 R1(gzip -c 压缩)
{outdir}/{sample_id}/{sample_id}_2.fq.gz   # 修剪/去重后 R2
logdir/{sample_id}/demultiplexer.log       # 运行日志
```

jla-trim 在输入文件旁产出 `{tmp}.trimmed.fastq`,压缩到声明输出后临时文件
被清理。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.demultiplexer.ranmer_len` | `-r` | 10 | read1 3' 端 randommer 长度 |
| `Params.demultiplexer.max_ham` | `-maxHam` | 1 | randommer 匹配允许的最大 Hamming 距离 |
| `indir` | — | `input` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

(本模块 json 的 `Params.demultiplexer` 为空 `{}`,上表默认值取自 .smk 的
`config.get(..., 默认)` 回退;`Procedure.demultiplexer` 出现在 json 但 .smk
未读取,引擎固定为 `jla-trim`。)

线程/内存经 `rule_res("demultiplex_trim_dedup")` 配置:threads 默认 1,
mem_mb 默认 8000。

引擎命令(每样本):

```bash
jla-trim -r1 {sample}_1.tmp.fq -r2 {sample}_2.tmp.fq -r 10 -maxHam 1
gzip -c {sample}_1.tmp.fq.trimmed.fastq > {outdir}/{s}/{s}_1.fq.gz
gzip -c {sample}_2.tmp.fq.trimmed.fastq > {outdir}/{s}/{s}_2.fq.gz
```

## 2. 算法原理与生物学意义

纯工具封装:jla-trim 从 read1 3' 端修剪 randommer + 固定接尾序列(AG),
并按 randommer 标签做 PCR 去重(同一原始分子的扩增拷贝只留一条)。
randommer 与 read 的配对允许 `-maxHam` 以内的 Hamming 距离,吸收合成/测序
错误。生物学意义:randommer 是 UMI 的一种,用于区分原始分子与 PCR 拷贝;
去重后的分子计数消除扩增偏好,是小 RNA 定量与低起始量文库分析的前提。

## 3. 参考文献

1. jla-demultiplexer(提供 `jla-trim`:去重 + randommer 修剪;randommer
   长度建议 10/11)— https://github.com/TimNicholsonShaw/jla-demultiplexer ;
   https://pypi.org/project/jla-demultiplexer/(安装来源见
   `demultiplexer.yaml`)
