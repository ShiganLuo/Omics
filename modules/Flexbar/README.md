# Flexbar — 按 barcode 拆分混样 FASTQ 模块

用 flexbar 把带 barcode 的混样 FASTQ 按样本拆分(`flexbar_demultiplex` 单规则)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `samples` | 形态 A:`{indir}/{sample_id}.fq.gz`,gzip 压缩的单端 FASTQ(read 5' 端带 barcode),由 input 函数拼接 |
| `Params.flexbar.barcode_file` | barcode 序列文件,FASTA 格式(每条记录名 = barcode ID,序列 = barcode),传给 flexbar `-b` |

### 1.2 输出

```
{outdir}/{sample_id}.demux.fq.gz                 # 拆分后 reads(pigz 压缩;声明为 temp)
{outdir}/{sample_id}/flexbar_demultiplex.sh      # 回执(单条引擎命令)
logdir/{sample_id}/flexbar_demux.txt             # 声明的日志输出
logdir/{sample_id}/flexbar_demux_run.txt         # 实际运行日志(log 定向,shell 重定向追加)
RECEIPT_DIR/flexbar_demultiplex_{sample}_{时间戳}.sh  # 同回执(带时间戳副本)
```

flexbar 以 `-t {outdir}/{sample_id}.demux` 为前缀产拆分文件
`{sample_id}.demux_1.fastq`,再由 `pigz -c` 压入 `{sample_id}.demux.fq.gz`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.flexbar` | 可执行文件 | null(回退 `flexbar`) | flexbar 可执行路径 |
| `Params.flexbar.barcode_file` | `-b` | `barcodes.fasta` | barcode FASTA 文件 |
| `Params.flexbar.barcode_min_overlap` | `-n` | 1 | barcode 与 read 的最小重叠碱基数 |
| `Params.flexbar.barcode_allow_mismatch` | `-m` | 0 | barcode 匹配允许的错配数 |
| `indir` | — | `output/raw_fastq` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

线程/内存经 `rule_res("flexbar_demultiplex")` 配置:threads 默认 4,mem_mb 默认 8000。

引擎命令(每个样本,回执脚本原文):

```bash
flexbar -r {indir}/{sample_id}.fq.gz -b barcodes.fasta \
    -n 1 -m 0 \
    -t {outdir}/{sample_id}.demux
pigz -c {outdir}/{sample_id}.demux_1.fastq > {outdir}/{sample_id}.demux.fq.gz
```

## 2. 算法原理与生物学意义

纯工具封装:flexbar 是短读长预处理工具(barcode 拆分 + 接头/质量修剪)。
本模块只用其 barcode demultiplex 能力:flexbar 把每条 read 5' 端与 FASTA
中的 barcode 序列比对,重叠达到 `-n` 且错配不超过 `-m` 时归入该 barcode
对应的拆分文件;容忍少量错配可吸收测序错误,避免样本间 reads 误分配。
生物学意义:多样本混池测序(multiplexing)后按样本还原 FASTQ,是按样本
开展下游分析的前提。

## 3. 参考文献

1. flexbar — https://github.com/seqan/flexbar
2. Dodt M, Roehr JT Jr, Ahmed R, Dieterich C. *FLEXBAR—Flexible Barcode and
   Adapter Processing for Next-Generation Sequencing Platforms*.
   Biology 1(3):895–905, 2012. DOI: 10.3390/biology1030895
