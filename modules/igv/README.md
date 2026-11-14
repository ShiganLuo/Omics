# igv — IGV 展示用去重 BAM 与 bigWig 轨道模块

为 IGV 浏览生成展示层文件:samtools markdup 去重 BAM + deepTools
bamCoverage bigWig(`samtools_dedup` → `wig` → `igv_result` 三规则)。

## 1. 模块输入输出与参数

### 1.1 输入

形态 A:`indir` + `bam_substring` 布局参数:

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}{<.bam_substring>}.bam` | 比对 BAM;`bam_substring` 非空时文件名为 `{sample_id}.{bam_substring}.bam`,空则 `{sample_id}.bam` |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.dedup.bam       # markdup 去重 BAM(-r 移除重复)
{outdir}/{sample_id}/{sample_id}.dedup.bam.bai   # 索引
{outdir}/{sample_id}/{sample_id}.bigwig          # 覆盖度 bigWig
{outdir}/{sample_id}/samtools_dedup.sh           # 回执(bamCoverage_wig.sh 同目录)
logdir/{sample_id}/samtools_dedup.log            # 运行日志(wig.log 同目录)
```

`igv_result` 规则无 output,只消费 bigWig 作为 DAG 终点锚点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.samtools` | 可执行文件 | null(回退 `samtools`) | samtools 可执行路径 |
| `Procedure.bamCoverage` | 可执行文件 | null(回退 `bamCoverage`) | deepTools bamCoverage 可执行路径 |
| `bam_substring` | — | "" | 输入 BAM 名中样本名后的后缀(如 `dedup` → `{s}.dedup.bam`) |
| `Params.bamCoverage.binSize` | `--binSize` | 50 | bigWig 分箱大小(bp) |
| `Params.bamCoverage.normalizeUsing` | `--normalizeUsing` | `CPM` | 归一化方式 |
| `Params.bamCoverage.offset` | `--Offset` | null | 只统计 read 起点偏移处碱基 |
| `Params.bamCoverage.extendReads` | `--extendReads` | false | true 加开关;正整数则按该长度延伸 read |
| `indir` | — | `output/hisat2` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

(json 的 `Procedure.deepTools` 键与 .smk 读取的 `Procedure.bamCoverage` /
`Procedure.samtools` 不一致,以 .smk 为准。)

线程/内存经 `rule_res(...)` 配置:samtools_dedup threads 12 / mem 16000,
wig threads 12 / mem 8000。

引擎命令(回执脚本原文):

```bash
samtools sort -n -@ 12 {in.bam} | samtools fixmate -m - - | samtools sort -@ 12 - | samtools markdup -r -@ 12 - {outdir}/{s}/{s}.dedup.bam && samtools index -@ 12 {outdir}/{s}/{s}.dedup.bam
bamCoverage --numberOfProcessors 12 --binSize 50 --normalizeUsing CPM -b {outdir}/{s}/{s}.dedup.bam -o {outdir}/{s}/{s}.bigwig
```

## 2. 算法原理与生物学意义

纯工具封装:去重走 samtools 标准 markdup 流程——按名字排序(`sort -n`)→
`fixmate -m` 补 mate 信息/插入片段标签 → 坐标排序 → `markdup -r` 标记并
移除 PCR/光学重复 → 索引。bigWig 由 bamCoverage 按 50 bp 分箱计覆盖度并
做 CPM 归一化,把不同测序深度的样本拉到同一尺度比较。生物学意义:重复
reads 会虚增覆盖峰,`-r` 移除后覆盖度反映真实分子丰度;bigWig 是 IGV/
UCSC 基因组浏览的轨道格式,用于人工核查比对与表达信号。

## 3. 参考文献

1. samtools — https://www.htslib.org/ ;Danecek P et al. *Twelve years of
   SAMtools and BCFtools*. GigaScience 10(2):giab008, 2021.
   DOI: 10.1093/gigascience/giab008
2. deepTools(bamCoverage)— https://deep-tools.readthedocs.io/ ;
   Ramírez F et al. *deepTools2: a next generation web server for deep
   sequence data analysis*. Nucleic Acids Res 44(W1):W160–W165, 2016.
   DOI: 10.1093/nar/gkw257
