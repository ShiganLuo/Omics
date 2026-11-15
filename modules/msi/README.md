# msi — 微卫星不稳定性(MSI)检测模块(msisensor-pro)

肿瘤-配对正常 MSI 分析:`msi_scan` → `msi_call` → `msi_plot` →
`msi_result`;绘图脚本为模块自带 `bin/bar_plot.r`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}{<.bam_substring>}.bam` + `.bai` | 每样本 BAM/索引;`bam_substring` 非空时为 `{s}.{sub}.bam`,索引为同 stem `.bai` |
| `pairs`(config) | 肿瘤样本 → 配对正常样本 ID 字典;键即待检样本列表(`samples = list(pairs)`) |
| `genome.fasta` / `genome.fai` | 参考基因组及索引(msi_scan 扫描用) |

### 1.2 输出

```
{outdir}/msi_scan/msisensor_scan.tsv          # 全基因组微卫星位点表(每 genome 一次,共享)
{outdir}/{sample_id}/{sample_id}_msi_pro.tsv  # 每肿瘤样本 MSI 结果(vs 配对正常)
{outdir}/msi_input.txt                        # 汇总表(sample<TAB>score,每数据行取第 3 列)
{outdir}/msi_barplot.<fmt>                    # 跨样本 MSI 柱状图(每种 formats 一个文件;默认 msi_barplot.jpg)
{outdir}/msi_scan.sh 等                       # 各规则回执
logdir_combine/msi/msi_scan.log               # scan/plot 日志(默认 logdir_combine=logdir)
logdir/{sample_id}/msi_call.log               # 每样本 call 日志
```

`msi_result` 规则无 output,消费上述产物作为 DAG 终点锚点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.msisensorpro` | 可执行文件 | null(回退 `msisensorpro`) | msisensor-pro 路径 |
| `Procedure.Rscript` | 可执行文件 | null(回退 `Rscript`) | Rscript 路径(msi_plot) |
| `Params.msi.formats` | `bar_plot.r --format`(可重复) | `["jpg"]` | 柱状图输出格式,每种格式一个文件 `msi_barplot.<fmt>`;合法值 png pdf svg ps eps tif tiff jpg jpeg(R 设备集;pgf/raw/rgba 为 matplotlib 专属,传入即报错) |
| `pairs` | — | {} | 肿瘤 → 配对正常样本字典 |
| `bam_substring` | — | "" | 输入 BAM 名后缀 |
| `genome.fasta` / `genome.fai` | `scan -d` | — | 参考基因组与 .fai |
| `logdir_combine` | — | `logdir` | scan/plot 日志目录 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录基路径 |

线程/内存经 `rule_res(...)` 配置:msi_scan 1/16000,msi_call 4/16000,
msi_plot 1/4000。

引擎命令(回执脚本原文):

```bash
msisensorpro scan -d genome.fa -o msi_scan/msisensor_scan.tsv
msisensorpro msi -d msi_scan/msisensor_scan.tsv -n {normal}.bam -t {tumor}.bam -o {s}_msi_pro.tsv
Rscript modules/msi/bin/bar_plot.r -i msi_input.txt -o msi_barplot --format jpg
```

直接调用绘图脚本(更多可选项:`-s` 样本列名、`-c` 分数列名、`-t` 分隔符、
`-x/-y` 轴标签、`-f/--format` 可重复出多格式;`-o` 为不带扩展名的输出前缀,
每种格式写 `<前缀>.<fmt>`,传入 `out.jpg` 这类旧式全名会自动去掉扩展名):

```bash
Rscript bin/bar_plot.r -i msi_input.txt -o out --format jpg --format pdf -y MSI
```

## 2. 算法原理与生物学意义

msisensor-pro 分两步:`scan` 在参考基因组上枚举微卫星(短串联重复)
位点建立位点表;`msi` 对每个位点比较肿瘤与配对正常 BAM 的重复长度分布,
聚合出样本级 MSI 分数(输出第 3 列,`msi_plot` 汇总为
`sample<TAB>score`)。生物学意义:错配修复缺陷(dMMR)导致微卫星长度
不稳定累积,MSI 是结直肠癌等肿瘤的重要标志物(dMMR/MSI-H 与免疫治疗
获益相关);柱状图给出队列内各样本 MSI 水平概览。

## 3. 参考文献

1. MSIsensor-pro — https://github.com/xjtu-omics/MSIsensor-pro
