# macs3/polygenomes — 多基因组峰调用子模块(MACS3)

按基因组分层的 MACS3 峰调用:每样本同时产出 narrow 与 broad 两套峰
(ChIP-seq/DIP-seq),支持 IP−Input 对照模式,并汇总所有 IP 样本的
cutoff 分析曲线。C 形态子模块,共享父目录 `../macs3.yaml` 环境与
`modules/macs3/bin/plot_cutoff_analysis.py` 脚本。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_substring` | 形态 A 拼接:`{indir}/{sample_id}/{sample_id}.{bam_substring}`,默认后缀 `sorted_markdup.bam`(注意:路径不含 genome 段,按样本取 BAM) |
| `ip_samples` | IP 样本列表(macS3_cutoff_plot 汇总用) |
| `input_samples` | Input 对照样本列表 |
| `sample_ip_input_map` | IP 样本 → Input 样本 ID(或 None);有映射时加 `-c` 对照 BAM `{indir}/{input}/{input}.{bam_substring}`,否则单样本模式 |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}_peaks.narrowPeak            # narrow 峰
{outdir}/{genome}/{sample_id}/{sample_id}_peaks.xls                   # narrow 明细
{outdir}/{genome}/{sample_id}/{sample_id}_cutoff_analysis.txt         # narrow cutoff 分析
{outdir}/{genome}/{sample_id}/{sample_id}_broad_peaks.broadPeak       # broad 峰
{outdir}/{genome}/{sample_id}/{sample_id}_broad_peaks.xls             # broad 明细
{outdir}/{genome}/{sample_id}/{sample_id}_broad_cutoff_analysis.txt   # broad cutoff 分析
{outdir}/{genome}/cutoff_analysis.<fmt>                               # 全 IP 样本 cutoff 曲线(<fmt> 取自 formats,默认 png)
```

日志:`{logdir}/{genome}/{sample_id}/macs3.log`、
`{logdir}/{genome}/../group/macs3/macs3_cutoff_plot.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `logs` | 输入/输出/日志目录 |
| `bam_substring` | 文件名后缀 | `sorted_markdup.bam` | BAM 文件名后缀 |
| `sample_ip_input_map` | `-c` | `{}` | IP→Input 映射;无映射即无对照模式 |
| `ip_samples` / `input_samples` | — | `[]` / `[]` | 样本列表 |
| `Procedure.macs3` | `macs3 callpeak` | `macs3` | MACS3 可执行路径 |
| `Params.macs3.bw` | `--bw` | `200` | 模型带宽(bp) |
| `Params.macs3.pvalue` | `-p` | `1e-5` | 峰调用 p 值阈值 |
| `Params.macs3.genome_size` | `-g` | `mm` | 有效基因组大小(可按 genome 给字典,取 `Params.macs3.genome_size.<genome>`) |
| `Params.macs3.broad_cutoff` | `--broad-cutoff` | `0.1` | broad 峰 q 值阈值 |
| `Params.macs3.formats` | `plot_cutoff_analysis.py -f/--format`(可重复) | `["png"]` | cutoff 汇总图图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集) |
| 固定 | `--seed 2346` / `--cutoff-analysis` | — | 随机种子固定,并输出 cutoff 分析 |
| 线程/内存 | — | `threads=4`,`mem_mb=16000` | `rule_res("macs3_callpeak")` 可覆盖 |

环境:`conda: "../macs3.yaml"`。

直接调用示例(等价引擎命令 + 父目录脚本):

```bash
macs3 callpeak --bw 200 -p 1e-5 -g mm --outdir out --name s1 --seed 2346 \
  -t s1.sorted_markdup.bam [-c input.sorted_markdup.bam] --cutoff-analysis
macs3 callpeak --broad --broad-cutoff 0.1 --bw 200 -p 1e-5 -g mm \
  --outdir out --name s1_broad --seed 2346 -t s1.sorted_markdup.bam --cutoff-analysis
python modules/macs3/bin/plot_cutoff_analysis.py \
  --input-files out/s1_cutoff_analysis.txt --sample-names s1 \
  --output out/cutoff_analysis -f png -f pdf
```

(`--output` 为不带扩展名的输出 base 路径;每种 `-f` 格式各写一份
`<output>.<fmt>`。)

## 2. 算法原理与生物学意义

MACS3(模型化分析 ChIP-Seq)用双峰 Poisson 模型估计局部背景,把 IP 相对
Input 的读段富集判为峰,并做 FDR 校正;narrow 峰适合转录因子/尖锐修饰,
broad 峰(`--broad`,q ≤ `broad_cutoff`)适合组蛋白修饰等宽域信号;
`--cutoff-analysis` 输出不同阈值下峰数/基因组覆盖曲线,是选阈值的实用诊断。
结果怎么读:`.narrowPeak`/`.broadPeak` 为 BED6+3(score、-log10 pvalue、
-qvalue),按 score 排序取头部即高置信峰;`cutoff_analysis.<fmt>` 找曲线拐点
权衡峰数与假阳性;narrow 峰可直接喂下游 deeptools peaks 热图与
peak_te_overlap 等模块。

## 3. 参考文献

1. MACS3 — https://github.com/macs3-project/MACS3 ,
   https://macs3-project.github.io/MACS/docs/callpeak.html
2. Zhang Y et al. *Model-based analysis of ChIP-Seq (MACS)*.
   Genome Biology 9(9):R137, 2008. DOI: 10.1186/gb-2008-9-9-r137
