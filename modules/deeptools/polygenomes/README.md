# deeptools/polygenomes — 多基因组 ChIP/ATAC 信号热图子模块(deeptools heatmap)

按基因组分层的 deepTools 信号热图:先把 IP bigWig 与 Input bigWig 相除得
富集比值 bigWig,再围绕 TSS / peak 中心 / 基因(TE)列表画富集热图。stem 为
`heatmap`。C 形态子模块,共享父目录 `../deeptools_heatmap.yaml` 环境与
`modules/deeptools_heatmap/bin/run_heatmap.py` 脚本(computeMatrix +
plotHeatmap 统一封装)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `bigwig_dir` + 布局参数 | 形态 A 拼接:IP bigWig `{bigwig_dir}/{sample_id}/{sample_id}.bigwig`;对照 Input bigWig 由 `sample_ip_input_map` 映射 `{bigwig_dir}/{<input_sample>}/{<input_sample>}.bigwig`(未映射时无对照输入) |
| `sample_ip_input_map` | sample → 对照样本 ID 字典 |
| `genome.<g>.gtf` / `genome.<g>.te_gtf` | 注释 GTF / TE GTF(genes 模式下 `regions.gtf="te"` 用 TE GTF)。注意本模块直接读 `config['genome'][<g>][...]` 平铺键(与其他模块的 `genome.references` 不同) |
| `genome.<g>.tss_bed` | 可选:TSS BED;tss 模式下存在即用之,否则从 GTF 推 TSS |
| peaks 模式输入 | `{indir}/{sample_id}/{sample_id}_peaks.narrowPeak`(MACS3 等峰调用输出) |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}_IP_over_Input.bigwig        # IP/Input 比值 bigWig
{outdir}/{genome}/{sample_id}/{sample_id}_tss_heatmap.{fmt}           # tss 模式
{outdir}/{genome}/{sample_id}/{sample_id}_peaks_heatmap.{fmt}         # peaks 模式
{outdir}/{genome}/{sample_id}/{sample_id}_genes_heatmap.{fmt}         # genes 模式
{outdir}/{genome}/{sample_id}/{sample_id}_{gene_name}_heatmap.{fmt}   # 单基因/TE 热图(heatmap_gene)
```

文件名后缀由 `Params.computeMatrix.regions` 决定(`tss`/`peaks`/`genes`,
未识别时 `_heatmap`);`{fmt}` 取自 `Params.plotHeatmap.formats`(每种格式
一个文件,默认 `png`,与旧版文件名一致)。日志:
`{logdir}/{genome}/{sample_id}/bigwig_ratio.log`、`heatmap.log`、
`heatmap_{gene_name}.log`。

### 1.3 参数

config 位于 `Params.computeMatrix`(矩阵)与 `Params.plotHeatmap`(出图):

| 参数 | CLI(run_heatmap.py / bigwigCompare) | 默认 | 说明 |
|---|---|---|---|
| `bigwig_dir` | — | `""` | bigWig 根目录 |
| `sample_ip_input_map` | `-b2` 选取 | `{}` | IP→Input 样本映射 |
| `Params.computeMatrix.samples` | — | `samples` 键 | 顶层读取,当前规则未消费 |
| `Params.computeMatrix.regions` | 区域来源 | `"tss"` | `"tss"`/`"peaks"`/`{"genes":[...], "match_by":"gene_name", "gtf":"te"}` |
| `Params.computeMatrix.mode` | `--mode` | `reference-point` | 模式自动覆盖:tss/peaks→`reference-point`,genes→`scale-regions` |
| `Params.computeMatrix.referencePoint` | `--reference-point` | `center` | 自动覆盖:tss→`TSS`,peaks→`center`,genes→`TSS` |
| `Params.computeMatrix.before` / `after` | `--before` / `--after` | `3000` / `3000` | reference-point 模式上下(bp) |
| `Params.computeMatrix.upstream` / `downstream` | `--upstream` / `--downstream` | `3000` / `3000` | scale-regions 模式上下游(bp) |
| `Params.computeMatrix.bodyLength` | `--body-length` | `5000` | scale-regions 基因体长度(bp) |
| `Params.computeMatrix.binSize` | `--bin-size` | `10` | 分箱(bp) |
| `Params.computeMatrix.sortUsing` | `--sort-using` | `mean` | 行排序依据 |
| `Params.computeMatrix.top_n` | `--top-n` | `0` | 取前 N 行(0=全部) |
| `Params.computeMatrix.tss_flank` | `--tss-flank` | `1000` | GTF 推 TSS 时的侧翼(bp) |
| `Params.plotHeatmap.colorMap` | `--color-map` | `YlOrRd` | 色图 |
| `Params.plotHeatmap.heatmapHeight` / `heatmapWidth` | `--height` / `--width` | `15` / `8` | 图高宽(inch) |
| `Params.plotHeatmap.whatToShow` | `--what-to-show` | `plot, heatmap and colorbar` | 显示元素 |
| `Params.plotHeatmap.formats` | `-f/--format` → `--plotFileFormat` | `["png"]` | 图片格式列表,每种格式一个输出文件;合法值 `png` `eps` `pdf` `svg` `plotly`(deeptools plotHeatmap `--plotFileFormat` 工具集,解析期校验 fail fast) |
| 固定 | `bigwigCompare --operation ratio --pseudocount 1 -p` | — | 比值 bigWig 一律 ratio + 伪计数 1 |

`match_by`(genes 模式 `--match-by`,默认 `gene_name`)与 `gene-names`
(可重复)由 `regions` 字典控制;非 TE 的 genes 模式附加 `--merge`
(外显子合并为基因体),TE 位点不合并(每位点独立成行)。
环境:`conda: "../deeptools_heatmap.yaml"`。

直接调用示例(等价引擎命令):

```bash
bigwigCompare -b1 s1.bigwig -b2 input.bigwig --operation ratio --pseudocount 1 \
  -o s1_IP_over_Input.bigwig -p 4
python3 modules/deeptools_heatmap/bin/run_heatmap.py \
  --ip-bigwig s1_IP_over_Input.bigwig --output s1_tss_heatmap \
  -f png -f pdf \
  --mode reference-point --reference-point TSS \
  --before 3000 --after 3000 --bin-size 10 --sort-using mean \
  --regions tss.bed --threads 4 --title s1 \
  --color-map YlOrRd --height 15 --width 8 \
  --what-to-show "plot, heatmap and colorbar"
```

## 2. 算法原理与生物学意义

deepTools `bigwigCompare --operation ratio` 先把 IP 信号除以 Input(伪计数 1
抑噪),消除文库深度与开放性背景;`run_heatmap.py` 内部调 computeMatrix +
plotHeatmap:在 TSS/peak/基因体周围取信号矩阵,按均值排序画富集热图。
生物学读法:TSS 模式看信号是否在转录起始位点富集(启动子占位);peaks 模式
看峰中心对称富集(峰调用可信度的正交验证);genes/TE 模式看基因体或 TE
亚家族的信号分布差异(如衰老样本 TE 开放性变化)。IP/Input 比值 > 1 的行
即富集区域;`--top-n` 可聚焦最强位点。

## 3. 参考文献

1. deepTools — https://deeptools.readthedocs.io/ ,
   https://github.com/deeptools/deeptools
2. Ramírez F et al. *deepTools2: a flexible web suite for exploring
   genomic data*. Nucleic Acids Research 44(W1):W160–W165, 2016.
   DOI: 10.1093/nar/gkw257
