# fragment_size — 文库插入片段长度分布模块(自研脚本)

用 `samtools stats` 汇总每样本插入片段(insert size)分布,自研脚本
合并为多样本直方表并绘图(167 bp 单核小体红线标注),用于文库构建
质量评估。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | BAM 目录(上游 config 的 `outdir`),形态 A 布局:`{indir}/{sample_id}/{sample_id}{bam_suffix}.bam` |
| `bam_suffix` (config,顶层) | BAM 文件名后缀段;缺省空 |
| `samples` (config) | 样本 ID 列表 |
| `ROOT_DIR` (config) | 项目根目录,定位 `modules/fragment_size/bin/` 脚本 |
| `outdir` / `logdir` / `logdir_combine` (config) | 输出/日志(组合级日志走 `logdir_combine`) |

### 1.2 输出

```
{outdir}/stats/{sample_id}/{sample_id}_stats.txt     # samtools stats 全量输出
{outdir}/fragment/FragmentSize.txt                   # 插入长度 × 样本 计数表(TSV)
{outdir}/fragment/FragmentSize.<fmt>                 # 分布曲线图(<fmt> 取自 formats,默认 png)
```

回执:`{outdir}/stats/{sample_id}/samtools_stats.sh`、
`{outdir}/fragment/getFragmentSize.sh`、`plotFragmentSize.sh`。
`fragment_size_result` 为聚合规则。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `bam_suffix` | — | "" | 见输入 |
| `Procedure.samtools` | `samtools stats` | `samtools` | 统计引擎 |
| `Params.fragment_size.formats` | `-f/--format`(可重复) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集) |
| (threads) | — | — | 未设线程参数 |
| (mem_mb) | — | 8000 | 各规则 `rule_res` 可覆盖 |

自研脚本参数(`getFragmentSize.py` / `plotFragmentSize.py`):

| 参数 | 默认 | 说明 |
|---|---|---|
| `--input`(get,可多个) | 必填 | samtools stats 文件列表 |
| `--out`(get) | 必填 | 直方表输出路径 |
| `--input`(plot) | 必填 | 直方表路径 |
| `--output`(plot) | 必填 | 图片输出 base 路径(不带扩展名;实际文件为 `<output>.<format>`) |
| `-f/--format`(plot,可重复) | png | 输出格式(png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba) |
| `--highlight` | 167 | 高亮竖线的片段长度(bp) |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
samtools stats input/s1/s1.sorted.bam > out/stats/s1/s1_stats.txt

python modules/fragment_size/bin/getFragmentSize.py \
    --input out/stats/s1/s1_stats.txt out/stats/s2/s2_stats.txt \
    --out out/fragment/FragmentSize.txt

python modules/fragment_size/bin/plotFragmentSize.py \
    --input out/fragment/FragmentSize.txt \
    --output out/fragment/FragmentSize --highlight 167 -f png -f pdf
```

## 2. 算法原理与生物学意义

- `samtools stats` 输出中 `IS` 行给出插入片段长度分布
  (`IS\t<insert size>\t<pairs total>\t<inward>\t<outward>\t<other>`);
- `getFragmentSize.py` 流式读全部 stats 文件,取每条 `IS` 行的插入长度
  与计数字段(样本名取 stats 文件名去 `.txt`),pivot 成
  「插入长度 × 样本」计数表 `FragmentSize.txt`;
- `plotFragmentSize.py` 把每样本计数列归一为频率(除以总计数)后画
  分布曲线,并在 `--highlight`(默认 167)处画红色虚线。

生物学意义:插入片段分布是文库构建质量的直接读数——DNA 文库应呈
单峰窄分布(片段化/size selection 是否均一),峰位偏小提示降解或
过度片段化。167 bp ≈ 单核小体缠绕 DNA 长度,cfDNA/ATAC 类数据的
单核小体峰可作核小体定位健康的参照(多聚核小体梯 ~330/500 bp);
多样本曲线叠画便于横向比较批次差异。

## 3. 参考文献

1. getFragmentSize.py / plotFragmentSize.py 为本仓库自研脚本;
   统计引擎 samtools —— https://github.com/samtools/samtools
   (Li H et al. *The Sequence Alignment/Map format and SAMtools*.
   Bioinformatics 25(16):2078-2079, 2009. DOI: 10.1093/bioinformatics/btp352.)
