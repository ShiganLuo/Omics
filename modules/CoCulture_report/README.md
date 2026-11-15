# CoCulture_report — 共培养双基因组基因/TE 表达对比报告模块

对共培养(co-culture)双基因组(如人/鼠、宿主/共生)的 TEtranscripts
计数表出对比报告:每组过滤、取 top 特征、按 TE 分类层级画基因与 TE
双条形图。自研脚本 `bin/HRT.py`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | 双基因组计数目录(上游 TEtranscripts 的 `outdir`):`{indir}/{genome1}/all_TEcount_name.tsv` 与 `{indir}/{genome2}/all_TEcount_name.tsv` |
| `genome_pairs` (config) | 二元列表 `[genome1, genome2]`(基因组名,进路径与图例) |
| `genome1_samples` / `genome2_samples` (config) | 两基因组各自的样本名列表(用于分组取列;示例 config 写作 `{}`,应传列表) |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 输出/日志/项目根 |

### 1.2 输出

```
{outdir}/CoCulture_report/CoCulture_report.html   # .smk 声明的报告输出
{outdir}/CoCulture_report/{genome1}_2_{genome2}.{fmt}   # 双条形图(按 formats 每种格式一个文件)
{outdir}/CoCulture_report/{genome2}_2_{genome1}.{fmt}   # 同上(另一方向;仅对应样本列表非空时产出)
{outdir}/HRT.sh                                   # 回执
```

日志:`{logdir}/CoCulture_report/CoCulture_report.log`。

已知问题(接线前对齐):
- .smk 中脚本路径为 `modules/CoCulture_report/HRT.py`,实际脚本在
  `modules/CoCulture_report/bin/HRT.py`;
- .smk 传参(`--genome1_count/--genome2_count/--genome1/--genome2/
  --output`)与 `bin/HRT.py` 的 argparse(`-i1/-i2/-g1/-g2/-o` 等)不
  一致;且 HRT.py 的 `-o` 是输出目录(写图与中间表),与声明的 HTML
  输出不同。以下参数表按 `bin/HRT.py` 真实 argparse 记录。

### 1.3 参数(bin/HRT.py 真实 CLI)

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| (genome1 计数表) | `-i1` / `--infile1` | 必填 | TEtranscripts 输出,基因组 1 |
| (genome2 计数表) | `-i2` / `--infile2` | 必填 | TEtranscripts 输出,基因组 2 |
| `genome1_samples` | `-g1` / `--genome1_samples`(nargs +) | — | 基因组 1 样本名 |
| `genome2_samples` | `-g2` / `--genome2_samples`(nargs +) | — | 基因组 2 样本名 |
| (输出目录) | `-o` / `--outdir` | 必填 | 图与中间表输出目录 |
| — | `-mc` / `--min_count` | 5 | 特征最小计数过滤 |
| — | `--top` | 10 | 每组展示的 top 特征数 |
| — | `-m` / `--mode` | `TEcount` | 输入格式:`TEcount`(family:class:superfamily)或 `TElocal`(多一段 location) |
| — | `-l` / `--level` | `subfamily` | TE 分类层级(loci/subfamily/family/class) |
| `formats` | `-f` / `--format`(可重复) | `png` | 输出图格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集);非法值直接报错并列出合法值 |
| (mem_mb) | — | 8000 | `rule_res("HRT")` 可覆盖 |

直接调用示例(按脚本真实参数):

```bash
python modules/CoCulture_report/bin/HRT.py \
    -i1 input/genome1/all_TEcount_name.tsv \
    -i2 input/genome2/all_TEcount_name.tsv \
    -g1 s1 s2 -g2 s3 s4 \
    -o out/CoCulture_report \
    -mc 5 --top 10 -m TEcount -l subfamily \
    -f png -f pdf
```

## 2. 算法原理与生物学意义

HRT.py 读两基因组的 TEtranscripts 基因+TE 计数表,按样本分组
(genome1 vs genome2):每组按 `min_count` 过滤低计数特征、取计数最高的
`top` 特征,按 `level` 把 TE 特征名归并到指定分类层级(TEcount 特征名
形如 family:class:superfamily,TElocal 多一段 locus,由 `mode` 决定解析
方式),绘制基因与 TE 的双条形图并落中间表。

生物学意义:共培养体系里两个物种的表达互相影响(宿主应答、共生调控、
异种移植污染评估);基因与 TE 并排对比可观察 TE 去抑制等转录组重编程
是否物种特异。`level` 选择决定解读粒度——subfamily/family 级稳健,
loci 级更细但计数更稀疏。

## 3. 参考文献

1. HRT.py 为本仓库自研脚本(无外部文献);输入计数表由 TEtranscripts
   产出(https://github.com/mhammell-laboratory/TEtranscripts)。
