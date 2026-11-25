# DESeq2/polygenomes — 多基因组差异表达分析子模块(TEcount/TElocal + TE 位点图)

对 TEtranscripts/TElocal 定量矩阵做 DESeq2 差异分析:每个 (genome, contrast)
产出基因/TE/Gene_TE 三层差异表、PCA、火山图、热图,并在 TElocal 结果上出
Manhattan、位点 boxplot、亚家族 bar 三类 TE 位点图。C 形态子模块,共享父
目录 `../DESeq2.yaml`(分析)与 `../tefig.yaml`(出图)环境,以及父目录
`modules/DESeq2/bin/` 脚本(DESeq2.r、write_group_tsv.py、gene_id2name.py、
te_manhattan.py、te_locus_boxplot.py、te_subfamily_bar.py)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 布局参数 | 形态 A 拼接,计数矩阵 `{indir}/{genome}/{genome}_TEcount.tsv`、`{indir}/{genome}/{genome}_TElocal.tsv`(TEtranscripts/TElocal 定量产出) |
| `comparisons` | 顶层 config 键(注意:json 示例把对比定义放在 `sample_params.comparisons`,本 .smk 从顶层 `comparisons` 读)。结构 `{genome: {contrast: {control_group_name, experimental_group_name, control_samples, experimental_samples}}}` |
| `genome.references.<g>.geneIDAnno` | 基因 ID 注释 CSV(`-a`);.smk 读 `genome.references`(复数),json 示例平铺 `genome.<g>.*`,以 .smk 为准 |
| `genome.references.<g>.gtf` / `.gene_map` | 注释 GTF 与 gene_map pkl(`gene_id2name.py` 用),缺失报错 |
| `genome.references.<g>.TEind` | TElocal 位点索引 `.locInd`(TE 位点图的坐标来源),缺失报错 |

### 1.2 输出(每 genome/contrast;`{mode}` = TEcount 或 TElocal)

```
{outdir}/{genome}/{contrast}/{mode}/DESeq2.done
{outdir}/{genome}/{contrast}/{mode}/group.tsv
{outdir}/{genome}/{contrast}/{mode}/PCA/{contrast}.cpmPCA.{fmt}
{outdir}/{genome}/{contrast}/{mode}/volcano/{contrast}.{mode}_Gene_volcano.{fmt}
{outdir}/{genome}/{contrast}/{mode}/heatmap/{contrast}.{mode}_Gene_updown.{fmt}
{outdir}/{genome}/{contrast}/{mode}/upDown/{contrast}.{mode}_{Gene,TE,Gene_TE}_updown.tsv
{outdir}/{genome}/{contrast}/{mode}/{contrast}.{mode}_{Gene,TE,Gene_TE}.tsv
{outdir}/{genome}/{contrast}/{mode}/{contrast}.{mode}_{Gene,TE,Gene_TE}.name.tsv
```

TElocal 专属图：

```
{outdir}/{genome}/{contrast}/TElocal/manhattan/{contrast}.manhattan.{fmt}
{outdir}/{genome}/{contrast}/TElocal/locus_box/{contrast}.locus_box.{fmt}
{outdir}/{genome}/{contrast}/TElocal/subfamily_bar/{contrast}.subfamily_bar.{fmt}
```

`{fmt}` 由 `formats` 参数展开（每种格式一个文件；默认 `png`，与历史一致）。

日志:`{logdir}/DESeq2/{genome}/{contrast}.log`、`.{manhattan|locus_box|subfamily_bar}.log`。

### 1.3 参数

对比定义(顶层 `comparisons.<genome>.<contrast>`):`control_group_name`、
`experimental_group_name`、`control_samples`、`experimental_samples`
(任一为空即报错)。

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.DESeq2.manhattan.label_top` | `--label-top` | `10` | 标注前 N 位点 |
| `Params.DESeq2.manhattan.ylabel` | `--ylabel` | `log2(1+\|delta TPM\|)` | y 轴标签 |
| `Params.DESeq2.manhattan.highlight_families` | `--highlight-family`(可重复) | `{}` | 高亮 TE 家族,按 `[genome][contrast]` 给列表 |
| `Params.DESeq2.locus_boxplot.families` | `--family`(可重复) | 必需 | 目标 TE 家族(`[genome][contrast]`,为空报错) |
| `Params.DESeq2.locus_boxplot.match_field` / `match_mode` | `--match-field` / `--match-mode` | `subfamily` / `exact` | 家族匹配 |
| `Params.DESeq2.locus_boxplot.normalize` / `kind` | `--normalize` / `--kind` | `log2tpm` / `box` | 归一化与图形 |
| `Params.DESeq2.subfamily_bar.value` / `normalize` | `--value` / `--normalize` | `total` / `tpm` | 统计口径 |
| `Params.DESeq2.subfamily_bar.match_field` / `match_mode` | 同上 | `subfamily` / `exact` | 家族匹配 |
| `Params.DESeq2.subfamily_bar.yscale` / `ylabel` | `--yscale` / `--ylabel` | `linear` / `TPM` | 轴 |
| `Params.DESeq2.subfamily_bar.families` | `--family` | 必需 | 同 locus_boxplot |
| `Params.DESeq2.manhattan.formats` / `locus_boxplot.formats` / `subfamily_bar.formats` | `--format`（可重复） | `["png"]` | TE 位点图出图格式（每种格式一个文件 `<base>.<fmt>`）；缺省回退 `Params.DESeq2.formats`；合法值 matplotlib 全集 png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba |
| `Params.DESeq2.formats` | DESeq2.r `--format`（可重复） | `["png"]` | R 图（cpmPCA/火山/热图）格式；R 侧合法值仅 png pdf svg ps eps tif tiff jpg jpeg，pgf/raw/rgba 传入即报错 |
| 线程/内存 | — | 分析 `threads=4, mem_mb=16000`;出图 `threads=2, mem_mb=8000` | `rule_res(...)` 可覆盖 |

图参数解析顺序(.smk 内置 `fig_param`):`Params.DESeq2.<rule>.<key>` →
`Params.DESeq2.<key>` → 内置默认;家族列表按 `[genome][contrast]` 取。
固定调用:DESeq2.r `-m <mode> -i <matrix> -g group.tsv -p <ctrl> <exp>
-f heatmap volcano pca -o <dir> -a <geneIDAnno> -Tcm all -px <contrast>
--format <fmt>...`(每种 formats 一个 `--format`);Manhattan 固定
`--signif-col padj --signif-thresh 0.05 --value-thresh 0
--burst-sign up --y-transform log2 --y-break-auto --hline 0
--format <fmt>...`。

直接调用示例(等价引擎命令):

```bash
python modules/DESeq2/bin/write_group_tsv.py \
  -o group.tsv -c s1,s2 -t s3,s4 -p youth -e aged
Rscript modules/DESeq2/bin/DESeq2.r \
  -m TEcount -i counts/GRCm39_TEcount.tsv -g group.tsv -p youth aged \
  -f heatmap volcano pca -o out/GRCm39/aged_vs_youth/TEcount \
  -a geneIDAnnotation.csv -Tcm all -px aged_vs_youth --format png
python modules/DESeq2/bin/gene_id2name.py \
  -i aged_vs_youth.TEcount_Gene.tsv aged_vs_youth.TEcount_TE.tsv \
     aged_vs_youth.TEcount_Gene_TE.tsv \
  -g annotation.gtf -m gene_map.pkl
python modules/DESeq2/bin/te_locus_boxplot.py \
  --matrix counts/GRCm39_TElocal.tsv --annotation group.tsv \
  --match-field subfamily --match-mode exact --normalize log2tpm --kind box \
  --group-order youth,aged --loci-tsv TE.gtf.locInd \
  -o out/.../locus_box --out-prefix aged_vs_youth.locus_box --format png --family ERVL
```

## 2. 算法原理与生物学意义

DESeq2 用负二项广义线性模型 + 收缩估计做 RNA-seq 差异表达,输出
log2FC 与 BH 校正 padj;本模块按 `comparisons` 的 control/experimental
对比三层定量矩阵分别分析:Gene(基因表达)、TE(TE 亚家族表达)、
Gene_TE(基因-TE 融合计数,TEtranscripts 口径),并把基因 ID 注释为
gene name(`gene_id2name.py`)。TE 位点图基于 TElocal 逐位点结果:
Manhattan(`--burst-sign up` 只看上调爆发位点,log2 变换 + 断轴)找显著
上调的 TE 位点;`locus_boxplot` 按家族画逐位点表达分布;`subfamily_bar`
画亚家族总量对比——用于衰老/肿瘤中 TE 去抑制(如 ERVL/ERV 复活)的
可视化验证。结果怎么读:`upDown/*.tsv` 为差异基因/TE(看 log2FC、padj);
火山图上点=显著基因,横轴 log2FC;Manhattan 每点一个 TE 位点,y 为效应
量,标注 top 位点可直接进实验验证。

## 3. 参考文献

1. DESeq2 — https://bioconductor.org/packages/DESeq2
   (Love MI, Huber W, Anders S. *Moderated estimation of fold change and
   dispersion for RNA-seq data with DESeq2*. Genome Biology 15:550, 2014.
   DOI: 10.1186/s13059-014-0550-8)
2. TEtranscripts / TElocal — https://github.com/mhammell-laboratory/TEtranscripts
   (Jin Y et al., *TEtranscripts: differential expression of transposable
   elements*, 2015;TElocal 为同组逐位点实现)
