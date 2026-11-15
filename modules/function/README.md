# function — 功能富集分析模块(GO/KEGG、GSEA、GSVA)

对差异表达结果做功能解读:`function_go_kegg`(上下调基因 GO/KEGG
过表征,clusterProfiler)、`function_gsea`(排序基因集富集,fgsea)、
`function_gsva`(样本级通路活性,自研 Python GSVA)。同目录
`pseudobulk_go/`、`sv_go_kegg/` 为独立子模块(各有 README,不在本
文档范围);polygenomes 变体见 `polygenomes/`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | 差异表达目录(上游 DESeq2 的 `outdir`,默认 `output/diff_expression`):go_kegg 取 `{indir}/{contrast}/TEcount/{contrast}.TEcount_Gene.name.tsv`,gsea 取 `{indir}/{contrast}/TEcount/{contrast}.TEcount_Gene.tsv`(缺失即报错) |
| `counts_dir` (config) | featureCounts 计数目录(gsva 用):`all_paired_featureCounts.tsv` / `all_single_featureCounts.tsv`,按 `paired_samples`/`single_samples` 是否非空决定取哪些 |
| `group_pairs` (config) | `{contrast: {"control_group_name", "experimental_group_name"}}` 对比定义;contrast 不在其中直接 `ValueError` |
| `Params.function.gsea.gmt` (config) | GSEA 基因集 GMT(必需,存在性校验;示例 json 写作 `Params.function.gmt`,以 .smk 为准) |
| `Params.function.gsva.gmt` (config) | GSVA 基因集 GMT(必需) |
| `genome.references.geneIDAnno` (config) | 基因 ID 注释(gsea 用;注意取的是 `genome.references` 下,错误信息/示例 json 写作 `genome.geneIDAnno`) |
| `genome.gtf` (config) | GTF(gsva 的基因 ID → 名转换,必需) |

### 1.2 输出

```
# function_go_kegg(每 contrast;图按 formats 每种格式一份,默认 png)
{outdir}/{contrast}/go_back_to_back.<fmt> / kegg_back_to_back.<fmt>
{outdir}/{contrast}/go_up.csv / go_down.csv / kegg_up.csv / kegg_down.csv
{outdir}/{contrast}/up_genes.txt / down_genes.txt

# function_gsea(每 contrast;图默认 jpeg)
{outdir}/{contrast}/GSEA/TEcount_Gene_GSEA.<fmt> / TEcount_Gene_GSEA.csv
(+ TEcount_Gene_GSEA.rnk 排序列表、GSEA/enrichment_plots/<通路>_enrichment.<fmt>,
脚本内产物)

# function_gsva(全队列;热图默认 png)
{outdir}/gsva/gsva_scores.tsv          # geneset × 样本富集评分
{outdir}/gsva/gsva_heatmap.<fmt>       # 聚类热图
{outdir}/gsva/expression_tpm.tsv       # 中间产物:TPM 标准化矩阵
```

回执:`{outdir}/{contrast}/go_kegg.sh`、`gsea.sh`、`{outdir}/gsva/gsva.sh`;
日志 `{logdir}/../group/function/{contrast}.*.log`、`{logdir}/function/gsva.log`。

### 1.3 参数

config 位于 `Params.function` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.function.species` | `-s`(go-kegg)/ `--title` 文案 | `mouse` | 物种(org.Mm/org.Hs 注释库选择;gsva 图标题) |
| `Params.function.go_kegg.lfc_cut` | `--lfc-cut` | 1 | 上下调基因 log2FC 阈值 |
| `Params.function.go_kegg.p_cut` | `--p-cut` | 0.05 | 显著性阈值 |
| `Params.function.go_kegg.top` | `--top` | 10 | 每侧展示条目数 |
| `Params.function.go_kegg.formats` | `--format`(可重复) | `[png]` | 输出图片格式(多种同时输出);合法值 png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集;R 侧 ggsave 设备集)。go-kegg.r 是 R 脚本,仅受理 ggsave 设备集,pgf/raw/rgba 传入即报错 |
| `Params.function.gsea.gmt` | `-g`(gsea) | — | 必需 |
| `Params.function.gsea.formats` | `--format`(可重复) | `[jpeg]` | GSEA 图输出格式(同上;gsea.r 为 R 脚本,ggsave 设备集) |
| `Params.function.gsva.gmt` | `-g`(gsva) | — | 必需 |
| `Params.function.gsva.formats` | `-f`/`--format`(可重复) | `[png]` | GSVA 热图输出格式(matplotlib 全集) |
| `genome.references.geneIDAnno` | `-a`(gsea) | — | 必需 |
| `genome.gtf` | `normalization.py --gtf` | — | 必需 |
| (threads) | — | 1 | 各规则 `rule_res` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致;`--format` 可重复传,
R 侧脚本仅受理 png pdf svg ps eps tif tiff jpg jpeg):

```bash
Rscript modules/function/bin/go-kegg.r \
    -i out/diff_expression/treat_vs_ctrl/TEcount/treat_vs_ctrl.TEcount_Gene.name.tsv \
    -o out/function/treat_vs_ctrl -s mouse \
    --lfc-cut 1 --p-cut 0.05 --top 10 --format png

Rscript modules/function/bin/gsea.r -m Gene \
    -g geneset.gmt \
    -i out/diff_expression/treat_vs_ctrl/TEcount/treat_vs_ctrl.TEcount_Gene.tsv \
    -o out/function/treat_vs_ctrl -a geneIDAnno -t treat_vs_ctrl --format jpeg

python modules/function/bin/normalization.py \
    -i out/counts/all_paired_featureCounts.tsv \
    -o out/function/gsva/expression_tpm.tsv \
    --method tpm --gtf annotation.gtf --convert-gene-name
python modules/function/bin/gsva.py \
    -i out/function/gsva/expression_tpm.tsv -g geneset.gmt \
    -o out/function/gsva --title "GSVA (mouse)" --format png
```

## 2. 算法原理与生物学意义

- **GO/KEGG(go-kegg.r,clusterProfiler `enrichGO`/`enrichKEGG`)**:按
  `lfc_cut` + `p_cut` 把 DESeq2 结果分为上调/下调基因集,对 GO
  (BP/MF/CC)与 KEGG 通路做超几何检验过表征,BH 校正;`*_up/down.csv`
  为条目表(ID、Description、pvalue、p.adjust、geneID、Count),
  back-to-back 双向条形图上下对照,`up/down_genes.txt` 是两个基因集
  本体。读法:优先看 p.adjust 显著且 Count 不小的条目,方向由上/下
  调侧决定;
- **GSEA(gsea.r,fgsea)**:基因按 DESeq2 统计量排序成 .rnk,fgsea
  计算每基因集的富集沿排序分布的加权 KS 统计量与 NES(标准化富集
  分数)、padj;`TEcount_Gene_GSEA.csv` 每行一个通路(NES、padj),
  顶部通路出富集曲线图。与过表征不同,GSEA 不设基因截断,能捕捉
  协同微弱变化的通路;
- **GSVA(gsva.py,自研 Python 实现)**:TPM 标准化(`normalization.py`)
  后,逐样本对每基因集计算 KS 型富集评分(rankdata 累积统计,符号约定
  与 GSVA/GSEA 一致),输出 geneset × 样本评分矩阵与聚类热图——把
  「基因差异」提升为「样本级通路表型」,可与分组/临床表型并列解读。

生物学意义:富集分析把基因列表翻译成生物学过程/通路语言;GO/KEGG
回答「哪些过程过表征」,GSEA 回答「哪些通路整体偏移」,GSVA 回答
「每个样本的通路活性谱」。注意 gsva 为自研实现(非 R GSVA 包),评分
仅供相对比较。

## 3. 参考文献

1. clusterProfiler —— https://github.com/YuLab-SMU/clusterProfiler
   (Yu G, Wang LG, Han Y, He QY. *clusterProfiler: an R package for
   comparing biological themes among gene clusters*. OMICS 16(5):284-287,
   2012. DOI: 10.1089/omi.2011.0118.)
2. fgsea —— https://github.com/ctlab/fgsea
3. Hänzelmann S, Castelo R, Guinney J. *GSVA: gene set variation analysis
   for microarray and RNA-seq data*. BMC Bioinformatics 14:7, 2013.
   DOI: 10.1186/1471-2105-14-7.(方法概念;本模块为自研 Python 实现)
4. Subramanian A et al. *Gene set enrichment analysis: a knowledge-based
   approach for interpreting genome-wide expression profiles*. PNAS
   102(43):15545-15550, 2005. DOI: 10.1073/pnas.0506580102.
