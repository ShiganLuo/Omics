# enrichment — GO/KEGG 富集分析模块(引擎脚本;规则待补)

GO 与 KEGG 过表达富集分析 + 柱状图(gseapy/Enrichr)。注意:`enrichment.smk`
目前是占位(`# need test`,无任何规则),模块尚不能进 DAG;可用部分是
`bin/` 下的引擎脚本(`go-kegg.py` Python 版、`go-kegg_back.r` R 备用实现)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| 基因列表(gene symbols) | `enrich_go_kegg_analysis(gene_list, ...)`,来自差异基因/差异峰关联基因表(如 `gene_name` 列去重) |
| `organism` | `'human'` / `'mouse'`(函数参数,默认 human) |

config 骨架(`enrichment.json`)仅有 `indir`/`outdir`/`logdir`/`ROOT_DIR`
与空的 `Params`/`Procedure` 段——规则未实现,暂无消费。

### 1.2 输出

```
{outdir}/GO_Biological_Process_2021_enrichment.csv   # GO BP 富集表
{outdir}/GO_Cellular_Component_2021_enrichment.csv   # GO CC 富集表
{outdir}/GO_Molecular_Function_2021_enrichment.csv   # GO MF 富集表
{outdir}/KEGG_enrichment.csv                          # KEGG 富集表
{outdir}/{key}_top10.png                              # Top10 富集柱状图(每库一张)
```

### 1.3 参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `gene_list` | — | 输入基因 symbol 列表 |
| `organism` | `human` | 物种;mouse 时 KEGG 库用 `KEGG_2019_Mouse`,human 用 `KEGG_2021_Human` |
| `outdir` | `enrichment_results` | 输出目录 |
| `top_n`(绘图) | 10 | 柱状图取 Adjusted P-value 最小前 N 条 |
| `kind`(绘图) | — | `GO`/`KEGG` 决定配色与轴标签 |
| 富集显著性截断 | 0.05 | `gseapy.enrichr(cutoff=0.05)` |

直接调用示例(`bin/go-kegg.py`;当前 `__main__` 为硬编码示例路径,无 CLI
参数解析——直接运行会读写脚本内写死的路径,建议 import 函数调用):

```python
from go_kegg import enrich_go_kegg_analysis, plot_enrichment_bar  # bin/go-kegg.py
results = enrich_go_kegg_analysis(gene_list, organism="mouse", outdir="out")
for key, df in results.items():
    plot_enrichment_bar(df, f"out/{key}_top10.png", top_n=10,
                        kind="KEGG" if key == "KEGG" else "GO")
```

## 2. 算法原理与生物学意义

过表达富集(ORA):对输入基因集,按超几何检验评估每个 GO term / KEGG
通路的基因重叠是否显著超出随机期望,Enrichr 后端返回 Adjusted P-value
(BH 校正 FDR)与 Overlap(命中数)。脚本固定用
`GO_Biological_Process_2021` / `GO_Cellular_Component_2021` /
`GO_Molecular_Function_2021` / `KEGG_2021_Human`(或 `KEGG_2019_Mouse`)
四库。柱状图取 Adjusted P-value 最小前 N 条,画 −log10(FDR) 横条
(FDR<0.01 深色),条上标注 n=命中数。生物学意义:把差异基因列表转译为
功能/通路层面的结论,是组学分析的收尾注释步。结果阅读:按 Adjusted
P-value 排序看 Top term/通路,Overlap 大且 FDR 小者为稳健结论。

## 3. 参考文献

1. gseapy —— https://github.com/zqfang/gseapy(Enrichr 接口封装);
   Fang Z, Liu X, Peltz G. *GSEApy: a comprehensive package for performing
   gene set enrichment analysis in Python*. Bioinformatics 2023, 39(1):btac757.
   DOI: 10.1093/bioinformatics/btac757
2. Enrichr —— https://maayanlab.cloud/Enrichr/(基因集与富集后端)
