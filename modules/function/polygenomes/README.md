# function/polygenomes — 多基因组差异表达功能富集子模块(GO/KEGG/GSEA/GSVA)

对 DESeq2 差异结果做功能富集:GO/KEGG 上下调富集、GSEA(基因集富集分析)、
GSVA(单样本变异打分)。C 形态子模块,共享父目录 `../function.yaml` 环境与
父目录 `modules/function/bin/` 脚本(go-kegg.r、gsea.r、normalization.py、
gsva.py)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(DESeq2 输出目录) | 形态 A 拼接:GO/KEGG 用 `{indir}/{genome}/{contrast}/TEcount/{contrast}.TEcount_Gene.name.tsv`;GSEA 用同目录 `{contrast}.TEcount_Gene.tsv` |
| `counts_dir` | featureCounts 合并矩阵 `{counts_dir}/{genome}/{genome}_featureCounts.tsv`(GSVA 输入) |
| `genome.references.<g>.geneIDAnno` | 基因 ID 注释 CSV(GSEA `-a`);缺失报错。注意 .smk 读 `genome.references`(复数),json 示例平铺 `genome.<g>.geneIDAnno`,以 .smk 为准 |
| `genome.references.<g>.gtf` | 注释 GTF(GSVA 基因 ID→名转换) |
| `Params.function.<genome>.gsea.gmt` | GSEA 基因集 GMT,必需 |
| `Params.function.<genome>.gsva.gmt` | GSVA 基因集 GMT,必需 |
| `comparisons` | 对比组定义(顶层读取,规则目标由 subworkflow 按 `{contrast}` 注册) |

### 1.2 输出

```
# 图按 formats 每种格式一份(go_kegg 默认 png,gsea 默认 jpeg,gsva 默认 png)
{outdir}/{genome}/{contrast}/go_back_to_back.<fmt>        # GO 上下调背靠背图
{outdir}/{genome}/{contrast}/kegg_back_to_back.<fmt>      # KEGG 上下调背靠背图
{outdir}/{genome}/{contrast}/go_up.csv / go_down.csv    # GO 富集表(上/下调)
{outdir}/{genome}/{contrast}/kegg_up.csv / kegg_down.csv
{outdir}/{genome}/{contrast}/up_genes.txt / down_genes.txt
{outdir}/{genome}/{contrast}/GSEA/TEcount_Gene_GSEA.<fmt>
{outdir}/{genome}/{contrast}/GSEA/TEcount_Gene_GSEA.csv
{outdir}/{genome}/gsva/gsva_scores.tsv                  # 基因集×样本 富集分数
{outdir}/{genome}/gsva/gsva_heatmap.<fmt>               # GSVA 聚类热图
```

日志:`{logdir}/function/{genome}/{contrast}.go_kegg.log`、`.gsea.log`、
`gsva.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `counts_dir` / `outdir` / `logdir` | — | `output/diff_expression` / `output/counts` / `output/function` / `logs` | 目录 |
| `Params.function.<genome>.go_kegg.species` | `-s` | `mouse` | 物种(orgdb 选择) |
| `Params.function.<genome>.go_kegg.lfc_cut` | `--lfc-cut` | `1` | 上下调基因 log2FC 阈值 |
| `Params.function.<genome>.go_kegg.p_cut` | `--p-cut` | `0.05` | 显著性阈值 |
| `Params.function.<genome>.go_kegg.top` | `--top` | `10` | 展示条目数 |
| `Params.function.<genome>.go_kegg.formats` | `--format`(可重复) | `[png]` | 输出图片格式(多种同时输出);合法值 png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集;R 侧 ggsave 设备集)。go-kegg.r 为 R 脚本,仅受理 ggsave 设备集,pgf/raw/rgba 传入即报错 |
| `Params.function.<genome>.gsea.gmt` | `-g` | 必需 | GSEA 基因集 |
| `Params.function.<genome>.gsea.formats` | `--format`(可重复) | `[jpeg]` | GSEA 图输出格式(gsea.r 为 R 脚本,ggsave 设备集) |
| `Params.function.<genome>.gsva.gmt` | `-g` | 必需 | GSVA 基因集 |
| `Params.function.<genome>.gsva.formats` | `-f`/`--format`(可重复) | `[png]` | GSVA 热图输出格式(matplotlib 全集) |
| `ROOT_DIR` | — | `.` | 拼父目录 `modules/function/bin/` 脚本 |

注意:各 genome 段的 `formats` 取值必须一致(规则 output 在 parse 期展开,
`{genome}` 通配符维度上不能各自不同);只在部分 genome 段设置时其余沿用该
取值,多个 genome 设置了不同取值会在 parse 期报错。

固定调用:go-kegg.r `-i <deseq2 tsv> -o <outdir> -s <species> --lfc-cut --p-cut
--top`;gsea.r `-m Gene -g <gmt> -i <tsv> -o <outdir> -a <geneIDAnno> -t <title>`;
GSVA 两步:`normalization.py -i <counts> -o <tpm> --method tpm --gtf <gtf>
--convert-gene-name` → `gsva.py -i <tpm> -g <gmt> -o <outdir> --title <title>`。
json 另声明 `Params.function.<genome>.enabled/species`(`func_params` 顶层读取,
本 .smk 各规则直接按 genome 取子键)。环境:`conda: "../function.yaml"`。

直接调用示例(等价引擎命令):

```bash
Rscript modules/function/bin/go-kegg.r \
  -i out/GRCm39/aged_vs_youth/TEcount/aged_vs_youth.TEcount_Gene.name.tsv \
  -o out/function/GRCm39/aged_vs_youth -s mouse --lfc-cut 1 --p-cut 0.05 --top 10
Rscript modules/function/bin/gsea.r -m Gene -g sets.gmt \
  -i out/GRCm39/aged_vs_youth/TEcount/aged_vs_youth.TEcount_Gene.tsv \
  -o out/function/GRCm39/aged_vs_youth -a geneIDAnnotation.csv -t "GSEA"
python modules/function/bin/normalization.py -i counts.tsv -o tpm.tsv \
  --method tpm --gtf annotation.gtf --convert-gene-name
python modules/function/bin/gsva.py -i tpm.tsv -g sets.gmt \
  -o out/function/GRCm39/gsva --title "GSVA (mouse) - GRCm39"
```

## 2. 算法原理与生物学意义

- **GO/KEGG**:按 `lfc_cut`/`p_cut` 切上下调基因,分别做过表征富集
  (超几何/Fisher),背靠背图并排展示上调/下调条目,回答"差异基因落在
  哪些通路/功能"。
- **GSEA**:不切基因集,按全基因排序统计量检验基因集是否富集于排序两端,
  输出 NES 与 FDR,灵敏度高于硬阈值法,适合微弱一致的通路信号。
- **GSVA**:先 TPM 归一化(`normalization.py`),再对每个样本独立计算基因集
  富集分数,得到 geneset × sample 矩阵与聚类热图,适合多样本/多分组比较
  通路活性(如衰老通路随年龄的分数变化)。

结果怎么读:go/kegg `*_up/down.csv` 为富集表(GeneRatio、pvalue、qvalue、
基因数),`top` 控制图上条目;GSEA `*.csv` 关注 NES(方向)与 FDR;
`gsva_scores.tsv` 行为基因集、列为样本,数值越高该通路越活跃。

## 3. 参考文献

1. clusterProfiler(Yu G et al., OMICS 2012.
   DOI: 10.1089/omi.2011.0118)—— GO/KEGG/GSEA 常用实现
2. GSVA — Hänzelmann S et al. *GSVA: gene set variation analysis for
   microarray and RNA-seq data*. BMC Bioinformatics 14:7, 2013.
   DOI: 10.1186/1471-2105-14-7
