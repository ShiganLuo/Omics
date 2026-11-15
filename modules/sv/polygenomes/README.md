# sv/polygenomes — 多基因组结构变异(SV)比较分析子模块

对成对(control vs experiment)结构变异 VCF 做实验特异 SV 提取、VEP 注释、
GO/KEGG 富集、Circos 图、基因模型图、跨比较差异分析与 OncoPrint。C 形态
子模块,共享父目录 `../sv.yaml` 环境与 `modules/sv/bin/` 脚本
(exp_specific.py、run_enrichment.py、run_circos.py、utils/gene_model.py、
pbsv_sv_diff_analysis.py、run_OncoPrint.py)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `comparisons` | 比较定义列表,每项 `{name, control_vcf, experiment_vcf}`(VCF 路径直接给在 config;wildcard `comparison` 取 `name`) |
| `genome.<g>.fasta` / `.gtf` / `.fai` | 每基因组参考 FASTA、注释 GTF、FASTA 索引(平铺键,直接 `config["genome"][<g>][...]` 读取) |
| `Params.oncoprint.cosmic_file` | COSMIC 癌基因列表(OncoPrint 输入) |
| `Params.oncoprint.deseq2_files` | 条件 → DESeq2 结果路径字典(OncoPrint 叠加表达) |

### 1.2 输出

```
{outdir}/{genome}/{comparison}/{comparison}_merged.vcf       # 合并 VCF
{outdir}/{genome}/{comparison}/{comparison}_only.vcf         # 实验特异 SV
{outdir}/{genome}/{comparison}/{comparison}_annotated.tab    # VEP 注释表
{outdir}/{genome}/{comparison}/enrichment/                    # GO/KEGG 富集(directory)
{outdir}/{genome}/{comparison}/circos/                        # Circos 图(directory)
{outdir}/{genome}/{comparison}/gene_model/                    # 基因模型图(directory)
{outdir}/{genome}/sv_diff_analysis/                           # 跨比较差异分析(directory)
{outdir}/{genome}/sv_oncoprint/oncoprint_matrix.csv            # OncoPrint 矩阵
{outdir}/{genome}/sv_oncoprint/                                # OncoPrint 图(directory)
```

日志:`{logdir}/{genome}/{comparison}/sv_*.log`、`{logdir}/{genome}/all/sv_*.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `comparisons` | `-c`/`-e`(exp_specific) | `[]` | 比较定义(见输入) |
| `Params.formats` | `-f`(可重复) | `["png"]` | 出图格式列表(matplotlib 全集;旧键 `Params.image_formats` 兼容读取,以 `formats` 为准) |
| `Params.genes` | — | `[]` | gene_model 目标基因列表 |
| `Params.exp_specific.dist` | `-d` | `500` | SV 邻近判定距离(bp) |
| `Params.exp_specific.vep_cache` | `--vep_cache` | `~/.vep` | VEP 缓存目录 |
| `Params.exp_specific.species` | `--species` | `mus_musculus` | VEP 物种 |
| `Params.exp_specific.assembly` | `--assembly` | `GRCm39` | VEP 组装版本 |
| 固定(exp_specific) | `--annotate_format tab` | — | 注释输出为 TSV |
| `Params.circos.genome` | `--genome` | `mm39` | Circos 基因组配置 |
| `Params.circos.cytoband` | `--cytoband` | `""` | 非空才传 |
| `Params.circos.ins_bin_size` | `--ins_bin_size` | `100000` | 插入分箱(bp) |
| `Params.gene_model.threads` | `-j` | `1` | gene_model 并行度 |
| `Params.diff_analysis.large_sv_threshold` | `-s` | `10000` | 大 SV 阈值(bp) |
| `Params.oncoprint.cosmic_file` | `--cosmic_file` | `""` | 见输入 |
| `Params.oncoprint.deseq2_files` | 条件参数 | `{}` | 见输入 |
| `ROOT_DIR` | — | `.` | 拼 `modules/sv/bin/` 脚本 |

环境:`conda: "../sv.yaml"`。

直接调用示例(等价引擎命令):

```bash
python modules/sv/bin/exp_specific.py \
  -c control.vcf -e experiment.vcf -o out/cmp1 -d 500 \
  --vep_cache ~/.vep --species mus_musculus --assembly GRCm39 --annotate_format tab
python modules/sv/bin/run_enrichment.py \
  -a out/cmp1_annotated.tab -o out/enrichment -g annotation.gtf -fi genome.fa.fai
python modules/sv/bin/run_circos.py \
  --vcf out/cmp1_only.vcf --fasta genome.fa --outdir out/circos \
  --genome mm39 --ins_bin_size 100000
python modules/sv/bin/pbsv_sv_diff_analysis.py -o out/sv_diff_analysis -s 10000
python modules/sv/bin/run_OncoPrint.py --cosmic_file cosmic.txt -o out/sv_oncoprint
```

## 2. 算法原理与生物学意义

围绕「实验特异 SV」的下游解读链:`exp_specific.py` 合并 control/experiment
VCF,提取只在实验组出现的 SV(距离阈值 `dist` 内判同一位点),用 VEP 注释
断点邻近基因与功能影响;`run_enrichment.py` 对受影响基因做 GO/KEGG 富集并
找突变热点;`run_circos.py` 全基因组环视图展示 SV 类型与密度;`gene_model.py`
在目标基因(`Params.genes`)结构图上叠加 SV 位置;`pbsv_sv_diff_analysis.py`
跨比较组做差异 SV(按 `large_sv_threshold` 区分大小 SV);`run_OncoPrint.py`
把 SV × 基因事件矩阵与 COSMIC 癌基因、DESeq2 表达联合可视化。结果怎么读:
`_only.vcf` 是核心(实验特异 SV);`_annotated.tab` 看断点基因;OncoPrint 行
为基因、列为比较,色块 = 事件类型,适合找复发驱动 SV。

## 3. 参考文献

1. VEP(Ensembl Variant Effect Predictor)—— https://www.ensembl.org/info/docs/tools/vep/
   (McLaren W et al., Genome Biology 17:122, 2016. DOI: 10.1186/s13059-016-0974-4)
2. Circos — http://circos.ca/(Krzywinski M et al., Genome Research
   19(9):1639–1645, 2009. DOI: 10.1101/gr.092759.109)
