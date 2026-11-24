# TransEffect — 转录效应分析模块（WES somatic × RNA-seq）

对 DRclass 分类结果做"变异的转录效应"分析（《3.2.2.3节写作构思.md》第四节）：
变异等位基因表达、RNA/DNA VAF 比值、变异基因表达偏离（组内对照 Δ' + 随机
null）、CNV 表达剂量一致性、稳定性程序基因反向检验（P5），并输出 4 张图。
引擎为本模块自研脚本 `bin/transeffect.py`（analyze 仅标准库，figures 需 matplotlib）。

## 1. 模块输入输出与参数

### 1.1 输入（全部为具体文件路径，无目录扫描）

| 输入 | 说明 |
|---|---|
| `samples` (config) | 参与分析的样本 ID 列表（DRclass 已分类的样本） |
| `sample_group_dict` (config) | sample → 分组名；Δ' 组内对照按此分组，组内需 ≥2 样本 |
| `drclass_dir` | DRclass 输出根，消费 `<sample>/<sample>.dna_anchored.tsv` 与 `<sample>/anchors.tsv` |
| `counts_matrix` | featureCounts 合并矩阵（`{genome}_featureCounts.tsv`），列名含样本 BAM 路径 |
| `cnv_dir` | CNVkit 输出根，消费 `<sample>/<sample>.call.cns` |
| `gmt`（可选） | 基因集 GMT；缺省跳过 P5 反向检验及对应图 |
| `gene_sets`（可选） | GMT 内参与交集的集合名（\| 分隔）；缺省 12 个基因组稳定性通路 |
| `gtf`（可选） | gene_id→gene_name 映射来源（图内标注 / P5 符号映射） |

### 1.2 输出（{outdir}/ 下）

```
rna_support_sites.tsv          # 项1: RNA 支持位点(D1R1+D1R2) rna_AF/dna_AF 明细
vaf_ratios.tsv           # 项2: dna_AF/DP 下限过滤后的 RNA/DNA VAF 比值
variant_gene_delta.tsv   # 项6: 变异基因 Δ' + 随机 null(set/rep 列区分)
cnv_dosage_summary.tsv   # 项5: gain/loss/neutral 分类表达偏移 + 常染色体 cn≠2 计数
cnv_segments.tsv         # 项5 底层: CNVkit 段
program_overlap.tsv      # P5: 稳定性程序基因 ∩ 变异基因(需 gmt)
teffect_report.md        # 汇总报告(含判读口径)
fig_rna_support_vaf.png / fig_variant_gene_delta.png / fig_cnv_overview.png / fig_program_overlap.png
```

### 1.3 参数（config 位于 `Params.TransEffect.teffect_analyze`）

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `dna_af_min` | `--dna-af-min` | 0.1 | 项2 比值过滤的 dna_AF 下限 |
| `dna_dp_min` | `--dna-dp-min` | 10 | 项2 比值过滤的 dna_DP 下限 |
| `min_cpm` | `--min-cpm` | 1.0 | 表达基因过滤(CPM) |
| `min_expressed` | `--min-expressed` | 2 | 至少 N 样本达 CPM 阈值才算表达基因 |
| `null_reps` | `--null-reps` | 20 | 随机 null 重复次数 |
| `null_plot_reps` | `--null-plot-reps` | 3 | 入图/入表的 null 重复数 |
| `seed` | `--seed` | 42 | 随机种子 |
| `sex_chroms` | `--sex-chroms` | chrX,chrY,chrM,... | 性染色体名单(常染色体 cn≠2 判定排除) |
| `gene_sets` | `--gene-sets` | 12 个稳定性通路 | P5 基因集名(\| 分隔) |

`Procedure` 可覆盖可执行路径:`python3`。

直接调用示例（脱离流程）：

```bash
python3 bin/transeffect.py analyze \
  --samples KSR-1,KSR-2,TLSCS-1,TLSCS-2 \
  --sample-group KSR-1=KSR --sample-group KSR-2=KSR \
  --sample-group TLSCS-1=TLSCS --sample-group TLSCS-2=TLSCS \
  --drclass-dir .../DRclass --counts .../GRCm39_featureCounts.tsv \
  --cnv-dir .../DNA/variation/cnv --gmt .../GSI_mouse.gmt --gtf .../gencode.gtf.gz \
  --outdir .../TransEffect

python3 bin/transeffect.py figures --outdir .../TransEffect
```

经 run_module 调用（模块内直跑，自动注入 env 映射）：

```bash
python3 run_module.py -m TransEffect            # 默认目标 = 全部表 + 4 图
```

容器环境：`env.TransEffect` 指向 `/data/16T/Database/env/TransEffect/TransEffect.sif`
（由 `src/common/util/EnvUtil.py apptainer all -y modules/TransEffect/TransEffect.yaml
-o /data/16T/Database/env/TransEffect` 构建）。

## 2. 算法原理与统计口径

- **Δ'（组内对照）**：log2CPM(携带者) − mean(log2CPM(同组其余样本))。设计文档
  原口径（携带者 − 其余 3 样本）混入组间整体表达偏移（随机基因对照即有 ±0.3
  的假偏移），本模块实现为组内对照版；同组多携带的基因剔除（无法归因）。
- **随机 null**：等量随机表达基因、同剔除规则，`null_reps` 次重复，给 Δ' 分布
  的零假设参照。判读：变异基因 Δ' 与 null 一致 → 变异不伴随表达偏离。
- **CNV 剂量一致性**：表达基因按中点落入 call.cns 段分类(gain cn>2 / loss
  cn<2 / neutral)，逐样本算类内 Δ' 分布；另分开报告常染色体
  （^(chr)?数字+$）与未定位支架的 cn≠2 段数（雄性 chrX/chrY 半合子与支架
  伪段不是体细胞 CNV）。
- **P5 反向检验**：变异基因（anchors.tsv dna 锚落基因体）∩ 基因集成员；
  期望 = 程序基因数/总基因数 × 变异基因数，判富集与否。
- **图**：fig_rna_support_vaf（rna_AF vs dna_AF + 比值）、fig_variant_gene_delta
  （Δ' 箱线 vs null）、fig_cnv_overview（全基因组段拷贝数轨道）、
  fig_program_overlap（观测 vs 期望）。重绘只需 `figures` 子命令读 TSV。

## 3. 参考文献

1. 分析设计：`output/Totipotency/DR/3.2.2.3节写作构思.md` 第四节。
2. 上游分类：modules/DRclass（O1/O2/O3 输出为本模块的位点/标签来源）。
3. Danecek P et al. *Twelve years of SAMtools and BCFtools*. GigaScience
   10(2):giab008, 2021（上游 RNA 定量）。
