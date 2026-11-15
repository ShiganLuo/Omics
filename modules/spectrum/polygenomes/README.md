# spectrum/polygenomes — 多基因组体细胞突变谱与签名分析子模块

polygenomes(每基因组独立目录)布局下的突变谱三规则:SBS96 谱计数 →
COSMIC 已知签名 refit(NNLS)→ de novo NMF 签名提取,每基因组一套
产出。共享父目录 `../spectrum.yaml` 环境与 `modules/spectrum/bin/`
脚本(`spectrum.py`、`signature_refit.py`、`signature_nmf.py`);
原理、影响因素与签名矩阵库见 [spectrum/README.md](../README.md) 与
同目录《突变签名分析_原理与影响因素.md》。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `sample_somatic_vcf_dict` (config,顶层) | 样本 → 体细胞 VCF 完整路径字典(形态 B,node.py 从 tumor/normal 配对派生;仅读 VCF 前 6 列) |
| `sample_group_dict` (config,顶层) | 样本 → 分组标签(谱计数按组聚合) |
| `genome.<g>.fasta` (config) | 每基因组参考 FASTA(三核苷酸上下文提取);缺失 build 期 `ValueError` |
| `genome.<g>.signature_matrix` (config) | COSMIC 签名矩阵 TSV(SBS 96 通道,按 build 选,如 `/data/16T/Database/cosmic_signatures/mm39/COSMIC_v3.6_SBS_mm39.txt`);缺失 build 期 `ValueError` |
| wildcards `genome` | 基因组 ID(输出按基因组分目录) |

### 1.2 输出(每基因组,`{outdir}/{genome}/`)

```
# somatic_spectrum
somatic_spectrum_per_sample.csv          # 96 通道 × 样本计数(refit/NMF 输入契约)
somatic_spectrum_sample_stats.csv        # 每样本 SNV 计数与丢弃账目
somatic_spectrum.csv                     # 96 通道 × 分组计数
somatic_spectrum_stacked_bar.<fmt> / _heatmap.<fmt>   # <fmt> 取自 plot.formats(默认 png)
experiment_vcf_map.tsv / experiment_group_map.tsv   # 运行时落盘的映射表

# somatic_signature_refit
signature_refit_exposures.csv / _exposure_fractions.csv / _family_exposures.csv
signature_refit_summary.csv
signature_refit_exposure_stacked_bar.<fmt> / _exposure_heatmap.<fmt>
signature_refit_signature_cosine.<fmt>

# somatic_signature_nmf
signature_nmf_signatures.csv / signature_nmf_exposures.csv
signature_nmf_rank_diagnostics.csv / .<fmt>
signature_nmf_cosmic_similarity.csv / _cosmic_best_match.csv / _cosmic_match.<fmt>
```

(refit/NMF 图类产出的 `<fmt>` 取自 `signature.formats`,默认 `["png"]`,
可多格式同时输出。)

回执:`somatic_spectrum.sh`、`somatic_signature_refit.sh`、
`somatic_signature_nmf.sh`;日志 `{logdir}/{genome}/spectrum/*.log`。
`spectrum_result` 聚合全部基因组关键产出。

### 1.3 参数

config 位于顶层 `signature` 段(全基因组共享):

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `plot.formats` | `spectrum.py -f/--format`(可重复) | `["png"]` | 谱图(`_stacked_bar`/`_heatmap`)输出格式;合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集) |
| `signature.formats` | refit `--format` / nmf `-f/--format`(可重复) | `["png"]` | 签名图(refit 3 图 + NMF 3 图)输出格式,同上 |
| `signature.min_mutations` | `--min-mutations`(refit) | 200 | 低于此 SNV 数样本标 unreliable |
| `signature.min_exposure_frac` | `--min-exposure-frac`(refit) | 0.02 | refit 曝光剪枝阈值(逐轮剔除后重拟合) |
| `signature.family_map` | `--family-map` | null | 签名家族映射 TSV;null 用内置表 |
| `signature.max_rank` | `--max-rank`(nmf) | 8 | NMF 秩扫描上限 |
| `signature.restarts` | `--restarts`(nmf) | 20 | 每秩随机重启次数 |
| `signature.seed` | `--seed`(nmf) | 42 | 随机种子 |
| `signature.rank` | `--rank`(nmf) | null | 强制秩;null 自动选(KL 拐点 + cophenetic 下限) |
| `signature.mode` / `enabled` | — | `both` / true | .smk 未消费(`mode`/`enabled` 的注册门控在 node.py 侧,见父 README) |
| (threads) | — | 1 | `rule_res` 可覆盖;回执内 export OMP/OPENBLAS/MKL/NUMBA 线程数 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python modules/spectrum/bin/spectrum.py \
    --vcf-map experiment_vcf_map.tsv --group-map experiment_group_map.tsv \
    --ref-fasta genome.fa --output-prefix out/genome1/somatic_spectrum -f png -f pdf

python modules/spectrum/bin/signature_refit.py \
    --counts out/genome1/somatic_spectrum_per_sample.csv \
    --signatures COSMIC_v3.6_SBS.txt \
    --output-prefix out/genome1/signature_refit \
    --min-mutations 200 --min-exposure-frac 0.02 --format png --format pdf

python modules/spectrum/bin/signature_nmf.py \
    --counts out/genome1/somatic_spectrum_per_sample.csv \
    --signatures COSMIC_v3.6_SBS.txt \
    --output-prefix out/genome1/signature_nmf \
    --max-rank 8 --restarts 20 --seed 42 -f png -f pdf
```

## 2. 算法原理与生物学意义

VCF → SBS96 通道计数(突变类型 × 三核苷酸上下文)→ 双路签名分析:

- **refit(NNLS)**:把每样本谱非负最小二乘拟合到 COSMIC 已知签名,
  `min_exposure_frac` 剪枝后重拟合,输出签名曝光(突变数/比例/家族
  级)与重构余弦;签名两两余弦高的对存在简并,应读家族级
  (Clock-like/APOBEC/MMR/POLE_POLD1/HRD/UV 等);
- **NMF(de novo)**:对计数矩阵做非负矩阵分解提取新签名,`max_rank`
  扫描 + restarts 稳定性评估自动选秩,再与 COSMIC 求余弦注释最佳
  匹配(`cosmic_match.<fmt>` 可目检证据);
- **polygenomes 差异**:同一 VCF 集对每个基因组各跑一遍——上下文与
  签名矩阵随参考 build 变化(如 mm39 小鼠矩阵),跨物种 COSMIC 拟合
  仅作参考,不构成同源赋值(详见父 README 跨物种提示与原理文档 §5.6)。

生物学意义:突变签名反映致突变过程(老化 Clock-like、APOBEC、错配
修复缺陷、POLE/POLD1 校对缺陷、HRD、UV),是病因归因与治疗敏感性
推断(HRD → 铂类/PARP)的关键读出;`min_mutations` 以下样本的曝光
估计采样噪声大,须谨慎。

## 3. 参考文献

1. COSMIC Mutational Signatures — https://cancer.sanger.ac.uk/signatures/
2. Alexandrov LB et al. *The repertoire of mutational signatures in human
   cancer*. Nature 578:94-101, 2020. DOI: 10.1038/s41586-020-1943-3.
3. 方法细节与影响因素见本仓库 `modules/spectrum/突变签名分析_原理与影响因素.md`;
   签名矩阵由 SigProfilerAssignment 分发包解包(见父 README)。
