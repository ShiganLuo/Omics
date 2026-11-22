# teprofile — TE 表达三层刻画(本源特征 / 异常亚群 / 分组重塑)

从 scRNAseq 工作流 combine 阶段的 h5ad 对(scTE merged + scanpy auto)
出发,把"TE 与细胞群体"的问题拆成三层递进分析,替代"选出一个 TE 主导
cluster"的旧叙事:TE 表达是连续谱 + 谱系依赖的,硬切二元标签会被分母
效应与塌缩细胞污染(见 tedominant 模块 README §2.2 的 Uterus 案例)。

- **Layer 1 本源特征**:各细胞类型的 TE 表达长什么样(基线)
- **Layer 2 异常亚群**:哪些 cluster 显著偏离自己的谱系背景(不是冠军选举)
- **Layer 3 分组重塑**:组间(如 young vs aged)基线平移与亚群比例变化

## 1. 模块输入输出与参数

### 1.1 输入(与 tedominant 相同的 h5ad 对)

| 输入 | 说明 |
|---|---|
| `{indir}/{tissue}/{tissue}_scTE_merged.h5ad` | 基因+TE 特征矩阵,须含整数 `counts` 层与 `var["gene_type"]`/`var["family_id"]` |
| `{indir}/{tissue}/{tissue}_scTE_auto.h5ad` | 注释对象(`obs["leiden"]`/`cell_type`/`sample_id`) |

### 1.2 输出

```
te_profile_baseline.csv         layer 1:每 (tissue, cell_type) 的特征摘要
te_profile_clusters.csv         layer 2:每 (tissue, cluster) 的偏离画像 + 标签
te_profile_remodeling.csv       layer 3:每 (tissue, cell_type) 的组间平移(groups 配置时)
baseline_te_fraction_<t>.png    每组织:TE 占比 + HK 归一化 TE 产量 violin(按细胞类型)
baseline_family_comp_<t>.png    每组织:TE 家族组成堆叠条(按细胞类型)
cluster_deviation_<t>.png       每组织:三指标偏离面板(全部 cluster,不选冠军)
remodeling_shift_<t>.png        每组织:layer-3 组间平移摘要(groups 配置时)
te_profile.sh / te_profile.done 回执与标记
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `tissues` | `--tissues` | — | 组织名列表 |
| `Params.teprofile.fe_cut` | `--fe-cut` | 1.5 | layer-2 TE 富集倍数门槛 |
| `Params.teprofile.comp_cut` | `--comp-cut` | 0.25 | layer-2 家族组成 Aitchison 门槛 |
| `Params.teprofile.groups` | `--groups` | "" | layer-3 分组 `'组A:s1,s2;组B:s3'`;空则跳过 |

## 2. 三层分析的指标、生物学含义与计算方式

### 2.0 核心归一化:HK 校正的 TE 产量(解决 v1/v2 之争的指标)

```
te_norm = TE_counts / max(sum(HK gene counts), 1) × 100
```

- **为什么不用 te_pct**(TE/总 counts):分母含全部基因,基因程序塌缩的
  细胞分母暴跌,占比被被动推高(分母塌缩,Uterus Leiden 2 实测:占比全
  组织第一而绝对 TE 量仅为同谱系 22%)
- **为什么不用原始 te_abs**:受测序深度混杂
- **HK 归一化同时抗两者**:管家基因(ACTB/GAPDH/RPL 家族等 16 个候选与
  数据基因符号取交集,≥5 个才启用,否则回退全部基因并标注)深度稳定、
  且不是转录组的主体——塌缩细胞 HK 降幅远小于总 counts 降幅。te_norm
  的语义是"每 100 条管家转录本对应的 TE 转录本数",即**深度校正后的
  TE 转录产量**
- 生物学含义:te_norm 高 = 该细胞在单位转录活性下生产更多 TE RNA,
  即 TE 转录去抑制的强度,可比地跨细胞、跨 cluster、跨样本

### 2.1 Layer 1 — 细胞类型本源特征

按 `cell_type`(注释层,非 Leiden)聚合,每类输出:

| 指标 | 计算 | 生物学含义 |
|---|---|---|
| `te_pct_median/q25/q75` | 占比分位数 | 该类型转录组的 TE 份额(构成描述,保留与文献可比性) |
| `te_norm_median/q25/q75` | HK 归一化产量分位数 | 该类型真实的 TE 转录强度(机制指标) |
| `n_gene_median` | 检出基因数中位 | 转录复杂度/健康度(QC 语境) |
| `top_family` | 家族份额均值最大者 | 该类型 TE 谱系偏好(如内皮 ERVK 高) |

配套 violin 图(占比 + log 轴 te_norm)与家族组成堆叠条。**分布形态
本身是结论**:同类型内 TE 特征的宽度反映状态异质性,类型间的差异即
"TE 表达受细胞身份调控"的证据。

### 2.2 Layer 2 — 亚群离群检测(与同谱系背景比)

每个 Leiden cluster 与**同 cell_type 的其余细胞**(n<30 时回退全组织)
比三个正交指标:

| 指标 | 计算 | 含义 |
|---|---|---|
| `fe_vs_lineage` | median(te_norm_cluster) / median(te_norm_谱系背景) | TE 产量倍数富集(≥ fe_cut 为异常) |
| `family_aitchison` | 家族份额(CLR 变换)均值的 Aitchison 距离 | TE 家族组成偏移(哪些家族异常活跃;≥ comp_cut 为异常) |
| `n_gene_ratio` | median(n_gene_cluster) / median(n_gene_背景) | 基因程序塌缩检测(<0.7 判 collapsed) |

三指标组合给标签(优先级从上到下):

- `collapsed`:基因程序塌缩(n_gene_ratio<0.7)——**先于 TE 判定**,
  塌缩簇的 te_pct 不可信,即便 fe 高也标塌缩而非 TE 异常
  (Uterus Leiden 2 在此框架下正确归位:fe 2.9 但 n_gene_ratio 0.26)
- `te_outlier`:fe ≥ cut 且组成偏移 ≥ cut(量与质同时异常 = 真去抑制)
- `te_enriched`:仅 fe ≥ cut(量高但家族构成正常 = 均匀升高)
- `typical`:其余

**没有冠军选举**:全部 cluster 都进表与图(cluster_deviation 三面板按
统一行序展示三指标),允许"该组织无 te_outlier"作为结论。

### 2.3 Layer 3 — 分组重塑(groups 配置时)

按 `sample_id` 分组(如 young/aged),每 cell_type 比较两组的
te_pct / te_norm 中位数平移(`Δ`),并列出两组细胞数。回答"衰老是
整体平移还是特定类型变化";亚群比例变化的完整检验应走 scCODA
(diff/abundance),此层只做描述性平移。样本数 <2 组或组内 <20 细胞
的类型自动跳过。

### 2.4 统计语义与限制

- HK 基因集是人类 curate 列表与猕猴符号的交集,依赖符号保守;若交集
  <5 个自动回退"全部基因"归一化并在日志/README 标注(hk_mode 列)
- `fe_cut`/`comp_cut` 是描述性门槛,报告应附全表(te_profile_clusters.csv)
  供审稿人按其他阈值复核
- Layer 3 是两组描述性比较,不做显著性;n=1/组的样本无法区分生物与批次
- Aitchison 距离对稀有家族敏感;小 cluster(n<100)的组成偏移可能是
  噪声,解释时结合 n_cells 列

## 3. 引擎与实现

- 引擎:`bin/te_profile.py`(scanpy 栈,纯 Python)。真实数据验证:
  ovaries 17 类 / Uterus 12 类基线;layer 2 ovaries 检出 1 个
  te_outlier(Ovarian surface epithelial Leiden 12:fe 1.67、Aitchison
  4.26、n_gene_ratio 1.03——健康背景上的真 TE 异常),Uterus 0 个
  te_outlier、6 个 collapsed(含 Leiden 2:fe 2.9 但 n_gene_ratio 0.26,
  正确识别为塌缩而非 TE 主导)
- HK 归一化模式(housekeeping/all_genes)写进日志与图轴标签
- 环境:teprofile.yaml(scanpy 栈,复用 tedominant 同级依赖)

## 4. 参考文献

1. scTE:Yan Z et al. *Nat Commun* 2021 — TE 特征的 gene_type/family_id
   注入。https://doi.org/10.1038/s41467-021-25733-2
2. Aitchison J. *J R Stat Soc B* 1982 — 成分数据 CLR/Aitchison 距离。
3. HK 归一化与三层框架为本仓库实现;案例数据 luanchao_mulatta
   scRNAseq(ovaries/Uterus)。
