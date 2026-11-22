# tedominant — TE 主导 cluster 鉴定(单细胞,v1/v2 双方法)

从 scRNAseq 工作流 combine 阶段的 h5ad 对(scTE merged + scanpy auto)出发,
在每个组织内鉴定"TE 主导"的 Leiden cluster。两种方法(`method` 参数,
默认 **v1**):

- **v1(默认,旧口径)**:直接取 mean TE% 最高的 cluster。简单、与既往
  文献口径一致;但会被"基因程序塌缩 → TE% 被动抬高"的假阳性欺骗。
- **v2(反分母塌缩)**:先过绝对 TE 量资格门槛,再在合格簇中按 TE% 排名
  (§2.2)。判定更严格,图件口径同步切换为绝对量。

一次 group 级运行覆盖全部组织,产出鉴定表 + 三类图(burden 柱状图 /
UMAP / 衰老 marker dotplot);v1/v2 的图件差异见 §2.4。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{indir}/{tissue}/{tissue}_scTE_merged.h5ad` | scTE 定量合并对象:基因+TE 特征矩阵,须含整数 `counts` 层与 `var["gene_type"]`(TE 特征标记为 `"TE"`) |
| `{indir}/{tissue}/{tissue}_scTE_auto.h5ad` | scanpy 自动注释对象:须含 `obs["leiden"]`、`obs["cell_type"]`、`obsm["X_umap"]`,以及 `raw` 槽(全基因 log1p 矩阵,dotplot 用) |

两对象按 obs_names 取交集后配对使用。

### 1.2 输出

```
{outdir}/te_dominant_clusters.csv       # 全 cluster 指标表 + v2 判定(见 §2.3 列定义)
{outdir}/cluster_burden_{tissue}.png    # 每组织:per-cluster 绝对 TE 负荷柱状图
{outdir}/umap_te_{tissue}.png           # 每组织:UMAP 双面板(绝对 TE counts 上色 + 主导簇高亮)
{outdir}/markers_dotplot.png            # 衰老/谱系 marker dotplot(TE vs siblings vs rest)
{outdir}/te_dominant.sh                 # 回执脚本
logdir/te_dominant.log                  # 运行日志(含 [gate] 拒绝记录)
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `tissues` | `--tissues` | — | 组织名列表(输入子目录名) |
| `Params.tedominant.method` | `--method` | v1 | v1 = 纯 mean TE% 最高(旧口径);v2 = 绝对量门槛 + TE% 排名(§2.2) |
| `Params.tedominant.min_cells` | `--min-cells` | 100 | cluster 规模门槛 |
| `Params.tedominant.gate_frac` | `--gate-frac` | 0.5 | 绝对量门槛系数(两道门共用) |
| `Params.tedominant.skip_dotplot` | `--skip-dotplot` | false | 跳过 marker dotplot |

## 2. 指标定义、生物学含义与计算方式

### 2.1 细胞级与 cluster 级指标

所有指标基于 merged 对象的**整数 counts 层**(未归一化原始计数;
`layers["counts"]`,无则回退 `X`)。TE 特征由 `var["gene_type"] == "TE"`
界定(scTE 注入,不做名称匹配)。

对单个细胞 c,记:

- `TE_counts(c)` = 该细胞所有 TE 特征的 counts 之和(绝对 TE 转录量)
- `gene_counts(c)` = 非 TE 特征 counts 之和
- `total_counts(c)` = 两者之和
- `n_genes(c)` = counts > 0 的非 TE 特征数(检出基因数)

由此派生四个细胞级指标:

| 指标 | 计算 | 生物学含义 | 读法 |
|---|---|---|---|
| **te_abs**(绝对 TE 量) | `TE_counts(c)` | 该细胞内 TE 来源转录本的绝对丰度,直接反映 TE 转录去抑制的强度;与测序深度正相关,但**不受基因表达塌缩影响** | 高 = TE 转录真实升高 |
| **te_pct**(相对 TE 占比) | `TE_counts / total × 100` | TE 转录在本细胞总输出中的份额;**分母含基因表达**,基因程序沉默会被动推高它 | 高 = TE 主导,但需先过 te_abs 门槛排除假阳性 |
| **gene_abs**(绝对基因量) | `gene_counts(c)` | 蛋白编码基因程序的整体输出水平 | 低 = 基因程序塌缩(衰竭/退行信号) |
| **n_genes**(检出基因数) | `>0 的非 TE 特征数` | 转录复杂度/细胞健康度的经典 QC 代理 | 低 = 深度衰竭或低质量 |

cluster 级聚合:**中位数**(te_abs / gene_abs / n_genes,对离群细胞稳健)
与**均值**(te_pct,保持与旧文献的可比性),记入 `te_dominant_clusters.csv`。

### 2.2 v2 鉴定算法(反分母塌缩;`method: v2` 时启用)

**动机**:TE 主导的朴素定义 = mean TE% 最高的 cluster。但 te_pct 的分母
是总 counts——当一群细胞整体转录塌缩(基因输出暴跌、核糖体停滞、检出
基因数减半)时,即便 TE 转录量不变甚至下降,te_pct 也会被动冲到全组织
第一。这类"塌缩簇"是衰竭/低质量细胞,不是 TE 去抑制。

**实测案例**(本课题 Uterus):Leiden 2 的 TE% 38.5 全组织第一,但其
median te_abs = 918,仅为同谱系 siblings(4129)的 22%,基因检出数只有
siblings 的 30%——分母塌缩的典型指纹。旧算法会把它选为"TE 主导"。

**v2 两步判定**:

```
第一步 资格门槛(eligibility,三道门全过):
  gate_lineage: cluster median te_abs >= gate_frac × 同谱系 siblings median te_abs
                (siblings = 同 cell_type 的其他 Leiden cluster 细胞;n<30 时不设此门)
  gate_global : cluster median te_abs >= gate_frac × 组织内全部细胞 median te_abs
  gate_n      : cluster n_cells >= min_cells

第二步 排名:在 eligible 的 cluster 中取 mean te_pct 最高者为 TE 主导;
  全部不合格时回退纯 TE% 冠军并在日志 WARNING(不静默)。
```

**为什么门槛用 te_abs(中位数)而排名用 te_pct(均值)**:门槛回答
"TE 转录是否真的升高"(必须用绝对量,且中位数抗离群);排名回答
"在真升高的 cluster 里谁最 TE 主导"(占比衡量主导性)。两者缺一不可:
只用 te_abs 会把高测序深度的大簇误选,只用 te_pct 会中分母塌缩的招。

被拒的旧冠军会在日志打印 `[gate] ... TE abs X < 0.5×sib Y
(denominator-collapse guard)`,柱状图上以箭头标注——审计痕迹完整。

### 2.3 `te_dominant_clusters.csv` 列定义

| 列 | 含义 |
|---|---|
| `tissue` / `cluster` / `cell_type` / `n_cells` | 组织、Leiden ID、主要注释(众数)、细胞数 |
| `te_pct_mean` | cluster 内 te_pct 均值(v2 排名指标) |
| `te_abs_median` | cluster 内 te_abs 中位数(v2 门槛指标) |
| `gene_abs_median` / `n_gene_median` | 基因输出与检出基因数中位数(塌缩诊断) |
| `sib_te_abs_median` | 同谱系 siblings 的 te_abs 中位数(gate_lineage 参照;独占型 cell_type 为 inf) |
| `tissue_te_abs_median` | 组织全局 te_abs 中位数(gate_global 参照) |
| `gate_lineage` / `gate_global` / `gate_n` | 三道门各自是否通过 |
| `eligible` | 三门交集(v2 资格) |
| `dominant` | 是否被选为该组织的 TE 主导 cluster |

### 2.4 图件与指标口径的一致性

图件口径随 `method` 切换:v1 图回到相对 TE% 口径(排序、虚线、UMAP 上色),
v2 图用绝对量口径并画门槛与拒绝标注:

| 图 | v1 口径 | v2 口径 | 设计意图 |
|---|---|---|
| `cluster_burden_*.png` | y 轴 = te_pct_mean,虚线 = 全组织均值 TE%,红条 = TE% 冠军 | y 轴 = te_abs_median,虚线 = gate_frac × tissue median,箭头标注被拒旧冠军 | v1 复现旧图;v2 把 §2.2 判定过程画出来,被拒簇的低矮条 + 拒绝理由一目了然 |
| `umap_te_*.png` 左面板 | 颜色 = te_pct(0–60) | 颜色 = te_abs(98 分位截尾) | v1 与旧口径一致;v2 与门槛同源,相对 TE% 上色会"点亮"塌缩细胞造成误导 |
| `umap_te_*.png` 右面板 | 灰底 + 蓝色高亮主导簇(两种方法相同) | 同左 | 空间定位 |
| `markers_dotplot.png` | size = 检出率(% expressing),color = mean log1p(auto.raw.X 全基因矩阵,非 HVG 子集);不随 method 变 | 同左 | 衰老级联定位;raw 槽是 log1p 归一化前的全基因快照,低表达 SASP 基因不在 HVG 里,必须用 raw 才不静默丢基因 |

dotplot 的 marker 分组:SASP/炎症、CDKi 与检查点、DNA 损伤与衰老
(注意 LMNB1 是反向 marker,FOXO/SIRT 是抗应激方向)、增殖(衰老的反面)、
应激/凋亡、谱系身份对照(granulosa / glandular epithelial)。

### 2.5 统计语义与限制

- 每 `gate_frac` / `min_cells` 是启发式默认(0.5 / 100),阈值变化±0.1
  内判定稳健即可采信;报告建议附 `te_dominant_clusters.csv` 全表供复核。
- 主导簇 vs siblings 的表达对比是**描述性**的;组间(n≥2 生物学重复)
  显著性需另走 pseudobulk DEG。
- te_pct 跨组织不可直接比绝对色值(分母结构不同);跨组织只比
  "te_abs 中位数 + 判定结论"的模式。
- markers_dotplot 的 color 用 `auto.raw.X`:它是 log1p(target-sum) 归一化
  值,未做 batch 校正或 scale,跨组织绝对色值比较同样需谨慎。

## 3. 引擎与算法实现

- 引擎:`bin/te_dominant.py`(纯 Python,scanpy 栈);判定核心
  `dominant_v2()` 与制图在容器内以真实数据验证:ovaries Leiden 18
  (Granulosa, TE% 65.9, te_abs 3370 ≥ 门槛)与 Uterus Leiden 1
  (Glandular, TE% 37.6, te_abs 4142 = siblings 3.5×)当选;
  Uterus Leiden 2 被门槛拒绝(918 vs 4129)。
- 环境:`tedominant.yaml`(scanpy 栈);SIF 复用规则见 modules.md。
- 接线:group 级单规则 `te_dominant`(无样本通配符,产物含 done 标记)
  + 收尾 `tedominant_result`。

## 4. 参考文献

1. scTE:Yan Z et al. *Nat Commun* 2021(single-cell TE 定量;merged
   对象的 gene_type 注入)—— https://doi.org/10.1038/s41467-021-25733-2
2. Ghandi M et al. *Nucleic Acids Res* 2016(te BuccSC 框架中 TE 占比
   定义的先例与分母效应讨论)
3. v2 门槛判定(§2.2)为本仓库实现,无外部算法依赖;
   案例数据:luanchao_mulatta scRNAseq(ovaries/Uterus)。
