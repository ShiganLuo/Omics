# scRNAseq 自动注释：设计与踩坑

## 背景

单细胞注释的传统做法是人工看 marker genes、查文献、画 UMAP，一个 dataset 来回折腾几天。我们做了一套全自动流程——聚类、LLM 注释、质量检查、迭代优化，跑一次出结果。

核心思路：用 Step 0 查组织已知细胞类型作为 reference，Step 2 让 LLM 看着 top DEGs 和 specificity table 注释每个 cluster，Step 3 做 QC，不干净的过滤掉重跑。听起来简单，实际踩了不少坑。

## Wilcoxon DEG 的本质问题

`sc.tl.rank_genes_groups(method="wilcoxon")` 是 scanpy 的默认 DEG 方法。它对每个基因做 Wilcoxon rank-sum test，比较 cluster 内 vs 其他所有 cluster 的表达分布，按统计显著性排名。

问题出在排名逻辑上。以 Oocyte cluster 为例：

| Gene | Oocyte mean | Oocyte % | 其他 max mean | Wilcoxon 排名 |
|------|-------------|----------|---------------|--------------|
| MKI67 | ~3.0 | ~50% | ~0.01 | top 1 |
| TOP2A | ~2.5 | ~40% | ~0.005 | top 2 |
| NLRP5 | 0.082 | 4.3% | 0.010 | 不在 top 50 |
| DDX4 | 0.022 | 2.1% | 0.008 | 不在 top 50 |

MKI67 在 oocyte 里 50% 细胞高表达，其他 cluster 几乎不表达——effect size 极大，Wilcoxon score 碾压一切。NLRP5 只在 4.3% 细胞里表达，mean=0.082，统计上不够显著。

但生物学上，MKI67 是 oocyte 减数分裂的伴随信号，NLRP5 才是定义 oocyte 身份的基因。Wilcoxon 分不清这两者。

引入 cross-cluster specificity：

```
specificity = this_cluster_mean / max(other_clusters_mean)
```

NLRP5: 0.082 / 0.010 = 8.3x。DDX4: 0.022 / 0.008 = 2.9x。虽然绝对表达量低，但相对其他 cluster 有明显富集。

关键不是加了 specificity 这个指标，而是怎么把它喂给 LLM。最初 specificity table 放在 prompt 末尾，标记为"tie-breaker"。但 LLM 看到 50 个 DEGs 和 5 个 specificity markers 时，不觉得这是 tie——它顺着 dominant signal 走，注释为 Proliferating_Cell。

改成 specificity table 放在 DEGs 前面，标记为 PRIMARY EVIDENCE：

```
## PRIMARY EVIDENCE: Top specific markers for THIS cluster
specificity = this cluster's mean / max mean in other clusters.

  Oocyte | NLRP5: specificity=8.3x, 4.3% cells express
  Oocyte | LHX8: specificity=7.4x, 1.8% cells express
  ...

## SECONDARY EVIDENCE: Top 50 DEGs (Wilcoxon)
NOTE: DEGs ranked by statistical significance, NOT by cell-type specificity.

MKI67, TOP2A, CENPF, CDK1, ...
```

LLM 第一眼看到的是 NLRP5 8.3x，而不是 MKI67。Oocyte 回来了。

## 空间不连续的检测

v7 的 UMAP 上 Cluster 18 有 4 个断开的碎片。代码里有 `_check_cluster_continuity()` 做检测，算法是沿 UMAP 两个维度排序坐标，算相邻点的 gap，max_gap > 2.0 就标记为不连续。

但这个函数从来没触发过。原因是执行顺序：

```python
# Step 3: 质量分析
flagged = [r for r in quality_reports if r["should_filter"]]

# Early exit — 0 flagged 直接跳出
if not flagged:
    break  # ← 连续性检查在后面，永远执行不到

# Step 3.6: 连续性检查（被跳过了）
_check_cluster_continuity(...)
```

把连续性检查挪到 break 前面。跑了一版，检测到 10 个不连续 cluster。看数据：

```
Cluster 5:  1849+1 cells, max_gap=3.79  → 1个outlier触发
Cluster 9:  2221+3 cells, max_gap=2.83  → 3个outlier触发
Cluster 12: 1942+1 cells, max_gap=7.78  → 1个outlier触发
Cluster 18: 67+370 cells,  max_gap=4.25  → 真正不连续
```

大部分是 1-7 个 outlier cells 导致的假阳性。加了 `min_side_pct=0.05`：gap 两侧各需 >= 5% 细胞才判定不连续。过滤后只剩 Cluster 11（44+510）和 Cluster 18（67+370）。

## LLM 注释的几个问题

**被 dominant markers 劫持。** Top DEGs 里 40 个 cell cycle genes、10 个 oocyte-specific genes 时，LLM 顺着多数选了 Proliferating_Cell。prompt 里加了 Generic marker warning 明确列出哪些基因是 generic 的（MKI67/TOP2A/FOS/JUN/COL1A1），不能用来定义 cell type。

**Cell type 列表缺项。** Step 0 查卵巢返回 40 个类型，没有 Myofibroblast。Cluster 4 同时表达 SMC markers（MYH11/ACTA2）和 stromal markers（IGFBP5/SFRP1），验证了下：IGFBP5 在 Cluster 4 的 mean=3.8 是全 dataset 最高，SFRP1 91% 细胞表达而其他 SMC cluster 只有 10%。这不是 ambient RNA。被注释成 Smooth_Muscle_Cell 是因为列表里没有更合适的选项。Step 0 prompt 加了要求包含 transitional/mixed-phenotype types 后，查询结果从 40 增加到 50+。

**JSON 解析失败。** MiniMax-M3 偶尔返回 trailing comma、markdown fences、unterminated string。Step 0 失败返回 0 个 cell types，后续全废。`_call_anthropic` 加了 `max_retries=3` + 指数退避。实测 3 次重试基本能解决偶发问题。

## Marker 可视化

最初按 cluster 画 marker 图，每个 cluster 一张，top 3 markers。同一个 cell type 有多个 cluster 时（如 Endothelial 的 C1/C9/C12），看图要来回翻。

改成按 cell type 分组。同一 cell type 的所有 cluster markers 合并到一张图，LLM 返回的 5-10 个 key_markers 全部画出。每张图最后一个子图是 UMAP cluster reference——标注所有 cluster 编号，高亮当前 cell type 的 clusters，看 marker 表达时能直接对照 cluster 边界判断特异性。

## 迭代策略

auto mode 的迭代循环：聚类 → 注释 → 质检 → 不干净就过滤重跑。最多 max_iterations 轮。

连续性检查和分离检查会影响 resolution：不连续 → 提高 resolution 拆分；over-clustering → 降低 resolution 合并。但第一轮不调 resolution——用户的初始选择有生物学依据，QC 过滤改变数据分布后，下一轮再评估。

如果同时有质量 flagged clusters 和不连续 clusters，先过滤再处理不连续。因为过滤会改变数据分布，过滤后的连续性可能不一样。
