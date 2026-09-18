# scRNAseq 自动注释：从问题到设计

单细胞 RNA-seq 的细胞类型注释，传统做法是人工看 marker genes、查文献、画 UMAP，来回折腾几天。我们做了一个全自动流程——聚类、AI 注释、质量检查、迭代优化，一轮跑完。但做下来发现，让 AI 注释准确比想象中难得多。这篇文章记录我们踩过的坑和背后的设计思考。

## 一、Wilcoxon DEG 的陷阱

最开始，我们把 Wilcoxon rank-sum test 排出的 top 50 DEGs 直接丢给 LLM 注释。看起来很直觉：哪些基因在这个 cluster 里显著高表达，就用它们判断细胞类型。

然后 Oocyte 给了我们一记暴击。

Oocyte cluster 的 top 50 DEGs 全是 cell cycle genes——MKI67、TOP2A、CENPF、CDK1。LLM 看了一眼，注释为 "Proliferating_Cell"。但这个 cluster 里真正定义它身份的是 DDX4、NLRP5、FIGLA 这些卵母细胞特异基因，它们只在 2-4% 的细胞里低表达，Wilcoxon 排名根本排不进去。

问题的本质：**Wilcoxon 按统计显著性排名，不按生物学特异性排名。** 一个基因在 oocyte 里 50% 细胞高表达、其他 cluster 里几乎不表达（如 MKI67），它的 Wilcoxon score 碾压一个只在 2% 细胞里表达但完全特异的基因（如 NLRP5）。

**解法：引入 cross-cluster specificity。**

在注释之前，先算每个 tissue-specific marker 的 specificity：

```
specificity = this_cluster_mean / max(other_clusters_mean)
```

NLRP5 在 oocyte 里 mean=0.082，其他 cluster 最高 mean=0.010，specificity=8.3x。虽然表达量低，但它是真正的 oocyte marker。

**关键设计：specificity table 放在 prompt 最前面，作为 PRIMARY EVIDENCE。** DEGs 降为 SECONDARY EVIDENCE。LLM 第一眼看到的是 "NLRP5: specificity=8.3x"，而不是 "MKI67, TOP2A, CENPF..."。

这个改变让 Oocyte 从 "Proliferating_Cell" 变回了 "Oocyte (high confidence)"。

## 二、聚类空间不连续

UMAP 上看到 Cluster 15 和 18 很分散，一开始以为是聚类参数问题。但仔细查发现更严重——它们不是分散，是**空间不连续**。Cluster 18 有 4 个断开的碎片，散布在 UMAP 不同位置。

代码里其实有 `_check_cluster_continuity()` 函数检测不连续 cluster，但它从来没触发过。为什么？

因为它的执行位置在 Step 3.6——**在 early exit 之后**。流程是：

```
Step 3: 质量分析 → 0 flagged
→ if not flagged: "All clusters clean" → break  ← 直接跳出
→ Step 3.5 / 3.6 永远执行不到
```

Cluster 18 的 UMAP max_gap=4.25（阈值2.0），本来应该被检测到，但代码直接跳过了。

**解法：把连续性检查移到 early exit 之前。**

但这还不够。第一次跑的时候，连续性检查触发了，但检测到 10 个"不连续" cluster——大部分只有 1-7 个 outlier cells。比如 Cluster 5 有 1849 个主群 + 1 个离群细胞，就被判定为不连续，导致 resolution 被错误提升。

**解法：增加 `min_side_pct` 过滤。** gap 两侧各需 >= 5% 细胞才判定为真正不连续。这样 1 个 outlier 不会触发假阳性，但 Cluster 18 的 67+370 细胞（15%+85%）能正确检测。

## 三、LLM 注释的局限性

### 被 generic markers 误导

Step 0 查询细胞类型时，如果列表里有 "Proliferating_Cell"，LLM 就会把 oocyte 注释成它——因为 cell cycle genes 占据了 top DEGs。如果没有这个选项，LLM 只能选 "Oocyte"。

这揭示了一个根本问题：**LLM 倾向于被 dominant markers 主导**。如果 top DEGs 里 40 个是 cell cycle genes、10 个是 oocyte-specific，LLM 会选 Proliferating_Cell 而不是 Oocyte。

我们在 prompt 里加了 Generic marker warning：

```
Cell cycle genes (MKI67, TOP2A, CENPF...)、stress genes (FOS, JUN, EGR1...)
和 ECM genes (COL1A1, COL1A2, DCN) 出现在多种细胞类型中。如果这些基因
占据 top DEGs，但 specificity table 显示 tissue-specific markers，
必须基于 tissue-specific markers 注释。
```

### 只看 top DEGs，不看 specificity table

即使把 specificity table 放在 prompt 里，LLM 有时候还是不看。特别是当 top DEGs 有 50 个而 specificity table 只有 5 个时，LLM 觉得这不是"tie"，不需要看 specificity。

**解法：把 specificity table 标记为 PRIMARY EVIDENCE，DEGs 标记为 SECONDARY EVIDENCE。** 在 Rules 第一条写 "SPECIFICITY FIRST"。这样 LLM 被明确告知：先看 specificity，再看 DEGs。

### Cell type 列表不完整

Step 0 查询卵巢的细胞类型，返回了 40 个，但没有 Myofibroblast。结果 Cluster 4（同时表达 SMC markers 和 stromal markers）被注释为 "Smooth_Muscle_Cell"——因为列表里没有更合适的选项。

**解法：在 Step 0 prompt 里明确要求包含 transitional/activated/mixed-phenotype 类型。** 加了一句：

```
IMPORTANT: Include transitional, activated, and mixed-phenotype cell types
(e.g. Myofibroblast, Activated_Stromal_Cell, Tip_Endothelial_Cell,
Proliferating_cell). These express markers from TWO or more canonical
lineages simultaneously.
```

查询结果从 40 个增加到 50+ 个。

## 四、Marker 可视化的问题

### 按 cluster 还是按 cell type？

最初 marker 图按 cluster 输出——每个 cluster 一张图，画 top 3 markers。但同一个 cell type 可能有多个 clusters（如 Endothelial 有 C1/C9/C12），看图时要来回翻才能理解一个 cell type 的完整 marker profile。

**改成按 cell type 分组。** 同一 cell type 的所有 clusters markers 合并到一张图，标题标注对应的 cluster 编号。如 `Endothelial_Cell [C1, C9, C12]`。

### 只画 top 3 还是全部？

LLM 返回 5-10 个 key_markers，但最初只画 top 3。用户看不到完整的 marker profile。

**改成画全部 markers。** 每个 cell type 的图包含 N 个 marker 表达子图 + 1 个 UMAP cluster reference 子图。

### 缺少 cluster 边界参考

看 marker 表达图时，不知道表达的位置对应哪个 cluster。需要一张 UMAP 聚类图做参照。

**在每张 marker 图的最后加一个 UMAP cluster reference 子图。** 标注所有 cluster 编号，高亮当前 cell type 的 clusters。

## 五、LLM 调用的稳定性

MiniMax-M3 偶尔返回格式错误的 JSON——trailing comma、markdown fences、unterminated string。Step 0 查询失败后返回 0 个 cell types，后续所有注释都失去参考。

**解法：`_call_anthropic` 增加 `max_retries=3`。** JSON 解析失败时指数退避重试。实测中 3 次重试基本能解决偶发的格式问题。

## 六、迭代优化策略

auto mode 不是一次跑完，而是迭代循环：

1. **聚类 → 注释 → 质量检查**
2. 如果有 flagged clusters → 过滤 → 重新聚类
3. 如果有不连续 clusters → 提高 resolution → 重新聚类
4. 如果有 over-clustering → 降低 resolution → 重新聚类
5. 最多跑 `max_iterations` 轮

每一轮的 resolution 可能变化。第一轮尊重用户的初始 resolution，后续轮次根据连续性和分离检查结果自动调整。

关键保护：**第一轮不调 resolution。** 因为用户的初始选择有生物学依据，QC 过滤后的数据可能需要不同的 resolution。

## 七、从 v7 到 v10 的演进

| 版本 | 问题 | 修复 |
|---|---|---|
| v7 | Cluster 18 不连续未检测到 | 连续性检查移到 early exit 前 |
| v8 | 10 个 cluster 被误判为不连续 | 增加 min_side_pct=0.05 过滤 |
| v9 | Oocyte 被注释为 Proliferating_Cell | Step 0 增加过渡态类型 |
| v10 | LLM 被 generic DEGs 误导 | Specificity table 升级为 PRIMARY EVIDENCE |
| v10 | Marker 图按 cluster 输出 | 改为按 cell type 分组 + UMAP reference 子图 |
