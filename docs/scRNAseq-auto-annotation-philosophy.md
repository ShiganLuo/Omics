# scRNAseq 自动注释的坑

做了一套全自动的单细胞注释流程，聚类到注释到质检一轮跑完。跑通不难，注释准确才是大坑。记录一下踩过的问题。

## Oocyte 去哪了

最早的版本，把 Wilcoxon top 50 DEGs 直接扔给 LLM 注释。跑完一看，Oocyte 没了——被注释成了 Proliferating_Cell。

翻了下 DEGs，全是 MKI67、TOP2A、CENPF 这些 cell cycle genes。Oocyte 在减数分裂，转录本里 cell cycle 信号确实强。但真正定义 oocyte 身份的 DDX4、NLRP5、FIGLA 只在 2-4% 的细胞里低表达，Wilcoxon 排不进 top 50。

Wilcoxon 是按统计显著性排的，不是按生物学特异性排的。MKI67 在 oocyte 里 50% 高表达、其他 cluster 几乎不表达，Wilcoxon score 天然碾压 NLRP5。这是方法本身的局限，不是 bug。

加了一层 cross-cluster specificity 计算：`specificity = this_cluster_mean / max(other_clusters_mean)`。NLRP5 在 oocyte 里 mean=0.082，其他 cluster 最高 0.010，specificity 8.3x。虽然表达量低，但够特异。

然后把 specificity table 放到 prompt 最前面，标记为 PRIMARY EVIDENCE，DEGs 降为 SECONDARY。LLM 第一眼看到的是"NLRP5: 8.3x"而不是"MKI67, TOP2A, CENPF..."。Oocyte 就回来了。

但 marker 图上 oocyte 的表达还是很弱——只有 2-4% 细胞有信号，UMAP 上看起来到处都有零星红点。这是生物学事实，oocyte-specific genes 就是低频低表达。specificity ratio 高是因为其他 cluster 几乎不表达，不是因为 oocyte 里表达高。

## 空间不连续

UMAP 上看到 Cluster 15 和 18 很分散，一开始以为是聚类参数问题。仔细查了下，比分散严重——Cluster 18 有 4 个断开的碎片，散布在 UMAP 不同位置。

代码里有 `_check_cluster_continuity()` 检测不连续 cluster，但从来没触发过。看了下执行顺序，发现问题：

```
Step 3: 质量分析 → 0 flagged
→ if not flagged: break  ← 直接跳出
→ 连续性检查在 break 后面，永远执行不到
```

把连续性检查挪到 break 前面。跑了一版，检测到 10 个"不连续" cluster。大部分只有 1-7 个 outlier cells——Cluster 5 有 1849 个主群 + 1 个离群细胞，就被判定不连续，resolution 被错误提升。

加了个 `min_side_pct` 过滤，gap 两侧各需 >= 5% 细胞才算真正不连续。1 个 outlier 不触发，但 Cluster 18 的 67+370 细胞（15%+85%）能正确检测到。

## LLM 的毛病

LLM 注释有几个坑。

第一个是被 generic markers 误导。Step 0 查卵巢细胞类型时，如果列表里有 Proliferating_Cell，LLM 就把 oocyte 注释成它。因为 cell cycle genes 占了 top DEGs 的大多数，LLM 顺着 dominant signal 选了。如果列表里没有这个选项，LLM 只能选 Oocyte。所以加了过渡态类型（Myofibroblast、Activated_Stromal_Cell 等）的同时，也要在 prompt 里警告 generic markers 不能定义细胞类型。

第二个是 LLM 不看 specificity table。即使放在 prompt 里，top DEGs 有 50 个而 specificity 只有 5 个时，LLM 觉得这不是 tie，不需要看后面。把 specificity 标记为 PRIMARY、DEGs 标记为 SECONDARY 后好了不少，但偶尔还是会忽略。

第三个是 cell type 列表不完整。Step 0 查卵巢返回 40 个类型，没有 Myofibroblast。Cluster 4 同时表达 SMC markers（MYH11, ACTA2）和 stromal markers（IGFBP5, SFRP1, DCN），被注释成 Smooth_Muscle_Cell。查了下，IGFBP5 在 Cluster 4 的表达是全 dataset 最高的，SFRP1 91% 细胞表达而其他 SMC cluster 只有 10%——这不是 ambient RNA，是真实的双重表型。在 Step 0 prompt 里加了要求包含 transitional/mixed-phenotype 类型后，查询结果从 40 增加到 50+，Myofibroblast 出来了。

LLM 返回的 JSON 也偶尔出错——trailing comma、markdown fences、unterminated string。Step 0 查询失败返回 0 个 cell types，后面全废。加了 3 次 retry + 指数退避，基本能解决偶发问题。

## Marker 图怎么画

最初 marker 图按 cluster 画，每个 cluster 一张，top 3 markers。同一个 cell type 有多个 cluster 时（如 Endothelial 的 C1/C9/C12），要来回翻才能看全。

改成按 cell type 分组，同一类型的所有 cluster markers 合并到一张图。LLM 返回 5-10 个 key_markers，全部画出来，不再截断 top 3。每张图最后加一个 UMAP cluster reference 子图，标注所有 cluster 编号，高亮当前 cell type 的 clusters——看 marker 表达时能直接对照 cluster 边界。

## 迭代逻辑

auto mode 不是一轮跑完。聚类、注释、质检之后，如果有 flagged clusters 就过滤重跑，有不连续 clusters 就提高 resolution 重跑，有 over-clustering 就降低 resolution 重跑。最多跑 max_iterations 轮。

第一轮不调 resolution——用户的初始选择有生物学依据，QC 过滤后的数据可能需要不同的 resolution。后续轮次再根据检查结果自动调整。

连续性检查里有个细节：如果同时有质量 flagged clusters 和不连续 clusters，先过滤再处理不连续。因为质量过滤会改变数据分布，过滤后的连续性可能不一样。
