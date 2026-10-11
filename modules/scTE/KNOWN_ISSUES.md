# scTE 已知问题：细胞数与 Cell Ranger 不一致（--expect-cells 硬截断）

记录时间：2026-10-10
样本案例：output/luancao/scRNAseq/common/3_raw_h5ad/luanchao-11238-10XSC3

## 问题现象

scTE 输出的细胞数不等于 Cell Ranger 调用的细胞数。luanchao-11238-10XSC3 样本中
scTE_out.csv.gz 正好 10,000 个细胞（19,141 个 TE/基因行），而
luanchao-11238-10XSC3_cellranger.h5ad 有 9,792 个细胞。两者并非包含关系：

- 404 个 Cell Ranger 调用的真细胞被 scTE 丢弃
- 612 个 Cell Ranger 未调用的 barcode（空滴/低质量）反而混入 scTE 矩阵
- scTE_out.stat 显示通过初筛（≥400 reads）的 barcode 共 18,907 个，
  其中约 8,900 个被截断丢弃

## 根因

scTE 的 `--expect-cells` 参数（默认 10000）在 count_expression() 中是硬截断，
不是 Cell Ranger 式的期望值校正：

```python
sortcb = sorted(whitelist.items(), key=读数, reverse=True)  # 按 reads 数降序
for n, k in enumerate(sortcb):
    if k[1] < genenumber: break   # --min_genes，默认 200
    if n >= cellnumber: break     # --expect-cells，默认 10000，到此为止
```

即：按 reads 数取 top N barcode。reads 数 ≠ UME 数 ≠ 真细胞，高 RNA 空滴、
死细胞、双胞体都排在前面，所以既截断了数量，选出来的细胞也不对。

## 作者设计意图（He et al., Nat Commun 2021, 12:1456）

scTE 定位为独立 pipeline，内置 barcode demultiplexing + quality filtering，
不依赖 Cell Ranger 的 barcodes.tsv.gz（兼容 STARsolo、Smart-seq/C1 等）。
10000 是工程默认值：作者自己的数据集单样本多 ≤1 万细胞，从未触发截断，
论文也从未验证过细胞调用与 Cell Ranger 的一致性（只验证了基因表达
Pearson>0.95 与 UMAP 拓扑）。社区早有同款反馈且无人维护：

- GitHub issue #96 "The number of cells is not consistent with the result of
  cellranger"：scTE 1,545 vs Cell Ranger 7,958
- GitHub issue #108 同题（无回复）
- 其他：#88 filterCRs IndexError、#106 超 3000 细胞报错、#117 IndexError

## 修复方向（未实施）

1. 精确修复：`--whitelist barcodes.tsv.gz` 喂 Cell Ranger 的 filtered
   barcodes，同时 `--expect-cells` 调大（如 50000，count_expression 的截断
   对 whitelist 同样生效）
2. 粗放修复：仅 `--expect-cells 20000` 以上，覆盖全部通过初筛的 barcode
   （本样本约 18,907 个），下游再自行过滤低质量滴

注意：封装脚本 modules/scTE/bin/scTE_quantify.py 目前未透传 --expect-cells
和 --whitelist，run_scTE() 的 cmd 需要加参数并经 scTE.smk 的 Params.scTE
配置暴露。

## 同赛道替代工具（交叉验证用）

- IRescue (2024)：亚家族水平，EM 分配多比对 + 规范 UMI 去重（含错误校正），
  专为 3' tag-based 数据设计，同赛道精度优于 scTE/SoloTE（自评测）
- SoloTE (2022, Commun Biol)：位点水平（unique reads），轻量 ≤15GB，
  2026 benchmark 推荐精确优先方案
- Stellarscope (2025, Cell Rep Methods)：位点水平贝叶斯混合模型 + EM，
  pooled 拟合降稀疏性
- 2026 benchmark (Finazzi et al., bioRxiv 2026.02.26.708244) 结论：
  年轻 TE 位点本质不可分辨，建议亚家族水平聚合；unique-mapper 优于
  EM 多比对；gene-TE 归属混淆必须显式检查
