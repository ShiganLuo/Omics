# fiberseq/epi_drift — 单纤维 m5CpG 表观漂移检验（fig9）

在随机窗口上比较两组样本的"每纤维 m5CpG 比例"方差（表观漂移），
输出窗口级 TSV、汇总 JSON 与 fig9 风格双面板图（无标题版）。

来源：`output/luanchao_mulatta/Fiberseq/conclusion_figures/`
的 `plot_fig9_drift.py`（数据部分）+ `plot_fig9_replot.py`（绘图部分）
（已参数化，无硬编码路径/样本）。

## 规则

| rule | 输出 |
|---|---|
| `fiberseq_epi_drift_run` | `epigenetic_drift_windows.tsv`（chrom/start/end、两组 mean/var）、`epigenetic_drift.json`（groups、n、signs、sign_test_p、median_var_ratio、params） |
| `fiberseq_epi_drift_fig` | `epigenetic_drift.<fmt>`（左：log2 方差比直方图 + 中位线；右：drift vs mean shift 散点） |

## 输入配置

```json
"inputs": {"epi_drift": {"bam": {"<sample>": "<aligned .fiberseq.fire.sorted.bam>"}}},
"Params": {"fiberseq_epi_drift": {
    "groups": {"youth": ["S1"], "aged": ["S2"]},
    "fiberseq_epi_drift_run": {"chrom_prefix": "NC_"}}},
"genome": {"fai_index": "/path/genome.fa.fai"}
```

- `groups` 必须恰好 2 组（键=组标签，值=样本列表；第二组是"被比较"组，
  方差比 = 第二组 / 第一组）。多样本组的纤维合并后再抽样。
- 窗口参数：`win`(5000)、`n_windows`(600)、`min_fibers`(20)、
  `subsample`(20)、`seed`(42)、`min_chrom_len`(30Mb)、
  `chrom_prefix`（NCBI 登录号前缀 `NC_`；Ensembl 命名的参考填 `""`）。
- 每纤维 m5CpG 比例 = 窗口内 m5CpG marks / 覆盖碱基数；覆盖 <50% 窗口长
  的纤维剔除；两组均 ≥min_fibers 才保留，之后各抽 subsample 条。

## RNG 语义

窗口抽样用 `--seed` 的全局 RNG；每窗口的纤维子抽样用
`Random(f"{seed}:{chrom}:{start}")`，与迭代顺序无关（原一次性脚本用全局
RNG shuffle，结果依赖顺序）。同 seed 可复现。

## Pitfalls

- 方差比方向：signs 计数的是"第二组方差 > 第一组方差"的窗口数；
  图 x 轴标签自动携带组名，勿再手写 aged/youth。
- 大 BAM 逐窗口 fetch 较慢（600 窗口串行）；属 IO 密集，threads 参数预留。
