# fiberseq/telocen_norm — 端粒/着丝粒区域深度归一化比较（fig2）

对两个样本组（如 youth / aged）在端粒-着丝粒区域（BED 第 4 列分类：
`cen_core` / `cen_flank` / `subtelomere` / `telomere_core`）上计算
深度归一化的 m5CpG 与 m6A 密度，并做区域配对比较，输出 fig2 风格双面板图。

来源：`output/luanchao_mulatta/Fiberseq/conclusion_figures/`
的 `compute_telocen_normalized_mp.py` + `plot_fig2_telocen_v5.py`
（已参数化，无硬编码路径/样本）。

## 规则

| rule | 输出 |
|---|---|
| `fiberseq_telocen_norm_run` | `region_normalized.tsv`（每区域×样本：aligned_bp、m5CpG/m6A marks、每 Mb 密度）、`class_diff_normalized.tsv`（每类×指标：组均值/SD、diff、log2fc、配对 Wilcoxon p） |
| `fiberseq_telocen_norm_fig` | `telocen_normalized.<fmt>`（左：m5CpG 密度；右：Hia5-m6A 密度；配对连线 + 星号） |

## 输入配置

```json
"inputs": {"telocen_norm": {"bam": {"<sample>": "<aligned .fiberseq.fire.sorted.bam>"}}},
"Params": {"fiberseq_telocen_norm": {"groups": {"youth": ["S1"], "aged": ["S2"]}}},
"genome": {"telocen_bed": "/path/Mmul_T2T.telocen_regions.bed"}
```

- `groups` 必须恰好 2 组（键=组标签，值=样本列表，顺序：基线组在前）；
  多样本组按区域取组内均值后再配对。
- `genome.telocen_bed`：4 列 BED（chrom start end region_class），
  Mmul_T2T 注释见 `/data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/Mmul_T2T.telocen_regions.bed`。
- `Params.fiberseq_telocen_norm.fiberseq_telocen_norm_fig.formats` 控制图片格式。

## 统计口径（诚实说明）

每组 n=1 个体时，Wilcoxon p 是"跨区域方向一致性"，不是群体显著性。
图注/解读时必须保留该限定。

## Pitfalls

- 比较是严格的两组配对设计；`groups` 不足 2 组时 run rule 报友好 ValueError。
- `mm_positions` 按 modules.md #0 自包含在 bin 脚本内（模块不得引用其他
  模块的脚本）。
- 空覆盖区域（aligned_bp=0）密度记 NA，不进配对统计。
