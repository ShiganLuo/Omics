# fiberseq/haplotype — 单倍型特异染色质可及性差异子模块

比较两条单倍型(hap1/hap2)在每个区域上的纤维可及性计数,做 Fisher
精确检验 + log2FC,找出单倍型特异的开放染色质区域。引擎脚本为父模块
共享的 `modules/fiberseq/bin/haplotype_diff.py`,conda 环境共享父目录
`../fiberseq.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `quant_dir`(config) | 上游 per-hap 区域定量表目录(纤维seq quant 步产出),缺省 `indir + "/quant"`;input 函数拼 `{quant_dir}/{sample_id}/{sample_id}.hap1.region_quant.tsv` 与 `.hap2.region_quant.tsv`(列:chrom start end name n_fibers n_accessible) |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.haplotype_diff.tsv    # 逐区域差异表
{outdir}/{sample_id}/{sample_id}.haplotype_summary.tsv # 检验区域数/显著数
```

`haplotype_diff.tsv` 列:`chrom start end name hap1_n_fibers
hap1_n_accessible hap1_frac hap2_n_fibers hap2_n_accessible hap2_frac
log2FC_hap1_vs_hap2 p_value q_value significant`。收尾规则
`fiberseq_haplotype_result` 聚合全部样本 diff 表。

### 1.3 参数

config 位于 `Params.fiberseq_haplotype` 下,与 `haplotype_diff.py` CLI 对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `min_fibers` | `-i/--min-fibers` | 10 | 每区域每单倍型最少纤维数,低于此跳过 |
| `alpha` | `-a/--alpha` | 0.05 | FDR(q 值)显著性阈值 |
| — | `-1/--hap1`, `-2/--hap2` | — | 两单倍型 region_quant.tsv |
| — | `-o` / `--sample` | — | 输出目录 / 样本名前缀 |
| — | `-l` | — | 日志文件 |
| — | `-n/--dry-run` | off | 只打印计划不产出 |
| threads / mem_mb | — | 4 / 8000 | `rule_res("fiberseq_haplotype_diff")` 配置 |

直接调用示例(脱离流程):

```bash
python modules/fiberseq/bin/haplotype_diff.py \
  -1 S1.hap1.region_quant.tsv -2 S1.hap2.region_quant.tsv \
  -o out/S1 --sample S1 -i 10 -a 0.05 -l haplotype_diff.log
```

## 2. 算法原理与生物学意义

对两单倍型共有区域,各建 2×2 列联表
`[[hap1 可及纤维, hap1 关闭纤维], [hap2 可及纤维, hap2 关闭纤维]]`,
做双侧 Fisher 精确检验;效应量为可及比例的 log2 倍数变化
`log2((acc1/n1)/(acc2/n2))`(加 1e-9 防除零);所有 p 值经
Benjamini-Hochberg 校正得 q 值,`q < alpha` 标记 `significant=yes`。
两侧任一 `n_fibers < min_fibers` 的区域不检验。

生物学意义:等位基因特异的染色质可及性(allelic imbalance)是顺式
调控变异(cis-regulatory variant)的直接读出——同一细胞核内两条单倍型
共享 trans 环境,差异必然来自 cis 层面(如增强子 motif 破坏)。结果
读法:先看 `haplotype_summary.tsv` 的显著区域数,再按 `log2FC` 与
`q_value` 排序取 top 区域做下游 motif/变异富集。

## 3. 参考文献

1. FIRE analysis pipeline — https://github.com/fiberseq/FIRE(单倍型分析
   约定与区域定量格式)。
2. 引擎脚本 `bin/haplotype_diff.py` 为本仓库自研(共享 `fiberseq/bin/`,
   见父模块 `fiberseq/README.md`;Fisher 检验用 scipy.stats.fisher_exact,
   FDR 用 BH 校正)。
