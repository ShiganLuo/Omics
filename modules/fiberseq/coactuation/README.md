# fiberseq/coactuation — 同纤维共激活(co-actuation)分析子模块

C 形态子模块(单规则 `fiberseq_coactuation_run` + result 锚点):对 FIRE
peak 区域对做"同一 DNA 分子上共同激活"检验(Fisher 精确检验 + 实测 vs
独立期望)。共享父目录环境 `../fiberseq.yaml`,引擎脚本
`modules/fiberseq/bin/coactuation.py`(经 `ROOT_DIR` 引用)。

## 1. 模块输入输出与参数

### 1.1 输入

元素文件约定见父目录 `modules/fiberseq/README.md`(indir = `common/4_elements`):

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}.fire_peaks.bed` | 候选区域/FIRE peak BED |
| `{indir}/{sample_id}/{sample_id}.fire.bed.gz` | 逐纤维 FIRE 元件 BED(`ft fire --extract`) |
| `{indir}/{sample_id}/{sample_id}.msp.bed.gz` | 逐纤维 MSP BED12 |
| `{indir}/{sample_id}/{sample_id}.nuc.bed.gz` | 逐纤维核小体 BED12 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.coactuation.tsv          # 区域对明细表
{outdir}/{sample_id}/{sample_id}.coactuation_summary.tsv  # 汇总(3 个指标)
{outdir}/{sample_id}/fiberseq_coactuation_run_{sample_id}.sh  # 回执(无时间戳,覆盖式)
logdir/{sample_id}/fiberseq_coactuation_run.log           # 运行日志
```

明细表列:`chrom site1 name1 site2 name2 n_covering_fibers n_fire_site1
n_fire_site2 n_both n_site1_only n_site2_only n_neither actual_frac
expected_frac diff p_value q_value significant`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.fiberseq_coactuation.max_distance` | `-d` | 50000 | 区域对最大间隔(bp) |
| `Params.fiberseq_coactuation.min_fibers` | `-i` | 10 | 区域对最少共同覆盖纤维数 |
| `Params.fiberseq_coactuation.max_pairs` | `-m` | 200000 | 受检区域对上限 |
| `Params.fiberseq_coactuation.alpha` | `-a` | 0.05 | FDR 显著性阈值 |
| `samples` | — | [] | 样本列表(result 锚点 expand) |
| `indir` / `outdir` / `logdir` | — | `common/4_elements` / `fiberseq_analysis/coactuation` / `log` | 目录基路径 |
| `ROOT_DIR` | — | `.` | 仓库根(拼 bin/ 脚本路径) |

直接调用示例(引擎脚本,模块外独立可用):

```bash
python modules/fiberseq/bin/coactuation.py \
  -p sample.fire_peaks.bed -f sample.fire.bed.gz -s sample.msp.bed.gz -u sample.nuc.bed.gz \
  -o out/ --sample s1 -d 50000 -i 10 -m 200000 -a 0.05
```

## 2. 算法原理与生物学意义

对同染色体、间隔 ≤ `max_distance` 的每对 peak:取同时覆盖两个区域的
纤维集合 `both_cov`(逐纤维去重,`msp ∪ nuc` 记为覆盖),其中在区域 1/2
携带 FIRE 元件的纤维集为 f1/f2,构造 2×2 表
`[[b, a1_only], [a2_only, n]]`(b = 两处都有 FIRE),做单侧 Fisher 精确
检验;并计算

      actual   = |f1 ∩ f2| / |both_cov|
      expected = (|f1| / |both_cov|) · (|f2| / |both_cov|)
      diff     = actual − expected

p 值经 BH-FDR 得 q 值,`q < alpha` 记 significant。生物学意义:若两个
调控元件在同一 DNA 分子上共同激活的频率显著高于独立期望(diff > 0 且
显著),提示二者可能受同一调控机制协调(如染色质域共开放),是单分子
分辨率下顺式调控协同的证据。

## 3. 参考文献

1. fiberseq 家族方法来源(含 co-actuation 独立性期望比较的方法出处):
   见 `modules/fiberseq/README.md` Sources —
   https://doi.org/10.1016/j.crmeth.2024.100911 、
   https://doi.org/10.1101/gr.279095.124 、
   https://doi.org/10.1038/s41588-025-02324-w
2. FIRE analysis pipeline — https://github.com/fiberseq/FIRE
