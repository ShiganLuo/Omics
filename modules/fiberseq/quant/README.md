# fiberseq/quant — Fiber-seq 区域定量与矩阵合并子模块

对 FIRE peak 区域逐样本定量四个 Fiber-seq 指标(可及性百分比、核小体占据、
m6A 密度、FIRE 占比),并合并为区域 × 样本矩阵;另支持单倍型(hap1/hap2)
拆分元素集的定量。C 形态子模块,共享父目录 `../fiberseq.yaml` 环境与
`modules/fiberseq/bin/quant_matrix.py` 脚本(run/merge 两个子命令)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 布局参数 | 形态 A 拼接的元素 BED:`{indir}/{sample_id}/{sample_id}.fire_peaks.bed`(FIRE 区域)、`.msp.bed.gz`、`.nuc.bed.gz`、`.m6a.bed.gz`、`.fire.bed.gz` |
| `hap_dir`(缺省 `indir`) | 单倍型拆分元素:`{hap_dir}/{sample_id}/{sample_id}.{hap1\|hap2}.msp.bed.gz` 等(由 modules/fibertools 的 ft_split_hap/ft_extract_hap 写在 phased 目录) |
| `samples` | 样本列表(merge 与聚合用) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.region_quant.tsv          # 逐区域四指标
{outdir}/{sample_id}/{sample_id}.{hap}.region_quant.tsv    # 单倍型版(hap1/hap2)
{outdir}/matrix/percent_accessible_matrix.tsv              # 可及性百分比矩阵
{outdir}/matrix/nuc_occupancy_matrix.tsv                   # 核小体占据矩阵
{outdir}/matrix/m6a_per_fiber_matrix.tsv                   # m6A 密度矩阵
{outdir}/matrix/fire_frac_matrix.tsv                       # FIRE 占比矩阵
```

日志:`{logdir}/{sample_id}/fiberseq_quant_run.log`、
`fiberseq_quant_hap_{hap}.log`、`{logdir}/group/fiberseq_quant_merge.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `common/4_elements` / `fiberseq_analysis/quant` / `log` | 目录 |
| `hap_dir` | — | `indir` | 单倍型元素目录 |
| `samples` | merge `-q` 列表 | `[]` | 样本列表 |
| `Params.fiberseq_quant.min_fibers` | — | `10` | json/schema 声明,当前 `quant.smk` 未消费 |
| `ROOT_DIR` | — | `.` | 拼 `modules/fiberseq/bin/quant_matrix.py` |
| 固定 | `run/merge` 子命令、`--sample`、`-l` | — | 见示例 |

环境:`conda: "../fiberseq.yaml"`。wildcard 约束:`hap = "hap[12]"`。

直接调用示例(等价引擎命令):

```bash
python modules/fiberseq/bin/quant_matrix.py run \
  -p s1.fire_peaks.bed -s s1.msp.bed.gz -u s1.nuc.bed.gz \
  -a s1.m6a.bed.gz -f s1.fire.bed.gz \
  -o out/s1 --sample s1 -l out/s1/run.log
python modules/fiberseq/bin/quant_matrix.py merge \
  -q out/s1/s1.region_quant.tsv,out/s2/s2.region_quant.tsv \
  -o out/matrix -l out/matrix/merge.log
```

## 2. 算法原理与生物学意义

对每个 FIRE peak 区域统计四类单分子指标:percent accessibility(纤维中
无核小体足迹/开放占比)、nucleosome occupancy(核小体占据占比)、m6A
density(每纤维 m6A 密度,裸露 DNA 量)、FIRE fraction(被判为 FIRE 调控
元件的占比)。merge 把逐样本 `region_quant.tsv` 拼成区域 × 样本矩阵,便于
跨样本/分组比较(如 aged vs youth 的调控元件可及性变化)。单倍型模式把
phased 纤维拆 hap1/hap2 后独立定量,可看等位基因特异的染色质组织。
结果怎么读:矩阵行 = 区域、列 = 样本;`percent_accessible_matrix.tsv`
高值 = 该区域开放;`fire_frac_matrix.tsv` 高值 = 该区域多为 FIRE 调控
元件;结合 `nuc_occupancy_matrix.tsv` 反向解读关闭/占位。

## 3. 参考文献

1. FIRE analysis pipeline — https://github.com/fiberseq/FIRE(FIRE 元件
   判定与 precision/q 约定)
2. fibertools-rs — https://github.com/fiberseq/fibertools-rs(元素抽取)
3. Peter CJ et al. *Single chromatin fiber profiling and nucleosome position
   mapping in the human brain*. Cell Reports Methods 4(12):100911, 2024.
   DOI: 10.1016/j.crmeth.2024.100911
