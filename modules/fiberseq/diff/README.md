# fiberseq/diff — 区域可及性差异分析子模块

对 Fiber-seq 定量矩阵(区域 × 样本)做组间差异可及性检验
(Mann-Whitney U + log2FC + BH-FDR),按区域类别分层汇总,输出差异表、
火山图与基因表(后者供 modules/function 富集)。属 fiberseq 家族
(布局 C),共享父目录 `../fiberseq.yaml` 环境与 `bin/` 脚本;
方法背景与家族输入约定见 [fiberseq/README.md](../README.md)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `quant_dir` / `indir` (config) | 定量矩阵目录:`{quant_dir or indir/quant}/matrix/percent_accessible_matrix.tsv`(上游 fiberseq/quant 产出,区域 × 样本) |
| `annot_dir` / `indir` (config) | peak 注释目录:`{annot_dir or indir/annotate}/{annot_sample}/{annot_sample}.peaks_annotated.tsv`(上游 fiberseq/annotate 产出) |
| `Params.fiberseq_diff.annot_sample` | 取哪个样本的注释表;缺省 `samples[0]` |
| `Params.fiberseq_diff.group_spec` | 组设计串 `name:gA:s1+s2\|gB:s3,...`(meta 推导,优先于脚本的 `--design` TSV) |
| wildcards `contrast` | 对比名(如 `aged_vs_youth`);`-C` 限定只跑该对比 |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 输出(默认 `fiberseq_analysis/diff`)/日志/项目根 |

### 1.2 输出

```
{outdir}/{contrast}/diff_{contrast}.tsv                          # 逐区域差异统计
{outdir}/{contrast}/volcano_{contrast}.pdf                       # 火山图
{outdir}/{contrast}/class_summary_{contrast}.tsv                 # 区域类别分层汇总
{outdir}/{contrast}/TEcount/{contrast}.TEcount_Gene.name.tsv     # 基因表(供 modules/function)
```

回执:`{outdir}/{contrast}/fiberseq_diff_run_{contrast}.sh`(另有
`{logdir}/cmd/...` 时间戳历史,见父 README 回执双写约定);日志
`{logdir}/group/fiberseq_diff_run_{contrast}.log`。

### 1.3 参数

config 位于 `Params.fiberseq_diff` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `quant_dir` | `-m` 矩阵路径 | `indir/quant` | 矩阵目录 |
| `annot_dir` | `-a` 注释路径 | `indir/annotate` | 注释目录 |
| `Params.fiberseq_diff.group_spec` | `--group-spec` | "" | 对比设计串(`name:gA:s1+s2\|gB:s3,...`),覆盖 `--design` |
| `Params.fiberseq_diff.annot_sample` | — | `samples[0]` | 注释表取样 |
| `Params.fiberseq_diff.alpha` | `-p` | 0.05 | BH-FDR 阈值 |
| (contrast) | `-C` | wildcard | 只跑指定对比(缺省脚本跑 design 里全部) |
| (threads) | — | 4 | `rule_res("fiberseq_diff_run")` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python modules/fiberseq/bin/diff_accessibility.py \
    -m fiberseq_analysis/quant/matrix/percent_accessible_matrix.tsv \
    --group-spec "aged_vs_youth:aged:s1+s2|youth:s3+s4" \
    -a fiberseq_analysis/annotate/s1/s1.peaks_annotated.tsv \
    -C aged_vs_youth \
    -o fiberseq_analysis/diff/aged_vs_youth \
    -p 0.05 -l log/group/fiberseq_diff_run_aged_vs_youth.log
```

## 2. 算法原理与生物学意义

`diff_accessibility.py` 对矩阵中每个区域做两组比较:

- **检验**:Mann-Whitney U(非参数,不假设正态;纤维级可及性比例
  分布偏态时稳健);
- **效应量**:log2FC(组间可及性比例之差的对数倍数);
- **多重检验**:Benjamini-Hochberg FDR,`-p`/alpha 决定显著阈值;
- **分层汇总**:给定 `peaks_annotated.tsv` 时按区域类别(如
  promoter/enhancer,或 Nat Genet 2025 的 CDR core vs flanking)汇总
  显著差异的分布。

结果怎么读:`diff_{contrast}.tsv` 每行一个区域(统计量、log2FC、
FDR),火山图纵轴 -log10(FDR)、横轴 log2FC,右上/左上为显著开放/
关闭区域;`class_summary` 回答"差异集中在哪类调控元件";基因表把
差异区域映射到基因,直接作为 modules/function(GO/KEGG/GSEA)输入。
生物学意义:比较不同条件(aged/youth、处理/对照)的染色质可及性
重编程,差异开放区域指向条件特异的调控元件与转录因子足迹变化。

## 3. 参考文献

1. Peter CJ et al. *Single chromatin fiber profiling and nucleosome
   position mapping in the human brain*. Cell Reports Methods 4(12):100911,
   2024. DOI: 10.1016/j.crmeth.2024.100911.
2. Genome Res 2024 — DOI: 10.1101/gr.279095.124.
3. Nat Genet 2025 — DOI: 10.1038/s41588-025-02324-w(CDR 分层方法出处;
   三条引文与 fiberseq 家族 README 的 Sources 一致)。
