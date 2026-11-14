# fiberseq/footprint — 单分子 TF 足迹(Fiber-seq footprint)子模块

围绕 motif 位点做单分子转录因子足迹分析:把每条覆盖 motif 的 DNA 纤维分入
三类(MSP 且 motif 核心无 m6A 足迹 / MSP 无足迹 / 核小体占据),按类聚合
m6A 比例轮廓,并对位点足迹得分做 Kolmogorov-Smirnov 检验。C 形态子模块,
共享父目录 `../fiberseq.yaml` 环境与 `modules/fiberseq/bin/footprint_profile.py`
脚本。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `motif_dir`(缺省 `indir/motif`) | motif 位点 BED:`{MOTIF_DIR}/{sample_id}/{sample_id}.motif_sites.bed`(上游 fiberseq motif 模块产出) |
| `indir` + 布局参数 | 形态 A 拼接的元素 BED:`{indir}/{sample_id}/{sample_id}.m6a.bed.gz`、`.msp.bed.gz`、`.nuc.bed.gz`(fibertools 元素抽取产出;BED12 块结构,多纤维记录) |
| `samples` | 样本列表(result 聚合用) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.footprint_categories.tsv   # 每纤维×位点 三类归属
{outdir}/{sample_id}/{sample_id}.footprint_profiles.tsv     # 每类 m6A 比例轮廓
{outdir}/{sample_id}/{sample_id}.footprint_scores.tsv       # 位点足迹得分 + KS 检验
{outdir}/{sample_id}/{sample_id}.footprint_profile.pdf      # 轮廓图
```

日志:`{logdir}/{sample_id}/fiberseq_footprint_run.log`。
规则注释:回执脚本在 log/cmd 留带时间戳副本,outdir 下另有无时间戳副本
(记录当前命令)。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `common/4_elements` / `fiberseq_analysis/footprint` / `log` | 目录 |
| `motif_dir` | `-m` | `indir/motif` | motif 位点目录 |
| `samples` | — | `[]` | 样本列表 |
| `Params.fiberseq_footprint.flank` | `-f` | `1000` | motif 侧翼窗口(bp) |
| `ROOT_DIR` | — | `.` | 拼 `modules/fiberseq/bin/footprint_profile.py` |
| 固定 | `--sample <id>`、`-l <log>` | — | 逐样本调用 |

脚本另支持 `--motif-filter`(如 `CTCF@,MYC@`,按 motif 名过滤)。
环境:`conda: "../fiberseq.yaml"`。

直接调用示例(等价引擎命令):

```bash
python modules/fiberseq/bin/footprint_profile.py \
  -m s1.motif_sites.bed -a s1.m6a.bed.gz -s s1.msp.bed.gz -u s1.nuc.bed.gz \
  -o out/s1 --sample s1 -f 1000 -l out/s1/run.log
```

## 2. 算法原理与生物学意义

Fiber-seq 单分子读出核小体足迹与调控元件占用(见 fiberbrowser 模块),本
模块把 motif 位点(如 CTCF)上每条纤维分三类 [1]:(i) MSP 覆盖且 motif 核心
有 m6A-free 足迹——TF 占位;(ii) MSP 覆盖无足迹——开放但未占位;
(iii) 核小体占据——关闭。按类聚合窗口内 m6A 比例轮廓(m6A 是"裸露 DNA"
标记,足迹 = 核心 m6A 缺口),并以位点级足迹得分做类别间 KS 检验,量化
TF 占位的统计显著性。结果怎么读:`*_footprint_scores.tsv` 中 KS p 值小 =
两类足迹得分分布可分(TF 占位真实存在);轮廓图看类 (i) 的中心 m6A 谷;
比较不同样本(如 aged vs youth)类 (i) 占比变化即 TF 占位差异。

实现要点:按染色体流式读取各元素 BED,每纤维块合并为 span 后用
searchsorted/isin 向量化分类,与旧逐条扫描实现输出一致但快数百倍。

## 3. 参考文献

1. Peter CJ et al. *Single chromatin fiber profiling and nucleosome position
   mapping in the human brain*. Cell Reports Methods 4(12):100911, 2024.
   DOI: 10.1016/j.crmeth.2024.100911 —— 三分类足迹方法出处 [1]
2. fibertools-rs — https://github.com/fiberseq/fibertools-rs
