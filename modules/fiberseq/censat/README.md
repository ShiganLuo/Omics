# fiberseq/censat — 着丝粒区域类别可及性/甲基化谱子模块

按 CenSat 注释的着丝粒区域类别(CDR/HOR/侧翼)分层统计 Fiber-seq 的
MSP 及性与 CpG 甲基化谱,输出表格 + 谱图(formats 可配,默认 PDF)。conda 环境共享父目录
`../fiberseq.yaml`(容器同源);引擎脚本共享父目录
`modules/fiberseq/bin/censat_profile.py`(经 `config["ROOT_DIR"]` 拼接)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `genome.censat_bed`(config,必填) | CenSat 区域注释 BED(着丝粒区域类别划分);缺失时 `_require()` 报错 |
| `indir`(默认 `common/4_elements`)+ 样本 ID 布局(形态 A) | `{indir}/{sample_id}/{sample_id}.msp.bed.gz`、`.nuc.bed.gz`、`.cpg.bed.gz`:上游 Fiber-seq 元素提取(4_elements)产物 |
| `genome.mcpg_pileup` | 可选 mCpG pileup 文件(提供时以 `-p` 传入,增强甲基化谱) |
| `samples` | 样本 ID 列表(供 result 规则聚合) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.censat_regions.tsv        # 区域级统计
{outdir}/{sample_id}/{sample_id}.censat_class_summary.tsv   # 类别汇总(CDR/HOR/flank)
{outdir}/{sample_id}/{sample_id}.censat_profile.tsv          # 分箱谱数据
{outdir}/{sample_id}/{sample_id}.censat_profile.<fmt>       # 谱图(formats 多格式,默认 pdf)
{logdir}/{sample_id}/fiberseq_censat_run.log
```

收尾规则 `fiberseq_censat_result` 聚合各样本 `censat_class_summary.tsv`。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.censat_bed` | `-r` | — | CenSat 注释 BED,必填 |
| `genome.mcpg_pileup` | `-p` | 不传 | 可选 mCpG pileup |
| `Params.fiberseq_censat.msp_min_len` | `-m` | 50 | 计入的长 MSP 最小长度(bp) |
| `Params.fiberseq_censat.flank` | `-f` | 10000 | 区域侧翼延伸(bp) |
| `Params.fiberseq_censat.bin` | `-b` | 500 | 谱图分箱宽度(bp) |
| `Params.fiberseq_censat.formats` | `--format`(可重复) | ["pdf"] | 图片输出格式(多种同时输出),合法值 png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba |
| `indir` / `outdir` / `logdir` | — | `common/4_elements` / `fiberseq_analysis/censat` / `log` | 输入/输出/日志根目录 |

资源默认:`threads = 4`、`mem_mb = 8000`(经 `rule_res("fiberseq_censat")`
可覆盖)。

直接调用示例(共享引擎脚本 `modules/fiberseq/bin/censat_profile.py`):

```bash
python $ROOT_DIR/modules/fiberseq/bin/censat_profile.py \
  -r censat.bed -s s1.msp.bed.gz -u s1.nuc.bed.gz -c s1.cpg.bed.gz \
  -o output/censat --sample s1 -m 50 -f 10000 -b 500 -l log.txt [-p mcpg_pileup]
```

## 2. 算法原理与生物学意义

人类着丝粒由高阶重复(HOR,alpha-satellite 大阵列)与散布的二倍体样区域
(CDR)构成,两侧为侧翼序列;CenSat 注释按序列/表观特征给区域分类。本模块
按类别分层统计 Fiber-seq 可及性指标:以长度 ≥ `msp_min_len` 的长 MSP
(甲基化敏感区 = 无核小体足迹的开放染色质)为单位,计算区域核心上的池化
覆盖占比与逐纤维平均覆盖占比(裁剪对齐区域核心,逐纤维覆盖用区间并集计算,
内存 O(核心区间)),按 `bin` 宽度分箱输出 CDR/HOR/侧翼的可及性与 CpG
甲基化 profile。生物学意义:活跃着丝粒(含 CENP-A 的 HOR 阵列)与失活
着丝粒/CDR 的染色质开放度与甲基化不同,分层 profile 可比较不同样本(如
aged/youth)着丝粒染色质组织差异。结果阅读:`class_summary.tsv` 看各类别
汇总指标;`profile.tsv/pdf` 看沿区域的分箱曲线。

## 3. 参考文献

1. 本模块引擎脚本:`modules/fiberseq/bin/censat_profile.py`;
   Fiber-seq 元素(MSP/核小体/CpG)定义见 `fiberseq/` 父目录
   (`FIBERSEQ_ANALYSIS.md` / fibertools-rs)
2. CenSat 区域注释:Altemose N, et al. *Complete genomic and epigenetic
   maps of human centromeres*. Science 2022, 376(6588):eabl4178.
   DOI: 10.1126/science.abl4178
