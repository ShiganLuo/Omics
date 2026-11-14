# fiberseq/mods — Fiber-seq 区域轴修饰长表与显著结果拆分子模块

按区域轴(FIRE peaks 或 TE 位点)聚合每纤维修饰(m6A/CpG)长表,交给
独立的 sigtest 模块做显著性检验,再后处理 sigtest 结果(显著对象拆分、
TE 轴 top-loci BED 供 regionplot 浏览)。conda 环境共享父目录
`../fiberseq.yaml`;引擎脚本在本目录 `bin/`(自包含)。与 sigtest 的跨
模块契约仅限文件(长表 → sigtest 输入;`significance_test.tsv` → 本模块
split/top_loci 输入,经 `sigtest_dir` 接线)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(默认 `common/4_elements`)per-sample | `{s}.{mod_type}.bed.gz`(默认 `m6a`)、`{s}.msp.bed.gz`;fire_peaks 轴另需 `{s}.fire_peaks.bed`(均为上游 4_elements 产物) |
| `te_bed` 或 `genome.repeat_bed`(TE 轴,必填) | RepeatMasker TE BED:第 4 列 subfamily、第 7 列 family、第 8 列 class |
| `sigtest_dir` | sigtest 模块输出根(默认 `{outdir}/sigtest`);split/top_loci 消费其中 `{contrast}/significance_test.tsv` |
| `contrasts` | 对比名列表(驱动 `{contrast}` wildcard) |
| `Params.fiberseq_mods.group_samples` | 组 → 样本列表(node.py 从 meta 设计派生,勿手填) |

### 1.2 输出

```
# longtable(fire_peaks 轴;metric 默认 m6a_per_fiber)
{outdir}/{metric}_long.tsv
# longtable(TE 轴)
{outdir}/te_{mod_type}_long.tsv
{outdir}/te_{mod_type}_loci.tsv          # 逐位点表(te_top_loci 消费)
# split(写在 sigtest 输出目录,与 significance_test.tsv 同目录)
{sigtest_dir}/{contrast}/significant_objects.tsv
{sigtest_dir}/{contrast}/significant_objects.bed   # 仅 fire_peaks 轴
{sigtest_dir}/{contrast}/nonsignificant_objects.tsv
# top_loci(仅 TE 轴)
{sigtest_dir}/{contrast}/top_loci.tsv / top_loci.bed
# 日志
{logdir}/mods/*.log
```

收尾规则 `fiberseq_mods_result` 聚合 longtable + 各 contrast 的
significant_objects.tsv(fire_peaks 轴含 BED;TE 轴含 top_loci.bed)。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.fiberseq_mods.mod_type` | `--mod-type` | `m6a` | 修饰类型(决定输入文件名与值列后缀 `{mod_type}_per_fiber`) |
| `Params.fiberseq_mods.region_axis.default` | — | `fire_peaks` | 区域轴:`fire_peaks` / `te`(= `te` 时启用 TE 分支) |
| `...region_axis.options.fire_peaks.metric` | `--metric` | `m6a_per_fiber` | 长表取值度量(兼作输出文件名) |
| `...region_axis.options.fire_peaks.top_n` | `--top-n`(split) | 0 | 显著 BED 取 p 最小前 N 个(0 = 全部) |
| `...region_axis.options.te.min_fibers` | `--min-fibers` | 3 | 位点最少纤维数 |
| `...region_axis.options.te.loci_top_n` | `--top-n`(top_loci) | 5 | 每显著组取 top 位点数 |
| `Params.fiberseq_mods.group_samples` | `--group-samples` / `--groups` | {} | 组 → 样本映射 |
| `Params.fiberseq_mods.groups_by_contrast` / `groups` | — | {} / [] | top_loci 比较组选择(contrast 映射 → 全局组 → 全部组取前二) |
| `contrasts` | — | [] | 对比列表 |
| `sigtest_dir` | — | `{outdir}/sigtest` | sigtest 输出根 |
| `te_bed` / `genome.repeat_bed` | `--te` | — | TE BED(TE 轴必填) |
| `tmp_dir` | `--tmp-dir` | `{outdir}/tmp` | TE 轴排序缓存目录 |

split 引擎另有 `--min-abs-effect`(默认 0.0,最小效应量过滤,.smk 未转发,
CLI 直调可用)。资源默认:longtable/split/top_loci 各 threads=1、
mem_mb 8000/4000/4000(经 `rule_res()` 可覆盖)。

直接调用示例(bin/ 引擎脚本):

```bash
# fire_peaks 轴长表
python bin/extract_modification_longtable.py \
  --mods "s1:common/4_elements/s1/s1.m6a.bed.gz,s2:.../s2.m6a.bed.gz" \
  --peaks "s1:.../s1.fire_peaks.bed,s2:.../s2.fire_peaks.bed" \
  --msp "s1:.../s1.msp.bed.gz,s2:.../s2.msp.bed.gz" \
  --group-samples '{"aged":["s1"],"youth":["s2"]}' \
  --metric m6a_per_fiber --out m6a_per_fiber_long.tsv -l log

# TE 轴长表 + 逐位点表
python bin/te_modification.py --mods ... --msp ... --te repeats.bed \
  --group-samples '{...}' --mod-type m6a --min-fibers 3 \
  --out-long te_m6a_long.tsv --out-loci te_m6a_loci.tsv --tmp-dir tmp -l log

# sigtest 结果拆分(contrast 维度)
python bin/split_significant.py -i significance_test.tsv \
  --out-tsv significant_objects.tsv --out-rest-tsv nonsignificant_objects.tsv \
  --out-bed significant_objects.bed --top-n 0

# TE 显著组 top 位点
python bin/te_top_loci.py --sig significance_test.tsv --loci te_m6a_loci.tsv \
  --groups '{"aged":["s1"],"youth":["s2"]}' \
  --value-suffix m6a_per_fiber --top-n 5 --min-fibers 3 \
  --out-tsv top_loci.tsv --out-bed top_loci.bed
```

## 2. 算法原理与生物学意义

Fiber-seq 的 m6A(足迹外 m6A = 调控蛋白触达位点)与 MSP(开放染色质)可
在两类区域轴上做组间比较:**fire_peaks 轴**在 FIRE 判定的调控元件上按
`metric`(默认每纤维 m6A 均值)聚合逐纤维长表;**te 轴**在 RepeatMasker
TE 位点上定量(subfamily/family/class 多层级都写入长表,sigtest 的 feature
列选择检验层级),位点需 ≥ `min_fibers` 条纤维支撑。长表交给 sigtest 模块
统计检验(文件契约);本模块随后:

- `split_significant.py` 按显著性把 `significance_test.tsv` 拆成显著/非
  显著 TSV,fire_peaks 轴再写 regionplot 可用的 BED(feature id 形如
  `chrom:start-end|name`,可解析为坐标;TE 轴 feature 为家族名无坐标,不产
  BED),`--top-n` 按 p 值截取前 N;
- `te_top_loci.py` 在显著 TE subfamily/family 内按 |组间均值差| 排序取
  top 位点(`loci_top_n`),产出 TSV + BED——这是 regionplot 浏览候选
  清单,家族级 p 值是唯一统计论断,位点只排序不检验。

生物学意义:比较不同条件下(如 aged/youth)调控元件或 TE 家族的单分子
修饰/可及性差异,定位差异调控的家族与候选位点。结果阅读:
`significant_objects.tsv` 为显著对象主表;`top_loci.bed` 直接送
regionplot 出图;`nonsignificant_objects.tsv` 留档对照。

## 3. 参考文献

1. 本模块引擎脚本:`bin/`(extract_modification_longtable.py、
   te_modification.py、split_significant.py、te_top_loci.py、_fiberlib.py);
   统计检验在仓库独立模块 sigtest(文件契约衔接)
2. FIRE(fire_peaks 来源)—— https://github.com/fiberseq/FIRE;
   RepeatMasker(TE BED 来源)—— https://www.repeatmasker.org/
