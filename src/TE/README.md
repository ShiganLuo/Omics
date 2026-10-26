# TE 绘图脚本（bulk 转录组，展示 TE 大爆发）

三个独立 CLI 脚本，"全局—家族—位点"三层展示 TE 爆发：

| 脚本 | 图 | 输入数据 |
|------|----|----------|
| `te_manhattan.py` | 全基因组 Manhattan 散点（爆发的全局性） | 位点级数值表（坐标在表内或经 rmsk TE GTF 解析） |
| `te_family_boxplot.py` | 家族表达分布箱线/小小提琴（阶段间分布上移） | TEcount 合并计数表或 DESeq2 结果表 + sample→group 注释 |
| `te_locus_coverage.py` | IGV 式堆叠覆盖度轨道（代表位点证据） | BAM / bigWig / bedGraph 轨道 + 基因组区间 |

公共约定：NumPy 风格英文 docstring + 全量类型标注；日志走
`modules/function/bin/LogUtil.py`（`setup_logger(__name__)`）；多格输出
`-f png -f pdf`（`action="append"`）；依赖不自动安装（matplotlib 必需、
pysam/pyBigWig 按输入格式惰性导入，缺失时报错提示手动安装）。
建议用 `DNA` 环境运行（pandas/matplotlib/pysam/pyBigWig 齐全）。

## 1. te_manhattan.py — 全基因组 TE 位点散点

x 轴按染色体依次排开全部 TE 位点坐标，y 轴为位点统计量
（如 DESeq2 log2FoldChange），爆发时整片点抬升。

输入（`--values` TSV，一行一个 TE 位点，带表头）：

- 坐标模式（默认）：表内含坐标列（`--chrom-col/--start-col/--end-col`，
  默认 `chrom/start/end`），可选 `name/family/class` 列
  （`--name-col` 等）。来源：对
  `<ref>_GENCODE_rmsk_TE.bed/.gtf` 做位点级计数（featureCounts
  `-g transcript_id` / TElocal）+ 差异分析。
- ID 模式（给 `--te-gtf` 时）：表内为位点 ID 列（`--locus-col`，默认
  `locus_id`），与 rmsk TE GTF 的 `transcript_id` 匹配，坐标和
  name/family/class 从 GTF 解析。

数值列 `--value-col`（默认 `value`）；显著性列 `--signif-col`（如
`padj`）与 `--value-thresh`/`--signif-thresh`/`--burst-sign` 联合定义
burst 位点（红色标注，`--label-top N` 标注最强的 N 个）；
`--highlight-family`（可重复）按 name/family/class 子串高亮家族。

注意：TEcount/TEtranscripts 输出是家族级、无坐标，不能喂给本图；
必须位点级计数。染色体顺序默认 chr1..N,X,Y,M，`--chrom-order` 可指定。

渲染参数：`--max-points`（默认 300000，仅对背景点抽样，burst/高亮点
全保留；3.8M rmsk 位点全画需 >15 分钟）、`--seed`、`--title`（默认无
标题）。染色体 x 标签按 45° 旋转所需间距自动抽稀（所有染色体保留刻度，
scaffold 等小染色体的标签可能不显示，日志会报告 kept N/M）；
`--label-top N` 标注 burst+highlight 中 |value| 最大的 N 个位点
（建议 ≤30）。`--hline` 参考线标注自动置顶。

```bash
python te_manhattan.py --values locus_deseq2.tsv --te-gtf GRCm39_GENCODE_rmsk_TE.gtf \
    --locus-col transcript_id --value-col log2FoldChange --signif-col padj \
    --value-thresh 1 --signif-thresh 0.05 \
    --highlight-family MT2_Mm --highlight-family IAP --label-top 8 \
    --ylabel "log2FC" -o outdir -f png -f pdf
```

## 2. te_family_boxplot.py — 家族表达分布（按组/阶段）

每个 `--family` 模式一个面板，面板内每个样本组一个小提琴/箱线，
展示该家族所有 TE 条目表达分布随阶段整体上移（爆发）。

输入：

- `--matrix` TSV：TE id 列 + 样本列。
  - TEcount 合并表（`<ref>_TEcount.tsv`，首列 `gene/TE`，id 形如
    `L1MD1:L1:LINE` = `name:family:class`，基因行参与 CPM 库大小但不
    参与家族匹配）；
  - DESeq2 结果表（`*.TEcount_TE.tsv`）：`baseMean/log2FoldChange/
    lfcSE/stat/pvalue/padj` 自动排除，样本列为 DESeq2 标准化计数，
    此时用 `--normalize none`。
- `--annotation` TSV：表头 `sample`/`group`（可选 `order` 列指定组序，
  `--group-order` 逗号分隔可覆盖）。

匹配：`--match-field {name,family,class,any}`（默认 any=全 id 串），
`--match-mode {substring,regex,exact}`（默认大小写不敏感子串）。
归一化 `--normalize {log2cpm,cpm,none}`（默认 log2cpm）。
分布口径 `--value-scope {group_mean_per_entry,all_samples}`（默认每条目
取组内重复均值、跨条目成分布）。图型 `--kind {box,violin}`；
`--test mannwhitney` 对每组标注 vs 参考组（`--test-ref-group`，默认首组）
的 Mann-Whitney p（需 scipy）。布局参数 `--ncols/--row-spacing/
--stat-headroom/--figsize-unit` 可调。

```bash
python te_family_boxplot.py --matrix GRCm39_TEcount.tsv --annotation samples.tsv \
    --family MT2 --family L1Md --family IAP --kind violin --ncols 2 \
    --test mannwhitney --test-ref-group mESC -o outdir -f png
```

## 3. te_locus_coverage.py — 位点覆盖度轨道（IGV 式）

每个区间一张图：样本覆盖度轨道竖直堆叠（共享 x），底部 TE/基因
注释 lane + 独立标签带，顶部组图例（按 `--annotation` 分组着色）。

输入：

- `--track NAME=PATH`（可重复）：按扩展名自动选择读取器 —
  `.bam/.cram`（pysam，需索引；`--normalize cpm` 用 1e6/mapped 校正）、
  `.bw/.bigWig`（pyBigWig，视为已归一化）、
  `.bedGraph/.bg[.gz]`（4 列区间，视为已归一化）。
- `--region "[NAME=]chr:start-end"`（可重复，1-based 闭区间，同 IGV）或
  `--regions-file` TSV（表头 `name/chrom/start/end`）。
- `--features-gtf` / `--features-bed`（二选一）：rmsk TE GTF
  （`gene_id/family_id/class_id` 属性）或 7 列 BED
  （`chrom start end strand gene_id family_id class_id`）。
- `--annotation` TSV：`sample/group`（可选 `order`），轨道填色按组 + 图例。

y 轴 `--share-y {max,track}`（默认 max 全轨道共享上限；track 则各自
缩放并在轨道标签里注明 max）。`--bin-size`（默认 50 bp，<=1 为逐碱基）、
`--flank` 侧翼、`-j` 多区间并行（ProcessPoolExecutor）。
布局参数 `--feature-lane-height/--label-band-unit/--label-fontsize/
--legend-row-height/--figsize` 可调。

```bash
python te_locus_coverage.py \
    --track mESC1=mESC1.bam --track TLSC1=TLSC1.bam \
    --annotation samples.tsv \
    --region MT2_Mm_chr1=chr1:3146900-3147900 \
    --features-gtf GRCm39_GENCODE_rmsk_TE.gtf --bin-size 25 \
    -o outdir -f png
```

## 参考数据路径（本项目）

- TEcount 合并计数：`output/pluripotent2totipotent/RNAseq/counts/TEtranscripts/<ref>/<ref>_TEcount.tsv`
- DESeq2 TE 结果：`output/pluripotent2totipotent/RNAseq/diff_expression/<ref>/<comp>/<comp>.TEcount_TE.tsv`
- RNA BAM：`output/pluripotent2totipotent/RNAseq/common/3_TEtranscripts_bam/<ref>/<sample>/<sample>.bam`
- rmsk TE 注释：
  - `/data/pub/zhousha/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf`（+ 同名 `.bed`）
  - `/data/pub/zhousha/Database/Reference/human/GENCODE/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf`
