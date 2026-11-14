# sexcheck — 样本性别判别模块(DNA/RNA 双模式)

从比对 BAM 推断样本遗传性别(DNA 模式:X/Y 染色体覆盖比;RNA 模式:
XIST/Y 标记基因表达),输出 `sexcheck.tsv`。形态 B:`sample_bams` 样本 →
BAM 路径字典由 subworkflow 按上游 config `outdir` 拼接(布局解耦)。conda
环境 `sexcheck.yaml`(容器同源);引擎脚本 `bin/sexcheck.py` 随模块携带。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `sample_bams`(config) | 样本 ID → BAM 完整路径字典(形态 B;键必须是纯 sample_id) |
| `genome.sex_chroms` | 性染色体名列表 `[X名, Y名]`(单基因组接线,node.py 可自动填充 NCBI 缩写组装);缺失时脚本按 chrX/X、chrY/Y 命名嗅探 |
| `genome.gtf` | RNA 模式标记基因注释(RNA 模式消费 BAM 计数时用) |
| `genome.gtf_map` / `genome.sex_chroms_map` / `genome.samples` | 多基因组接线:`{genome: gtf}`、`{genome: [X,Y]}`、`{sample: genome}` |

### 1.2 输出

```
{outdir}/idxstats/{sample_id}/{sample_id}.idxstats.tsv   # 单样本 samtools idxstats 表
{outdir}/idxstats_all.tsv                                # 全样本合并表(中间产物)
{outdir}/sexcheck.tsv                                    # 判型结果表
{outdir}/sexcheck.done                                   # 标记文件
{logdir}/{sample_id}/sexcheck_idxstats.log / sexcheck.log
```

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.samtools` | `--samtools` | `samtools` | samtools 可执行路径 |
| `Params.sexcheck.data_type` | `--data-type` | 脚本自适应 | `DNA`(覆盖比)/`RNA`(标记基因表达)选择判型逻辑 |
| `Params.sexcheck.ratio_xy_threshold` | `--ratio-xy-threshold` | 2.0 | X/Y 归一化覆盖比判型阈值(DNA) |
| `Params.sexcheck.y_min_frac` | `--y-min-frac` | 0.01 | Y 归一化覆盖存在下限(DNA) |
| `Params.sexcheck.xist_gene` | `--xist-gene` | `xist` | XIST 基因名(大小写不敏感匹配 GTF gene_name) |
| `Params.sexcheck.y_genes` | `--y-genes` | `ddx3y,eif2s3y,kdm5d,uty,usp9y` | Y 标记基因列表(逗号串或列表均可) |
| `Params.sexcheck.y_min_rpkm` | `--y-min-rpkm` | 0.5 | Y 标记 RPKM 下限(RNA) |
| `Params.sexcheck.xist_min_rpkm` | `--xist-min-rpkm` | 0.1 | XIST RPKM 下限(RNA) |
| `sample_bams` | `--bam-map` | {} | 样本 → BAM 字典(RNA 计数用) |
| `outdir` / `logdir` / `logdir_combine` | — | `output/QC/sexcheck` / `log` / `logdir` | 输出/日志目录 |

直接调用示例(`bin/sexcheck.py`):

```bash
python bin/sexcheck.py \
  --tsv idxstats_all.tsv --outdir output/QC/sexcheck --log sexcheck.log \
  --data-type DNA --ratio-xy-threshold 2.0 --y-min-frac 0.01 \
  --x-chrom chrX --y-chrom chrY

python bin/sexcheck.py \
  --tsv idxstats_all.tsv --outdir out --data-type RNA \
  --gtf annotation.gtf \
  --bam-map '{"s1":"s1.bam"}' \
  --xist-gene xist --y-genes ddx3y,eif2s3y,kdm5d,uty,usp9y \
  --y-min-rpkm 0.5 --xist-min-rpkm 0.1
```

## 2. 算法原理与生物学意义

性别判别是样本 QC 的标配(抓取样错配、移植瘤性别不符、细胞系污染)。

- **DNA 模式**:`samtools idxstats` 得每条染色体 mapped reads,按染色体
  长度归一化(裸 read 数受染色体长度偏倚)得 X/Y 覆盖。判定规则:
  Y 覆盖 ≥ `y_min_frac` 且 X/Y 比 < `ratio_xy_threshold` → **XY**;
  Y 覆盖低于下限(Y 基本为空)→ **XX**;参考缺 Y/常染色体缺失,或比值
  落在阈值附近的模糊带(阈值两侧半带)→ **unknown**。阈值可显式给定,
  未给定时脚本按队列自身的比值分布自适应。
- **RNA 模式**:转录组上 Y 基因大多稀疏表达,覆盖比不适用;改用标记基因
  表达:XIST(雌性 X 失活程序专属表达)与 Y 连锁标记(ddx3y、eif2s3y、
  kdm5d、uty、usp9y,雄性专属)。判定:Y 标记 RPKM ≥ `y_min_rpkm` 且
  XIST RPKM < `xist_min_rpkm` → **XY**;镜像条件 → **XX**;其余(含标记
  基因/BAM 缺失)→ **unknown**。基因名大小写不敏感,人(XIST)鼠(Xist)
  通用。

结果阅读:`sexcheck.tsv` 每样本一行给 XX/XY/unknown 与依据的比值/表达
值;unknown 需人工复核(参考缺 Y 是常见原因)。

## 3. 参考文献

1. 本模块引擎脚本:`bin/sexcheck.py`(判定逻辑与阈值定义);
   `samtools idxstats` —— https://www.htslib.org/doc/samtools-idxstats.html
