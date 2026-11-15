# arriba/polygenomes — 多基因组融合基因检测子模块(arriba)

按基因组分层的 Arriba 融合基因检测:每 (基因组, 样本) 从 STAR 嵌合比对
BAM 中调用融合转录本,再按基因组汇总为 HTML 报告 + 复发/高置信/in-frame
融合表 + 6 张汇总图。C 形态子模块,共享父目录 `../arriba.yaml` 环境与
父目录 `modules/arriba/bin/summarize_arriba_fusions.py` 汇总脚本。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_substring` | 形态 A 拼接:`{indir}/{genome}/{sample_id}/{sample_id}[.{bam_substring}].bam`(STAR 嵌合比对 BAM;`bam_substring` 默认 `""` 即裸 `{sample_id}.bam`) |
| `genome_samples` | genome → 样本列表字典,决定每基因组汇总哪些样本 |
| `genome.references.<g>.fasta` / `.gtf` | 每基因组参考 FASTA(`-a`)与注释 GTF(`-g`);缺失报错。.smk 读 `genome.references`(复数),json 示例平铺 `genome.<g>.*`,以 .smk 为准 |
| `Params.arriba.blacklist` | Arriba blacklist 文件(`-b`),必需且校验存在 |
| `Params.arriba.known_fusions` | 已知融合文件(`-k`),必需且校验存在 |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}_passed_fusions.tsv       # 通过过滤的融合
{outdir}/{genome}/{sample_id}/{sample_id}_discarded_fusions.tsv    # 被过滤的融合
{outdir}/{genome}/arriba_report/arriba_fusion_report.html          # HTML 报告
{outdir}/{genome}/arriba_report/per_sample_summary.tsv             # 每样本汇总
{outdir}/{genome}/arriba_report/recurrent_fusions.tsv               # 复发融合
{outdir}/{genome}/arriba_report/high_medium_confidence_fusions.tsv  # 高/中置信融合
{outdir}/{genome}/arriba_report/inframe_fusions.tsv                 # in-frame 融合
{outdir}/{genome}/arriba_report/figures/fig1_per_sample_counts.png
{outdir}/{genome}/arriba_report/figures/fig2_fusion_type_distribution.png
{outdir}/{genome}/arriba_report/figures/fig3_type_by_sample.png
{outdir}/{genome}/arriba_report/figures/fig4_reading_frame.png
{outdir}/{genome}/arriba_report/figures/fig5_recurrent_heatmap.png
{outdir}/{genome}/arriba_report/figures/fig6_inframe_support.png
```

日志:`{logdir}/{sample_id}/{genome}/arriba.log`、
`{logdir_combine}/arriba/{genome}/arriba_report.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `data/fastq` / `output` / `logs` | 输入/输出/日志目录 |
| `logdir_combine` | — | `logs` | 汇总规则日志目录 |
| `bam_substring` | 文件名插入段 | `""` | 见输入 |
| `genome_samples` | — | `{}` | genome → 样本列表 |
| `Procedure.arriba` | `arriba -x ... -o ... -O ... -a ... -g ... -b ... -k ...` | `arriba` | Arriba 可执行路径 |
| `Params.arriba.E` | `-E` | `0.3` | 固定传入;语义见 Arriba 官方文档 |
| `Params.arriba.t` | `-t` | `None` | 非空才传 |
| `Params.arriba.d` | `-d` | `None` | 非空才传 |
| `Params.arriba.p` | `-p` | `None` | 非空才传 |
| `Params.arriba.formats` | `--format`(单次,空格分隔多值) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集);`figures/<base>.<fmt>` 每种格式一个文件 |
| `ROOT_DIR` | — | `.` | 拼父目录汇总脚本 `modules/arriba/bin/summarize_arriba_fusions.py` |
| 线程/内存 | — | `threads=4`,`mem_mb=8000` | `rule_res(...)` 可覆盖 |

环境:`conda: "../arriba.yaml"`。

直接调用示例(等价引擎命令 + 父目录脚本):

```bash
arriba -x s1.bam -o s1_passed_fusions.tsv -O s1_discarded_fusions.tsv \
  -a genome.fa -g annotation.gtf -b blacklist.tsv -k known_fusions.tsv -E 0.3
python modules/arriba/bin/summarize_arriba_fusions.py \
  -p s1_passed_fusions.tsv,s2_passed_fusions.tsv \
  -d s1_discarded_fusions.tsv,s2_discarded_fusions.tsv \
  -o out/GRCh38/arriba_report
```

## 2. 算法原理与生物学意义

纯工具封装:Arriba 在 STAR 嵌合比对(需 `chimSegmentMin > 0`)基础上,用
过滤器排除常见假阳性(同源序列、read-through、随机嵌合、低复杂度),结合
blacklist 与已知融合库给出通过/丢弃两档融合调用;每条融合含基因对、断点、
支持读段、读框(in-frame)与置信度。融合基因是肿瘤驱动事件与部分发育/衰老
过程的关键信号。结果怎么读:`*_passed_fusions.tsv` 为核心结果(关注
`confidence` 高/中、`reading_frame` in-frame、`#split_reads`/#spanning_reads
支持数);`recurrent_fusions.tsv` 找多样本复发事件;fig4/fig6 汇总读框与
in-frame 支持,报告 HTML 供人工审读。

## 3. 参考文献

1. Arriba — https://github.com/suhrig/arriba
2. Uhrig S et al. *Accurate and efficient detection of gene fusions from
   RNA-seq data*. Genome Research 31(3):448–460, 2021.
   DOI: 10.1101/gr.257246.119
