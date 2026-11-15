# arriba — RNA-seq 融合基因检测模块(Arriba)

Arriba 融合基因检测封装:输入 STAR 比对 BAM,输出通过/被滤除的融合
候选表;自研脚本 `summarize_arriba_fusions.py` 汇总多样本结果出 HTML
报告、复现融合表与 6 张统计图。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | BAM 目录(上游 config 的 `outdir`),形态 A 布局:`{indir}/{sample_id}/{sample_id}{.bam_substring}.bam`;须为 STAR 比对产物(`-x` 主比对文件,`--chimOutType WithinBAM` 风格) |
| `bam_substring` (config,顶层) | BAM 文件名后缀段;缺省空 |
| `samples` (config) | 样本 ID 列表(报告规则 expand 用) |
| `genome.fasta` / `genome.gtf` (config) | 参考基因组 FASTA(`-a`)与基因注释 GTF(`-g`) |
| `Params.arriba.blacklist` (config) | 黑名单 TSV(`-b`);**文件不存在直接报错**(运行期校验) |
| `Params.arriba.known_fusions` (config) | 已知/复现融合 TSV(`-k`);**文件不存在直接报错** |
| `ROOT_DIR` (config) | 项目根目录,定位 `modules/arriba/bin/summarize_arriba_fusions.py` |

### 1.2 输出

```
# arriba(每样本)
{outdir}/{sample_id}/{sample_id}_passed_fusions.tsv       # 通过全部过滤器的融合
{outdir}/{sample_id}/{sample_id}_discarded_fusions.tsv    # 被滤除的候选
{outdir}/{sample_id}/{sample_id}_arriba.sh                # 回执

# arriba_report(多样本汇总)
{outdir}/arriba_report/arriba_fusion_report.html          # HTML 报告
{outdir}/arriba_report/per_sample_summary.tsv             # 每样本融合计数汇总
{outdir}/arriba_report/recurrent_fusions.tsv              # 复现(多样本共有)融合
{outdir}/arriba_report/high_medium_confidence_fusions.tsv # 高/中置信融合
{outdir}/arriba_report/inframe_fusions.tsv                # in-frame 融合
{outdir}/arriba_report/figures/fig1_per_sample_counts.png # 每样本融合计数
{outdir}/arriba_report/figures/fig2_fusion_type_distribution.png
{outdir}/arriba_report/figures/fig3_type_by_sample.png
{outdir}/arriba_report/figures/fig4_reading_frame.png
{outdir}/arriba_report/figures/fig5_recurrent_heatmap.png
{outdir}/arriba_report/figures/fig6_inframe_support.png
```

日志:`{logdir}/{sample_id}/arriba.log`、`{logdir}/arriba/arriba_report.log`。

### 1.3 参数

config 位于 `Params.arriba` 与 `Procedure` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.arriba.blacklist` | `-b` | — | 黑名单(artifacts 区段);必填且须存在 |
| `Params.arriba.known_fusions` | `-k` | — | 已知融合列表(提升灵敏度);必填且须存在 |
| `Params.arriba.t` | `-t` | null | 融合 tags TSV(前两列基因、第三列 tag,写入输出 `tags` 列) |
| `Params.arriba.d` | `-d` | null | WGS 结构变异坐标 TSV(佐证弱表达融合、消低置信) |
| `Params.arriba.E` | `-E` | 0.3 | `relative_support` 过滤器的 e-value 上限(支持读段数随机期望 > 阈值即弃;调高会显著增假阳性) |
| `Params.arriba.p` | `-p` | null | 蛋白结构域 GFF3(输出 `retained_protein_domains` 列) |
| `Params.arriba.formats` | `--format`(单次,空格分隔多值) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集);`figures/<base>.<fmt>` 每种格式一个文件 |
| `Procedure.arriba` | 引擎可执行 | `arriba` | — |
| (threads) | — | 4 | `rule_res("arriba")` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

自研汇总脚本参数(`summarize_arriba_fusions.py`,与 --indir 二选一):

| 参数 | 默认 | 说明 |
|---|---|---|
| `-p` / `--passed` | 必填(与 --indir 互斥) | passed TSV 列表(规则按逗号拼接传入) |
| `-d` / `--discarded` | null | discarded TSV 列表 |
| `-o` / `--outdir` | 必填 | 报告输出目录 |
| `-f` / `--format` | `png` | 出图格式(单次跟一个或多个值,如 `-f png pdf svg`;合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba) |
| `--indir` | — | 直接扫目录取 tsv 的替代入口 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
arriba -x sample.bam \
    -o sample_passed_fusions.tsv -O sample_discarded_fusions.tsv \
    -a genome.fa -g annotation.gtf \
    -b blacklist.tsv.gz -k known_fusions.tsv.gz \
    -E 0.3

python modules/arriba/bin/summarize_arriba_fusions.py \
    -p s1_passed_fusions.tsv,s2_passed_fusions.tsv \
    -d s1_discarded_fusions.tsv,s2_discarded_fusions.tsv \
    -o out/arriba_report
```

## 2. 算法原理与生物学意义

Arriba 从 STAR 比对中提取 split reads 与 discordant mates 候选,经
多级过滤器(黑名单、read-through、短锚、同源基因、多重比对、错配、
`relative_support` e-value 等)判别真实融合,除基因融合外还能检出病毒
整合、内部串联重复(ITD)、外显子重复、基因内倒位、增强子劫持等结构
事件。要点:

- **passed vs discarded**:两表同构,`discarded` 保留被滤候选及滤除原因,
  便于审计过滤严格度;
- **`-E`(e-value)**:按支持读段数估计随机碰撞期望,超过阈值判假
  (filter `relative_support`,默认 0.3);
- **`-d` WGS 佐证**:有基因组断点证据的弱表达融合可救回,反之降置信;
- **`-p` 蛋白结构域**:融合是否保留功能域(in-frame、retained domains)
  决定其致癌机制可能性,报告的 `inframe_fusions.tsv` 与 reading-frame
  图即此判读;
- 汇总侧:复现融合(多样本同一基因对)、高/中置信集、每样本计数是
  队列层面筛选候选驱动融合的入口。

生物学意义:融合基因是肿瘤(尤其血液与肉瘤)的关键驱动事件与用药
靶点(如 NTRK、ALK 融合);Arriba 临床研究导向的设计强调灵敏度与
假阳性控制,confidence 列(high/medium/low)指导后续验证优先级。

## 3. 参考文献

1. Arriba — https://github.com/suhrig/arriba
   (命令行手册 https://github.com/suhrig/arriba/wiki/07-Command-line-options)
2. Uhrig S, Ellermann J, Walther T, et al. *Accurate and efficient
   detection of gene fusions from RNA sequencing data*. Genome Research
   31:448-460, 2021. DOI: 10.1101/gr.257246.119.
