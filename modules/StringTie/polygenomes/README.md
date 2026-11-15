# StringTie/polygenomes — 多基因组转录本组装与 TE 嵌合转录本子模块

按基因组分层的 StringTie 转录本组装:每 (基因组, 样本) 用参考 GTF 引导组装,
同基因组样本合并为非冗余转录本集;另带 TE 嵌合转录本检测(转录本与 TE 注释
交联)与分组汇总绘图。C 形态子模块,共享父目录 `../StringTie.yaml` 环境与
父目录 `modules/StringTie/bin/` 脚本。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 布局参数 | 形态 A 拼接,多基因组布局:`{indir}/{genome}/{sample_id}/{sample_id}.bam`(坐标排序 BAM) |
| `genome_samples` | genome → 样本列表字典,决定每基因组参与组装/汇总的样本(注意:json 示例未列出该键,以 .smk 为准) |
| `genome.references.<g>.gtf` | 每基因组参考注释 GTF(`-G`);缺失报错。.smk 读 `genome.references`(复数),json 示例为 `genome.<g>.gtf` 平铺,以 .smk 为准 |
| `genome.references.<g>.TE_gtf` | 每基因组 TE 注释 GTF,TE 嵌合规则必需 |
| `sample_groups` | 分组 → 样本列表字典,写入 `sample_groups.tsv` 供分组汇总 |

### 1.2 输出

```
{outdir}/{genome}/raw/{sample_id}/{sample_id}.gtf                          # 样本级组装
{outdir}/{genome}/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt  # TE 嵌合转录本表
{outdir}/{genome}/TE_chimeric/TE_chimeric_group_stacked.png                # 分组堆叠图
{outdir}/{genome}/TE_chimeric/TE_chimeric_te_type_top.png                  # TE 类型 top 图
{outdir}/{genome}/TE_chimeric/TE_chimeric_te_type_by_group.png             # TE 类型×分组图
{outdir}/{genome}/TE_chimeric/TE_chimeric_sample_summary.tsv               # 样本汇总
{outdir}/{genome}/TE_chimeric/TE_chimeric_group_summary.tsv                # 分组汇总
{outdir}/{genome}/TE_chimeric/TE_chimeric_te_type_counts.tsv               # TE 类型计数
{outdir}/{genome}/stringtie_merged.gtf                                     # 合并转录本集
```

中间件:`{outdir}/{genome}/TE_chimeric/sample_groups.tsv`(模块生成的分组表)。
日志:`{logdir}/{sample_id}/{genome}/stringTie.log`、`TEChimericTranscripts.log`、
`{logdir_combine}/stringtie/{genome}/TEChimericPlot.log`、`stringTieMerge.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志目录 |
| `logdir_combine` | — | `log` | 汇总规则日志目录 |
| `genome_samples` | — | `{}` | genome → 样本列表 |
| `sample_groups` | — | `{}` | 分组 → 样本列表(绘图分组) |
| `Params.StringTie.formats` | `TEChimericPlot.py --format`(可重复) | `["png"]` | TE 嵌合图输出格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集);`TE_chimeric_<base>.<fmt>` 每种格式一个文件 |
| `Procedure.stringtie` | `stringtie -o ... -G ... -p` / `stringtie --merge` | `stringtie` | StringTie 可执行路径 |
| `ROOT_DIR` | — | `.` | 仓库根,用于拼 bin/ 脚本绝对路径 |
| 线程/内存 | `-p` | stringTie: `threads=5, mem_mb=16000`;merge: `threads=1, mem_mb=32000` | `rule_res(...)` 可覆盖 |

直接调用示例(等价引擎命令 + 父目录脚本):

```bash
stringtie -o out/GRCh38/raw/s1/s1.gtf in/GRCh38/s1/s1.bam -G annotation.gtf -p 5
python modules/StringTie/bin/TEChimericTranscripts.py \
  -s out/GRCh38/raw/s1/s1.gtf -t TE.gtf -o out/GRCh38/raw/s1/s1_TE_chimeric_transcripts.txt
python modules/StringTie/bin/TEChimericPlot.py \
  -i out/GRCh38/raw -g out/GRCh38/TE_chimeric/sample_groups.tsv \
  -o out/GRCh38/TE_chimeric/TE_chimeric
stringtie --merge s1.gtf s2.gtf -o out/GRCh38/stringtie_merged.gtf -G annotation.gtf
```

## 2. 算法原理与生物学意义

- **StringTie 组装**:`-G` 注释引导模式下,用读段覆盖图(流网络)拆分出转录本
  流量并估计表达量,输出样本级 GTF;`--merge` 把多样本 GTF 合并为非冗余
  转录本集,消除样本间碎片/冗余转录本,是 DESeq2 等下游定量的标准前置。
- **TE 嵌合转录本**(模块私有脚本 `TEChimericTranscripts.py`):把组装转录本
  与 TE 注释求交,找出与 TE 元件交联的嵌合转录本(如 TE 启动子驱动、TE-基因
  融合),逐样本输出候选表;`TEChimericPlot.py` 按 `sample_groups.tsv` 分组
  汇总三类图与 TSV 计数表。
- 结果怎么读:嵌合表逐行列出转录本与 TE 的交联关系;`*_summary.tsv` 给每
  样本/分组的嵌合计数,`te_type_counts.tsv` 给 TE 类型(LINE/SINE/LTR/DNA…)
  分布——TE 嵌合转录本在衰老、肿瘤与应激响应中常被激活,是本模块关注的
  生物学信号。

## 3. 参考文献

1. StringTie — https://github.com/gpertea/stringtie ,
   https://ccb.jhu.edu/software/stringtie/
2. Pertea M et al. *Transcript-level estimation of RNA-seq reveals
   isoform-level expression in the human brain*. Nature Methods
   12(10):975–977, 2015(原文注:DOI 10.1038/nmeth.3563;若需引用请以
   官方页面为准);Pertea G, Pertea M. *GFF Utilities: GffRead and
   GffCompare*. F1000Research 9:304, 2020.
