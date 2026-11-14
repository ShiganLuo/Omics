# TEtranscripts/polygenomes — 多基因组 TE 表达定量子模块

C 形态子模块(TEcount / TElocal 定量 + combine 汇总):按 genome × sample
定量转座元件(TE)表达,再按 genome 合并为表达矩阵。conda 环境共享父
目录 `../TEtranscripts.yaml`;合并脚本位于共享父目录
`modules/TEtranscripts/bin/`(`combineTE.py`、`geneId2Name.py`,经
`ROOT_DIR` 引用)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{indir}/{genome}/{sample_id}/{sample_id}.bam` | 比对 BAM(按 genome 分层;TEcount / TElocal 共用) |
| `genome.references.{genome}.gtf` | 基因注释 GTF(存在性校验,缺失抛 `ValueError`) |
| `genome.references.{genome}.TE_gtf` | TE 注释 GTF(TEcount `--TE`) |
| `genome.references.{genome}.TEind` | TE 位点索引 `TE.gtf.locInd`(TElocal `--TE`) |
| `genome_samples`(config) | genome → 样本列表(combine 规则按它聚合 cntTable,空则抛 `ValueError`) |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.TEcount.cntTable   # 单样本 TE 家族级计数
{outdir}/{genome}/{sample_id}/{sample_id}.TElocal.cntTable   # 单样本 TE 位点级计数
{outdir}/{genome}/{genome}_TEcount.tsv                       # 合并计数矩阵
{outdir}/{genome}/{genome}_TEcount_name.tsv                  # 基因 ID → 名称注释版
{outdir}/{genome}/{genome}_TElocal.tsv                       # 位点级合并矩阵
{outdir}/{genome}/{genome}_TElocal_name.tsv                  # 注释版
logdir/{sample_id}/{genome}/TEcount.log                      # 定量日志(TElocal.log 同目录)
logdir_combine/TEtranscripts/{genome}/combine_TEcount.log    # 合并日志(默认 log/combine/...)
```

`TEtranscripts_result` 无 output,消费两个合并矩阵作为 DAG 终点锚点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.TEcount` | 可执行文件 | `TEcount` | TEcount 路径 |
| `Procedure.TElocal` | 可执行文件 | `TElocal` | TElocal 路径 |
| `genome.references.{genome}.gtf` | `--GTF` | — | 基因 GTF(存在性强制) |
| `genome.references.{genome}.TE_gtf` | `--TE` | — | TE 注释 GTF(存在性强制) |
| `genome.references.{genome}.TEind` | `--TE` | — | TE 位点索引(存在性强制) |
| `genome_samples` | — | {} | genome → 样本列表(combine 输入) |
| `ROOT_DIR` | — | `./` | 仓库根,拼 `modules/TEtranscripts/bin/*.py` 路径 |
| `logdir_combine` | — | `log/combine` | combine 规则日志目录 |
| `indir` / `outdir` / `logdir` | — | `output/raw_fastq` / `output` / `log` | 目录基路径 |

注意:`TEtranscripts.json` 的参考键为 `genome.{genome}.TE_gtf/gtf/TEind`,
而 .smk 实际读取 `genome.references.{genome}.*`,接线以 .smk 为准。

线程/内存经 `rule_res(...)` 配置:TEcount/TElocal threads 默认 2、mem
32000;combine 规则 threads 默认 1、mem 8000。

引擎命令(回执脚本原文):

```bash
TEcount --sortByPos --format BAM --mode multi -b {s}.bam --GTF {gtf} --TE {TE_gtf} --project {s}.TEcount --outdir {outdir}/{g}/{s}
TElocal --sortByPos -b {s}.bam --GTF {gtf} --TE {TEind} --project {s}.TElocal   # 产物 mv 到输出名
python {ROOT_DIR}/modules/TEtranscripts/bin/combineTE.py -p TEcount -i {outdir}/{g} -o {g}_TEcount.tsv
python {ROOT_DIR}/modules/TEtranscripts/bin/geneId2Name.py -c {g}_TEcount.tsv -g {gtf} -o {g}_TEcount_name.tsv
```

## 2. 算法原理与生物学意义

TEcount(TEtranscripts 套件)把 reads 同时对基因模型与 TE 注释计数,
`--mode multi` 把多重比对 reads 按可能性分摊到基因/TE,`--sortByPos`
要求坐标排序 BAM;输出家族级 cntTable。TElocal 用 `TE.gtf.locInd` 位点
索引把 TE reads 落到具体基因组位点,输出位点级计数。combineTE.py 把各样
本 cntTable 拼成矩阵,geneId2Name.py 依 GTF 补基因名。生物学意义:TE
在发育、疾病与应激中可被重新激活,表达量变化反映基因组稳定性与调控
重塑;家族级与位点级两份矩阵分别支持"哪些 TE 家族变了"与"哪些位点
的 TE 被转录"两个层面的解读。

## 3. 参考文献

1. TE Toolkit(TEtranscripts / TEcount / TElocal)—
   https://github.com/mhammell-laboratory/TEtranscripts
2. Jin Y, Tam OH, Paniagua E, Hammell M. *TEtranscripts: a package for
   including transposable elements in differential expression analysis of
   RNA-seq datasets*. Bioinformatics 31(22):3593–3599, 2015.
   DOI: 10.1093/bioinformatics/btv422
