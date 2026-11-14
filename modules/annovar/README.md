# annovar — RNA-SNP 变异注释模块(ANNOVAR)

围绕 ANNOVAR 的 RNA-SNP 分析链:TE 表达过滤 → 公共表达 BED → VCF 区间
求交 → convert2annovar → table_annovar refGene 注释(6 条规则)。

## 1. 模块输入输出与参数

### 1.1 输入

各规则串联消费上游产物(路径均相对本模块 `outdir`,由 subworkflow 对齐):

| 输入 | 规则 | 说明 |
|---|---|---|
| `counts/humanTEcount.cntTable` | TEcoutCPM | TE 计数表(上游 TEcount 产物,不在本模块生成) |
| `filter/vcf/human/{sample_id}.vcf.gz` | vcfIntersectBed | 样本变异 VCF(gzip;上游过滤产物) |
| `bed.human.gtf` / `bed.human.TE_gtf` | getBed | 基因 GTF / TE 注释 GTF(config 传路径) |

### 1.2 输出

```
{outdir}/counts/humanTEcountCPM.cntTable            # CPM 归一化计数表
{outdir}/counts/humanTEcountCommon.cntTable         # 公共表达过滤后计数表
{outdir}/counts/humanTEcountCommon.bed              # 公共表达元件 BED
{outdir}/filter/vcf/human/{sample_id}Common.vcf     # 与 BED 求交的 VCF(bedtools -wa -wb 输出)
{outdir}/annovar/human/{sample_id}/{sample_id}.avinput          # ANNOVAR 输入格式
{outdir}/annovar/human/{sample_id}/{sample_id}.GRCh38_multianno.csv  # 最终注释表(CSV)
{outdir}/TEcoutCPM.sh 等                            # 各规则回执脚本(单条引擎命令)
{outdir}/log/human/...                              # 各规则日志(log 直接落在 outdir/log 下)
```

注意:最终输出文件名固定含 `GRCh38_multianno.csv`,而 table_annovar 按
`{out}.{buildver}_multianno.csv` 命名——`buildver` 需与文件名一致(如
GRCh38),否则声明输出与实际产物名不符。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.Rscript` | 可执行文件 | 无回退(json:`/usr/bin/Rscript`) | Rscript 路径(直接 `config["Procedure"]["Rscript"]`,必填) |
| `Procedure.perl` | 可执行文件 | 无回退(json:`/usr/bin/perl`) | perl 路径(必填) |
| `Procedure.table_annovar` | table_annovar.pl | `table_annovar.pl`(json:`/opt/annovar/table_annovar.pl`) | ANNOVAR 表格注释脚本 |
| `bed.human.gtf` | `--Gtf` | — | 基因 GTF(getBed) |
| `bed.human.TE_gtf` | `--TEGtf` | — | TE 注释 GTF(getBed) |
| `annovar.human.db` | table_annovar 第 2 位置参数 | — | ANNOVAR 数据库目录 |
| `annovar.human.buildver` | `-buildver` | — | 基因组版本(见上,输出名固定 GRCh38) |
| `conda.RNA-SNP` | — | null(回退 `annovar.yaml`) | commonExpression/getBed/vcfIntersectBed 的 conda 环境覆盖 |
| `outdir` | — | `output` | 输出目录(所有中间产物相对它) |
| `logdir` | — | `log` | 读取但日志实际写在 `{outdir}/log/...` |

固定写死的引擎选项:commonExpression `--threshold 5`;annovar_convert
`-format vcf4 -withfreq`;vcfIntersectBed `bedtools intersect -a -b -wa -wb`;
annovar_table `-remove -protocol refGene -operation g -nastring . -csvout`
(convert2annovar.pl 固定 `/opt/annovar/convert2annovar.pl`)。

注意:`annovar.json` 中的 `genome.gtf` / `genome.TE_gtf` / `Params.annovar.db`
/ `Params.annovar.buildver` 与 .smk 实际读取的 `bed.human.*` / `annovar.human.*`
键不一致,接线以 .smk 为准;.smk 引用的 `scripts/SNP/run-NormCountMat.R`、
`commonExpression.py`、`getBed.py` 当前仓库中不存在,启用前需补齐。

线程/内存经 `rule_res(...)` 配置:vcfIntersectBed / annovar_convert /
annovar_table threads 默认 4;mem_mb 默认 8000(annovar_table 16000)。

引擎命令(以 annovar_table 为例,回执脚本原文):

```bash
perl /opt/annovar/table_annovar.pl {s}.avinput {db} -buildver {buildver} \
    -out {outdir}/annovar/human/{s}/{s} -remove -protocol refGene \
    -operation g -nastring . -csvout
```

## 2. 算法原理与生物学意义

模块把变异注释限定在"有表达证据"的元件上:先对 TE 计数表做 CPM 归一
(消除文库深度差异),按表达阈值筛出在多样本中共同表达的元件
(`--threshold 5`),转成 BED;再用 bedtools intersect 把样本 VCF 中落在
这些元件内的变异切出;最后 ANNOVAR convert2annovar(`-withfreq` 保留频率
列)转格式,table_annovar 以 refGene 协议注释变异的基因功能影响(外显子/
内含子/剪接位点等)并输出 CSV。生物学意义:RNA 来源的变异需结合表达
上下文解读——只保留表达元件内的位点能压缩假阳性,refGene 功能注释给出
编码影响,用于 RNA-SNP/RNA 编辑候选位点筛选。

## 3. 参考文献

1. ANNOVAR — https://annovar.openbioinformatics.org ;Wang K, Li M,
   Hakonarson H. *ANNOVAR: functional annotation of genetic variants from
   high-throughput sequencing data*. Nucleic Acids Res 38(16):e164, 2010.
   DOI: 10.1093/nar/gkq603
2. bedtools — https://bedtools.readthedocs.io/ ;Quinlan AR, Hall IM.
   *BEDTools: a flexible suite of utilities for comparing genomic features*.
   Bioinformatics 26(6):841–842, 2010. DOI: 10.1093/bioinformatics/btq033
