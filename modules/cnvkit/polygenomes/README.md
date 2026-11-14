# cnvkit/polygenomes — 多基因组 CNV 拷贝数变异检测模块(CNVkit)

CNVkit 的 polygenomes 子模块:按基因组(`{genome}` 通配符)分别建参考、
跑 batch 分析,输出每样本 log2 拷贝数比值(`.cnr`)与分段(`.cns`)。
多基因组布局:基因组名进输出路径段,`config["genome"][<g>]` 提供该基因组
的 fasta 与 access BED。子模块共享父目录 conda 环境(`conda: "../cnvkit.yaml"`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | BAM 目录;`{indir}/{sample_id}/{sample_id}.bam`(形态 A 布局拼接) |
| `samples` / `control_samples`(config) | 待测样本 / 对照样本 ID 列表;`control_samples` 非空才建参考(`cnvkit_reference`) |
| `genome.<g>.fasta` | 每个基因组的参考 FASTA(`<g>` = `{genome}` 通配符取值) |
| `genome.<g>.access` | 每个基因组的 accessible 区域 BED(可比对区域) |
| `reference.cnn` | 上游 `cnvkit_reference` 产出,`control_samples` 非空时作为 `cnvkit_batch` 输入 |

### 1.2 输出

```
{outdir}/{genome}/reference/reference.cnn        # 对照样本构建的参考(有 control_samples 时)
{outdir}/{genome}/cnv/{sample_id}.cnr            # 每样本 log2 拷贝数比值(bin 级)
{outdir}/{genome}/cnv/{sample_id}.cns            # 每样本分段拷贝数(segment 级)
```

`cnvkit_batch` 同时产出 scatter/diagram 图(未声明为输出);收尾规则
`cnvkit_result` 聚合全部 `.cnr` + `.cns`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.cnvkit.method` | `--method` | `wgs` | CNVkit 取样方法,枚举 `hybrid`/`amplicon`/`wgs`(schema 声明必填;无对照走 fasta/access 路径时使用) |
| `Procedure.cnvkit` | — | `cnvkit.py` | CNVkit 可执行路径 |
| — | `--processes` | threads | 并行进程数(默认 8,`rule_res` 配置) |
| — | `--scatter --diagram` | 固定 | batch 同时出散点图与概览图 |
| threads | — | 8(reference/batch) | `rule_res` 配置 |
| mem_mb | — | reference=16000,batch=32000 | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
# 1) 有对照样本:建参考
cnvkit.py batch control1.bam control2.bam \
  --method wgs --fasta genome1.fa --access genome1.access.bed \
  --processes 8 --output-reference reference.cnn

# 2) batch 分析(有参考)
cnvkit.py batch S1.bam S2.bam --reference reference.cnn \
  --processes 8 --output-dir out/cnv --scatter --diagram

# 3) batch 分析(无参考,直接用 fasta/access)
cnvkit.py batch S1.bam S2.bam \
  --method wgs --fasta genome1.fa --access genome1.access.bed \
  --processes 8 --output-dir out/cnv --scatter --diagram
```

## 2. 算法原理与生物学意义

CNVkit 把测序 reads 在基因组 bin 中的覆盖深度对 GC 含量、可比对性等
系统偏差做校正,以对照样本(或 flat reference)归一化得 log2 拷贝数比值
(`.cnr`),再经 CBS 等分段算法聚合成拷贝数片段(`.cns`),可检出缺失/
扩增/杂合性缺失。生物学意义:体细胞拷贝数变异是肿瘤基因组的标志性
改变(polygenomes 场景下按各参考基因组分别评估);`method` 与建库类型
匹配(hybrid capture / amplicon / WGS)才能正确建模覆盖度偏差。
结果读法:`.cnr` 行为 bin、`log2` 列为相对拷贝数(0 ≈ 二倍体,负为缺失,
正为扩增);`.cns` 行为片段、`log2` 为片段均值,`cn` 为绝对拷贝数估计。

## 3. 参考文献

1. CNVkit — https://cnvkit.readthedocs.io/(官方文档);
   https://github.com/etal/cnvkit
2. Talevich E, Shain AH, Botton T, Bastian BC. *CNVkit: Genome-Wide Copy
   Number Detection and Visualization from Targeted DNA Sequencing*.
   PLoS Computational Biology 12(4):e1004873, 2016.
   DOI: 10.1371/journal.pcbi.1004873
