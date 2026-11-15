# trgt — 串联重复(TR)基因型分型与可视化模块

封装 PacBio TRGT(Tandem Repeat Genotyping Tool):对 HiFi 长读长
BAM 做指定串联重复位点的基因型分型,输出 VCF 与跨重复区的 spanning
reads BAM,并出单个重复位点的可视化图。两规则 `trgt_genotype` /
`trgt_plot` + 收尾 `trgt_result`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 比对目录;input 函数拼 `{indir}/{sample_id}/{sample_id}.bam` + `.bam.bai`(设 `bam_substring` 时为 `{sample_id}.{bam_substring}.bam(.bai)`) |
| `genome.fasta` / `genome.fai` | 参考基因组 FASTA 与索引(`trgt_plot` 严格校验存在) |
| `genome.repeat_bed` | TRGT 重复位点目录 BED(repeat catalog;`trgt_plot` 严格校验存在) |
| `bam_substring`(config) | 可选:BAM 文件名中缀 |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/genotype/{sample_id}/{sample_id}.trgt.vcf.gz                # TR 基因型 VCF
{outdir}/genotype/{sample_id}/{sample_id}.trgt.spanning.sorted.bam   # 跨重复区 reads BAM
{outdir}/plot/{sample_id}/{sample_id}.trgt.repeat.{fmt}              # 位点可视化图(每种格式一个文件)
```

`{fmt}` 取自 `Params.trgt.formats`(默认 `["png"]`,与旧版文件名
`{sample_id}.trgt.repeat.png` 一致)。收尾规则 `trgt_result` 聚合 vcf.gz + 全部格式的图。

### 1.3 参数

config 位于 `Params.trgt` 下,与 TRGT CLI 对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `karyotype` | `--karyotype` | `XX` | 性染色体核型(XX/XY,影响性染色体位点的倍性判定) |
| `repeat_id` | `--repeat-id` | `HTT` | 出图的重复位点 ID(如 HTT、FMR1) |
| `formats` | `plot --image <base>.<fmt>` | `["png"]` | 图片格式列表,每种格式一个输出文件;合法值 `png` `svg` `pdf`(trgt 5.1.0 `--image` 扩展名工具集,解析期校验 fail fast) |
| `Procedure.trgt` | — | `trgt` | TRGT 可执行路径 |
| — | `genotype --num-threads` | threads | 分型线程(默认 4,`rule_res("trgt_genotype")` 配置) |
| — | `genotype --genome --repeats --reads --output-prefix` | 固定 | 参考/目录/BAM/输出前缀 |
| — | `plot --genome --repeats --vcf --spanning-reads --image` | 固定 | 出图参数 |
| mem_mb | — | genotype 32000,plot 16000 | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
trgt genotype --num-threads 4 \
  --genome genome.fa --repeats repeats.bed --karyotype XX \
  --reads S1.bam --output-prefix out/S1.trgt
# 产物: out/S1.trgt.vcf.gz + out/S1.trgt.spanning.sorted.bam

trgt plot --genome genome.fa --repeats repeats.bed \
  --vcf out/S1.trgt.vcf.gz --spanning-reads out/S1.trgt.spanning.sorted.bam \
  --repeat-id HTT --image S1.trgt.repeat.png
# 一次请求多格式:重复调用,--image 换扩展名(png/svg/pdf 由扩展名决定格式)
#   trgt plot ... --image S1.trgt.repeat.svg
```

## 2. 算法原理与生物学意义

TRGT 针对串联重复(短串联重复 STR、VNTR 等)这类短读长难以分型的
区域,利用 HiFi 长读长直接跨过整个重复阵列:按 repeat catalog 中每个
位点的基序结构做逐 read 的重复次数/基序构成解析,汇总为每个位点的
基因型(VCF 中给出各等位基因的重复拷贝数与基序),并输出支持该位点
的 spanning reads 供可视化。生物学意义:串联重复扩增与数十种单基因
病(HTT→亨廷顿病、FMR1→脆性 X)及肿瘤相关扩增直接相关;TR 还是
表达变异的重要来源。结果读法:VCF 按位点给出等位基因重复长度,
`repeat.png` 逐 read 展示阵列结构(重复数/基序/甲基化),一眼看出
杂合扩增与嵌合。

## 3. 参考文献

1. TRGT — https://github.com/PacificBiosciences/trgt
2. *TRGT-ing the dark genome to accurately characterize tandem repeats at
   scale*. Nature Biotechnology 42:1504-1505, 2024.
   DOI: 10.1038/s41587-023-02073-3(TRGT 方法介绍)
3. English AC et al. *Benchmarking of small and large variants across
   tandem repeats*. bioRxiv, 2023. DOI: 10.1101/2023.10.29.564632
   (TRGT repeat catalog 来源,预印本)
