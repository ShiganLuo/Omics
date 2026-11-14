# msisensor-pro/tumor-normal — 肿瘤-正常配对微卫星不稳定性(MSI)检测模块

封装 MSIsensor-pro 的 tumor-normal 模式:扫描参考基因组微卫星位点,
用配对正常样本作基准判定肿瘤 MSI 状态,汇总为跨病例表。本模块为
MSIsensor-pro 作者(Peng Jia)提供的 baseline 构建风格实现。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `sample_path`(config) | 病例样本表 TSV(pandas `read_table(index_col=0)` 读取;索引 = case 名,列 `tumor_path` / `normal_path` 为配对 BAM 路径) |
| `reference`(config) | 参考基因组 FASTA(`scan`/`msi` 的 `-d`/`-g`) |
| `genome_version`(config) | 基因组版本名(如 `GRCh38`),进输出路径与文件名 |
| `msisensor-pro`(config) | MSIsensor-pro 可执行路径 |

### 1.2 输出

```
{outdir}/reference/{genome_version}.msisensor.scan.list               # 微卫星位点表(scan)
{outdir}/tumor_normal_output/{case}/{case}.{genome_version}.msisensor  # 每病例 MSI 结果
{outdir}/tumor_normal_output.{genome_version}.merge.tsv               # 全病例汇总
```

(`merge.tsv` 列:`Total_number_of_sites`、`Number_of_unstable_sites`、
`MSI_score`,取每病例 `.msisensor` 输出的统计行。)

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `msisensor-pro` | — | `msisensor-pro` | 可执行路径 |
| `genome_version` | — | `GRCh38` | 基因组版本名(输出命名) |
| `reference` | `scan -d` / `msi -g` | — | 参考基因组 FASTA |
| `sample_path` | `msi -t/-n` | — | 样本表(提供 tumor/normal BAM) |
| — | `scan -d <ref> -o <list>` | 固定 | 扫描微卫星 |
| — | `msi -d <list> -n <normal> -t <tumor> -g <ref> -o <out>` | 固定 | 配对 MSI 判定 |
| mem_mb | — | scan/merge 8000,msi 16000 | `rule_res` 配置 |

conda 环境按 `../msisensor_pro.schema.yaml` 解析(父目录;本子目录另有
同名 `msisensor_pro.schema.yaml`)。

直接调用示例(脱离流程):

```bash
msisensor-pro scan -d genome.fa -o genome.msisensor.scan.list
msisensor-pro msi -d genome.msisensor.scan.list \
  -n normal.bam -t tumor.bam -g genome.fa -o case1.msisensor
```

## 2. 算法原理与生物学意义

MSIsensor-pro 先在参考基因组上扫描微卫星(短串联重复)位点(`scan`),
再在配对样本中对这些位点做长度分布统计:以正常样本的分布为基准,
肿瘤中长度分布显著偏移的位点判为不稳定,`MSI_score` = 不稳定位点数 /
总位点数(百分比),超过阈值(通常 10~20%)判 MSI-H。生物学意义:MSI
是错配修复缺陷(dMMR)的标志,是结直肠癌等肿瘤的预后与免疫治疗
(PD-1 抑制剂)疗效预测 biomarker;tumor-normal 模式用配对正常样本消
除个体多态背景,判定更准。结果读法:`.msisensor` 给出总位点数、
不稳定位点数、MSI score;`merge.tsv` 横向比较全队列。

## 3. 参考文献

1. Jia P, Yang X, Guo L, et al. *MSIsensor-pro: Fast, Accurate, and
   Matched-normal-sample-free Detection of Microsatellite Instability*.
   Genomics Proteomics Bioinformatics 18(1):65-71, 2020.
   DOI: 10.1016/j.gpb.2020.02.001
2. MSIsensor-pro — https://github.com/xjtu-omics/msisensor-pro
3. Niu B, Ye K, Zhang Q, et al. *MSIsensor: microsatellite instability
   detection using paired tumor-normal sequence data*. Bioinformatics
   30(7):1015-1016, 2014. DOI: 10.1093/bioinformatics/btt755
