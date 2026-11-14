# msisensor-pro/tumor-only — 无配对正常样本的 MSI 检测子模块

MSIsensor-pro 的 tumor-only 工作流:先用一组正常样本构建微卫星基线
(baseline),再对肿瘤样本评估微卫星不稳定性(MSI)。该模块为独立形态
(config 顶层键,非 Procedure/Params 段);conda 环境共享父目录
`../msisensor-pro.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `reference`(config) | 参考基因组 FASTA |
| `genome_version` | 基因组版本名(如 GRCh38,用于产物命名) |
| `samples4baseline` | 基线样本表 TSV(pandas 读入,需含 `bam_path` 列;行索引为样本 ID) |
| `samples4test` | 待检样本表 TSV(需含 `tumor_path` 列;行索引为 case ID) |

### 1.2 输出

```
{outdir}/reference/{genome_version}.msisensor.scan.list            # 微卫星位点扫描表
{outdir}/baselines/details/{sample}.{genome_version}.baseline.out  # 逐基线样本分布
{outdir}/baselines/{genome_version}.baseline.samples.list          # 基线样本清单
{outdir}/baselines/{genome_version}.baseline.tsv                   # 基线模型
{outdir}/tumor_only_output/{case}/{case}.{genome_version}.msisensor-pro  # 逐 case MSI 结果
{outdir}/tumor_only_output.{genome_version}.merge.tsv               # 全 case 汇总(最终产物)
```

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `msisensor-pro` | 可执行文件本身 | — | MSIsensor-pro 可执行路径(顶层键,必填) |
| `reference` | `-g` / `scan -d` | — | 参考 FASTA |
| `genome_version` | — | — | 版本名(命名用) |
| `samples4baseline` | — | — | 基线样本 TSV(`bam_path` 列) |
| `samples4test` | — | — | 待检样本 TSV(`tumor_path` 列) |
| `outdir` | — | — | 输出根(日志也写在 outdir 下) |

各规则 `mem_mb` 默认 8000(`pro_pre` threads=2;经 `rule_res()` 可覆盖);
`baseline` 步固定 `-s 1`。模块无 logdir 概念,日志均在 outdir 内。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
msisensor-pro scan -d ref.fa -o ref.msisensor.scan.list
msisensor-pro pro -d ref.msisensor.scan.list -t normal.bam -g ref.fa -o s1.baseline.out
# merge 各基线样本路径 -> baseline.samples.list
msisensor-pro baseline -d ref.msisensor.scan.list -i baseline.samples.list -o baseline.tsv -s 1
msisensor-pro pro -d baseline.tsv -t tumor.bam -g ref.fa -o case1.msisensor-pro
# merge_msi_result 规则汇总为 tumor_only_output.GRCh38.merge.tsv
```

## 2. 算法原理与生物学意义

MSIsensor-pro 对微卫星(短串联重复)位点,统计 reads 支持的重复长度分布,
以泊松模型评估肿瘤相对基线的长度不稳定位点比例,给出 MSI score。tumor-only
场景(无配对正常组织)用一批正常样本(`samples4baseline`)先逐样本
`pro` 建分布(`pro_pre`)、合并(`merge`)、`baseline` 拟合位点基线,再对
每个肿瘤 case(`run_pro`)与基线比较;`merge_msi_result` 汇总每 case 的
`Total_number_of_sites`、`Number_of_unstable_sites`、`MSI_score` 三列。
生物学意义:MSI 是错配修复缺陷(dMMR)的标志,是结直肠癌等肿瘤的免疫
治疗伴随诊断指标;无配对正常时基线法是标准替代。结果阅读:merge.tsv 的
MSI_score 高低判读 MSI-H/MSI-L/MSS(阈值依临床/工具约定);不稳定位点数
占比与 score 一致可交叉验证。

## 3. 参考文献

1. MSIsensor-pro —— https://github.com/xjtu-omics/msisensor-pro(源码与文档)
2. Jia P, Yang X, Guo L, et al. *MSIsensor-Pro: fast, accurate, and
   matched-normal-sample-free detection of microsatellite instability*.
   Genomics, Proteomics & Bioinformatics 2020, 18(1):65–71.
   DOI: 10.1016/j.gpb.2020.02.001
3. Niu B, et al. *MSIsensor: microsatellite instability detection using paired
   tumor-normal sequence data*. Bioinformatics 2014, 30(7):1015–1016.
   DOI: 10.1093/bioinformatics/btt755(单端/配对 MSI 检测前作)
