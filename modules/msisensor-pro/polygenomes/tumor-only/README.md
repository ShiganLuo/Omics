# msisensor-pro/polygenomes/tumor-only — 肿瘤 MSI 检测子模块(MSIsensor-pro,无配对正常)

MSIsensor-pro 基线(baseline)模式:用一组正常/基线样本建立微卫星不稳定
(MSI)基线,再对肿瘤样本做无配对正常(tumor-only)MSI 打分。C 形态嵌套
子模块(msisensor-pro/polygenomes/tumor-only),自带 `rule all`,可独立
作 Snakefile 入口。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `samples4baseline`(config) | 基线样本 TSV(`pd.read_table(index_col=0)`),需含 `bam_path` 列;每样本一个基因组比对 BAM |
| `samples4test` | 待测肿瘤样本 TSV(`index_col=0`),需含 `tumor_path` 列 |
| `reference` | 参考基因组 FASTA(平铺键;配合 `genome_version` 构建 genome→fasta 映射)。json 另声明 `genome.<g>.fasta` 字典,本 .smk 不消费 |
| `genome_version` | 基因组版本名,进入输出命名与 `rule all` 目标 |

样本清单以 TSV 形态传入(非 `samples` 列表),模块直接用 pandas 索引列。

### 1.2 输出

```
{outdir}reference/{genome_version}.msisensor.scan.list               # 微卫星位点列表(scan)
{outdir}baselines/details/{sample}.{genome_version}.baseline.out      # 基线样本特征
{outdir}baselines/{genome_version}.baseline.samples.list              # 基线样本清单
{outdir}baselines/{genome_version}.baseline.tsv                       # MSI 基线模型
{outdir}tumor_only_output/{case}/{case}.{genome_version}.msisensor-pro # 逐肿瘤样本 MSI 结果
{outdir}tumor_only_output.{genome_version}.merge.tsv                  # 合并结果(rule all 目标)
```

日志:同名 `.log` 旁挂在输出旁。环境:`conda: "../../msisensor-pro.yaml"`
(conda 环境共享祖父目录)。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `outdir` | — | `output` | 输出目录(必需;尾部自动补 `/`) |
| `samples4baseline` | `pro -t <bam>`(逐样本) | — | 基线样本 TSV(`bam_path` 列) |
| `samples4test` | `pro -t <bam>`(逐 case) | — | 待测样本 TSV(`tumor_path` 列) |
| `reference` | `-d` / `-g` | — | 参考 FASTA |
| `genome_version` | 命名通配符 | `genome1`(json 示例) | 版本名 |
| `msisensor-pro` | — | `msisensor-pro` | 引擎可执行路径(顶层键,非 `Procedure` 字典) |

固定调用:`scan -d <fa> -o <list>`;基线 `pro -d <list> -t <bam> -g <fa>`;
`baseline -d <list> -i <samples.list> -o <baseline.tsv> -s 1`;测试
`pro -d <baseline.tsv> -t <tumor.bam> -g <fa>`。资源:`rule_res(...)` 可覆盖
(默认 `mem_mb=8000`,`pro_pre`/`run_pro` 默认 `threads=2`/按 rule_res)。

直接调用示例(等价引擎命令):

```bash
msisensor-pro scan  -d genome.fa -o genome.msisensor.scan.list
msisensor-pro pro   -d genome.msisensor.scan.list -t normal.bam -g genome.fa \
                    -o sample.baseline.out
msisensor-pro baseline -d genome.msisensor.scan.list -i baseline.samples.list \
                    -o baseline.tsv -s 1
msisensor-pro pro   -d baseline.tsv -t tumor.bam -g genome.fa \
                    -o case.genome1.msisensor-pro
```

## 2. 算法原理与生物学意义

MSIsensor-pro 扫描参考基因组微卫星(短串联重复)位点,统计每样本各位点的
重复长度分布;baseline 模式用一组基线样本学习各微卫星的种系长度分布,再对
无配对肿瘤样本计算体细胞长度偏移的似然比,输出 MSI score 与显著性。
MSI 是错配修复缺陷(dMMR)的标志物,是结直肠癌等肿瘤的免疫治疗伴随诊断
指标;tumor-only 模式免配对正常组织,适合回顾性队列。结果怎么读:
逐样本 `.msisensor-pro` 输出含 MSI score(score 越高越不稳定;常用阈值
20/35 分档,以本地规范为准),`merge.tsv` 汇总全队列便于排序筛选。

## 3. 参考文献

1. MSIsensor-pro — https://github.com/xjtu-omics/msisensor-pro
2. Jia P et al. *MSIsensor-pro: fast, accurate, and matched-normal-sample-free
   detection of microsatellite instability*. Genomics, Proteomics &
   Bioinformatics 18(1):65–71, 2020. DOI: 10.1016/j.gpb.2020.02.005
