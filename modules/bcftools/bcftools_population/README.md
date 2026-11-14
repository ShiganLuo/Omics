# bcftools/bcftools_population — 群体遗传 SNP 质控过滤子模块

从联合基因型 VCF 中提取双等位 SNP 并做群体遗传质控过滤(QUAL/DP/缺失率/
MAF),输出 `population_snps.vcf.gz`。conda 环境共享父目录 `../bcftools.yaml`
(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `input_vcf`(config,可选) | 联合基因型 VCF + `.tbi`;缺省为 `{outdir}/variants/joint/joint.vcf.gz`(上游 joint calling 产出) |

### 1.2 输出

```
{outdir}/variants/filtered/population_snps.vcf.gz      # 过滤后群体 SNP
{outdir}/variants/filtered/population_snps.vcf.gz.tbi  # tabix 索引
{logdir}/variants/bcftools_population_filter.log
```

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.bcftools` | 可执行文件本身 | `bcftools` | bcftools 可执行路径 |
| `Params.filtering.maf` | `-i "MAF>=X"` | 0.05 | 次等位基因频率下限 |
| `Params.filtering.max_missing` | `-i "F_MISSING<=1-X"` | 0.9 | 位点最大缺失率(实际以 F_MISSING ≤ 1−max_missing 实施) |
| `Params.filtering.min_qual` | `-e "QUAL<X"` | 30 | 位点质量下限 |
| `Params.filtering.min_dp` / `max_dp` | `-e "DP<X \|\| DP>Y"` | 8 / 200 | 深度下限/上限(过滤极端深度位点) |
| `input_vcf` | — | `{outdir}/variants/joint/joint.vcf.gz` | 联合 VCF 路径 |
| `outdir` / `logdir` | — | `output` / `logs` | 输出/日志根目录 |

资源默认:`threads = 4`、`mem_mb = 8000`(经 `rule_res()` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
bcftools view -m2 -M2 -v snps joint.vcf.gz \
  | bcftools filter -e "QUAL<30 || DP<8 || DP>200" \
  | bcftools view -i "F_MISSING<=0.1 && MAF>=0.05" -Oz -o population_snps.vcf.gz
bcftools index -t population_snps.vcf.gz
```

## 2. 算法原理与生物学意义

三段流式过滤:①`view -m2 -M2 -v snps` 只保留双等位 SNP(剔除 indel/
多等位位点,群体遗传 PCA/ADMIXTURE/关联分析通常只要干净的 biallelic
SNP);②`filter -e` 剔除低质量(QUAL)、深度异常(DP 过低或过高——过高
常是重复区比对堆积)的位点;③`view -i` 保留缺失率低(样本基因型完整度
高)且 MAF 达标的位点(低频位点在群体结构与关联分析中信息量低且易为测序
错误)。结果阅读:输出 VCF 的位点数、Ti/Tv(`bcftools stats`)是质控
金指标(全基因组 biallelic SNP 的 Ti/Tv ≈ 2.0);该产物即下游 plink2
群体分析的输入。

## 3. 参考文献

1. bcftools —— https://samtools.github.io/bcftools/howtos/index.html
2. Danecek P, et al. *Twelve years of SAMtools and BCFtools*. GigaScience
   2021, 10(2):giab008. DOI: 10.1093/gigascience/giab008
