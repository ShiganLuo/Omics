# plink2 — 群体遗传格式转换/PCA/GWAS 模块(PLINK 2)

把群体 SNP VCF 转成 PLINK 二进制格式(pgen → bed),可选跑 PCA 与
GWAS(线性/逻辑混合模型)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `input_vcf`(config,可选) | 群体 SNP VCF + `.tbi`;缺省为 `{outdir}/variants/filtered/population_snps.vcf.gz`(上游 bcftools_population 产出) |
| `Params.analysis.phenotype` | GWAS 表型文件(`gwas = true` 时必填,缺失直接抛 `ValueError`) |
| `Params.analysis.covariates` | 协变量文件列表(取第一个传 `--covar`) |

### 1.2 输出

```
{outdir}/analysis/plink/population_snps.{pgen,pvar,psam}   # PLINK2 二进制三件套
{outdir}/analysis/plink/population_snps.{bed,bim,fam}      # PLINK1 兼容三件套(admixture 等消费)
{outdir}/analysis/plink/population_snps.{eigenvec,eigenval}# PCA 结果(可选,pca=true)
{outdir}/analysis/gwas/gwas.done                           # GWAS 标记文件(可选,gwas=true)
{outdir}/analysis/gwas/population.*                        # GWAS 结果(--glm 输出)
{logdir}/analysis/plink2_population_*.log
```

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.plink2` | 可执行文件本身 | `plink2` | PLINK 2 可执行路径 |
| `Params.analysis.pca` | — | true | 是否跑 PCA(决定 `plink2_population_pca` 规则是否注册) |
| `Params.analysis.pca_components` | `--pca N approx` | 10 | 主成分个数 |
| `Params.analysis.gwas` | — | false | 是否跑 GWAS(决定 `plink2_population_gwas` 是否注册) |
| `Params.analysis.phenotype` | `--pheno` | null | 表型文件,gwas=true 必填 |
| `Params.analysis.phenotype_name` | `--pheno-name` | 不传 | 表型列名 |
| `Params.analysis.covariates` | `--covar` | [] | 协变量文件列表(取首个) |
| `Params.analysis.covariate_names` | `--covar-name` | 不传 | 协变量列名 |
| `input_vcf` | `--vcf` | `{outdir}/variants/filtered/population_snps.vcf.gz` | 输入 VCF |
| `outdir` / `logdir` | — | `output` / `logs` | 输出/日志根目录 |

各规则 `mem_mb` 默认 8000(经 `rule_res()` 可覆盖);conda 环境
`plink2_population.yaml`(容器同源)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
plink2 --vcf population_snps.vcf.gz dosage=DS --make-pgen --out population_snps
plink2 --pfile population_snps --make-bed --out population_snps
plink2 --pfile population_snps --pca 10 approx --out population_snps
plink2 --bfile population_snps --pheno pheno.txt --pheno-name BMI \
  --covar covar.txt --covar-name AGE,SEX \
  --glm hide-covar allow-no-covars --out analysis/gwas/population
```

## 2. 算法原理与生物学意义

- **格式转换**:`dosage=DS --make-pgen` 保留 DS 剂量信息转 PLINK2 三件套,
  再 `--make-bed` 导出 PLINK1 兼容 bed/bim/fam(下游 ADMIXTURE 等传统
  工具只认该格式);
- **PCA**:`--pca N approx` 对基因型矩阵做近似主成分分析,特征向量
  (eigenvec)是个体在群体结构主轴上的坐标,用于可视化群体分层、检测离群
  样本,亦作 GWAS 协变量校正群体分层;
- **GWAS**:`--glm` 对每个位点做广义线性模型关联检验(数量/二元表型),
  `hide-covar` 隐藏协变量行、`allow-no-covars` 允许无协变量;输出
  `population.<pheno>.glm.*` 逐位点效应量与 p 值,以 `gwas.done` 标记文件
  驱动 DAG。

结果阅读:PCA 散点图看群体聚类;GWAS 输出按 p 值排序画 Manhattan/QQ,
λ 偏离 1 提示残余群体分层。

## 3. 参考文献

1. PLINK 2.0 —— https://www.cog-genomics.org/plink/2.0/(官方文档)
2. Purcell S, et al. *PLINK: a tool set for whole-genome association and
   population-based linkage analyses*. American Journal of Human Genetics
   2007, 81(3):559–575. DOI: 10.1086/519795
