# admixture — 群体祖先成分分析模块(ADMIXTURE)

对 PLINK 格式 SNP 集跑 ADMIXTURE,对每个祖先群体数 K 输出个体祖先成分
矩阵 Q 与祖先等位基因频率 P。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{outdir}/analysis/plink/population_snps.{bed,bim,fam}` | PLINK 1 二进制三件套(bed 基因型矩阵、bim 位点注释、fam 样本注释),路径固定,由上游 PLINK 模块(plink2_population)产出 |

### 1.2 输出

```
{outdir}/analysis/structure/admixture.K{k}.Q   # 各 K 的个体祖先成分矩阵
{outdir}/analysis/structure/admixture.K{k}.P   # 各 K 的祖先等位基因频率矩阵
{logdir}/analysis/admixture.K{k}.log           # 日志
```

`{k}` 为 wildcard,取值来自 `Params.analysis.admixture_k` 列表(每 K 一个
规则实例)。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.admixture` | 可执行文件本身 | `admixture` | ADMIXTURE 可执行路径 |
| `Params.analysis.admixture_k` | 命令行末位的 K 值 | `[2,3,4,5,6]` | 要拟合的祖先群体数列表(驱动 wildcard) |
| `Params.analysis.admixture_seed` | `--seed` | 2026 | 随机种子(保证可重复) |
| `outdir` / `logdir` | — | `output` / `logs` | 输出/日志根目录 |

资源默认:`threads = 4`(拼为 `-j4`)、`mem_mb = 8000`(经
`rule_res("admixture_population")` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
admixture --seed 2026 -j4 output/analysis/plink/population_snps.bed 3
# 产物:population_snps.3.Q / population_snps.3.P
# (规则把 .Q/.P 重命名聚合为 analysis/structure/admixture.K3.Q/.P)
```

## 2. 算法原理与生物学意义

ADMIXTURE 用极大似然(块松弛 + 拟牛顿加速)估计模型:给定 K 个假想祖先
群体的等位基因频率 P(K × 位点),求每个个体的祖先成分比例 Q(个体 × K),
使观测基因型的似然最大;`--seed` 固定初始化随机性,`-j` 多线程。生物学
意义:Q 矩阵揭示个体的群体祖源构成(混合比例),用于群体分层校正、迁徙
与混合历史推断。结果阅读:`admixture.K{k}.Q` 每行一个样本、每列一个祖先
成分(行和为 1),常用 `pong`/`CLUMPAK` 或堆叠条形图展示;K 的选择通常
看交叉验证误差(ADMIXTURE `--cv`,本模块未启用)与生物学可解释性。

## 3. 参考文献

1. ADMIXTURE(官方主页与文档)—— https://dalexander.github.io/admixture/;
   源码 https://github.com/DAlexander/admixture
2. Alexander DH, Novembre J, Lange K. *Fast model-based estimation of
   ancestry in unrelated individuals*. Genome Research 2009, 19(9):1655–1661.
   DOI: 10.1101/gr.094052.109
