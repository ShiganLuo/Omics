# vcftools — 群体遗传多样性与 Fst 分析模块(vcftools_population)

用 vcftools 对群体 SNP VCF 计算每群体窗口核苷酸多样性 π、Tajima's D,以及
成对群体 Weir-Cockerham Fst。本目录 stem 为 `vcftools_population`
(`vcftools_population.smk/.json/.schema.json/.yaml`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `input_vcf`(config) | 群体 VCF,bgzip 压缩 `.vcf.gz` 且带 `.tbi` 索引(input 同时声明 `vcf` 与 `vcf.tbi` 保证 DAG 依赖);缺省 `{outdir}/variants/filtered/population_snps.vcf.gz`,即上游群体变异过滤产出 |
| `populations`(config) | 群体名字典(键 = 群体名);模块只取键名生成通配符目标。每群体的 keep 文件 `{outdir}/metadata/populations/{population}.txt`(vcftools `--keep` 约定:每行一个样本 ID)由上游 metadata 产出,作为 input 锚点 |

样本 → 群体的归属写在 keep 文件里,模块自身不解析样本清单。

### 1.2 输出

```
{outdir}/analysis/diversity/{population}.windowed.pi              # 窗口核苷酸多样性 π
{outdir}/analysis/diversity/{population}.Tajima.D                 # 窗口 Tajima's D
{outdir}/analysis/fst/{population_a}_vs_{population_b}.weir.fst   # 成对群体 Fst
```

日志:`{logdir}/analysis/vcftools_diversity_{population}.log`、
`{logdir}/analysis/vcftools_fst_{population_a}_vs_{population_b}.log`。
成对 Fst 的群体对通配符目标由 subworkflow/node.py 注册(模块含
`itertools.combinations` 导入)。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `input_vcf` | `--gzvcf` | `{outdir}/variants/filtered/population_snps.vcf.gz` | 群体 VCF |
| `populations` | `--keep`(每群体) | `{}` | 群体名字典;空则无目标 |
| `Params.analysis.diversity` | — | `true` | 关闭时 parse 期不注册 diversity 规则 |
| `Params.analysis.fst` | — | `true` | 关闭或群体数 ≤ 1 时不注册 fst 规则 |
| `Params.analysis.window_size` | `--window-pi` / `--TajimaD` | `100000` | 窗口大小(bp) |
| `Params.analysis.window_step` | `--window-pi-step` | `50000` | 窗口步长(bp) |
| `Procedure.vcftools` | — | `vcftools` | vcftools 可执行路径 |

环境:`conda: "vcftools_population.yaml"`(模块目录内)。

直接调用示例(等价引擎命令):

```bash
# 多样性(π + Tajima's D)
vcftools --gzvcf pop.vcf.gz --keep pops/A.txt \
  --window-pi 100000 --window-pi-step 50000 --out analysis/diversity/A
vcftools --gzvcf pop.vcf.gz --keep pops/A.txt \
  --TajimaD 100000 --out analysis/diversity/A
# 成对 Fst
vcftools --gzvcf pop.vcf.gz \
  --weir-fst-pop pops/A.txt --weir-fst-pop pops/B.txt \
  --out analysis/fst/A_vs_B
```

## 2. 算法原理与生物学意义

vcftools 是群体遗传学标准工具,本模块消费其三类统计:

- **核苷酸多样性 π**(`--window-pi`):滑动窗口内所有成对序列差异位点比例的
  期望值,反映群体近期有效群体大小与变异丰富度;π 低的区域可能受选择清除。
- **Tajima's D**(`--TajimaD`):窗口内 π 与分离位点数 S 推导的 θ_W 之差的
  标准化统计量,D ≈ 0 为中性;显著正值提示平衡选择/群体收缩,显著负值提示
  扫荡选择/群体扩张。
- **Weir-Cockerham Fst**(`--weir-fst-pop` × 2):成对群体间遗传分化指数,
  Fst → 0 两群体无分化,→ 1 完全固定差异;高 Fst 窗口是群体特异选择或
  局部适应的候选区。

结果均为 vcftools 原生表格式(`.windowed.pi` / `.Tajima.D` / `.weir.fst`),
首列为染色体与窗口起止,逐行读即可;`.weir.fst` 逐位点给出 Fst,可再聚合。

## 3. 参考文献

1. vcftools — https://vcftools.github.io/ ,https://github.com/vcftools/vcftools
2. Danecek P et al. *The variant call format and VCFtools*.
   Bioinformatics 27(15):2156–2158, 2011. DOI: 10.1093/bioinformatics/btr330
