# gatk/polygenomes/gatk_germline — 多参考基因组 GATK 胚系变异检测子模块

GATK HaplotypeCaller + VariantFiltration 的多参考基因组(polygenomes)变体:
按 `genome` wildcard 对每样本 × 每基因组调用胚系变异并硬过滤。conda 环境
共享 `../../gatk.yaml`(即 `modules/gatk/gatk.yaml`,容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 样本 ID 布局(形态 A) | `reference.known_sites` 与 `reference.interval` 同时提供时取 BQSR 产物 `{indir}/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam`(该分支在 input 函数内 `include: "../gatk_bqsr/gatk_bqsr.smk"` 引入 BQSR 规则);否则回退 `{indir}/{sample_id}/{sample_id}.sorted_markdup.bam` |
| `genome.<G>.fasta` / `.fai_index` / `.dict_index` | 每基因组参考 FASTA、.fai、.dict;三者都配置时直接使用,任一缺失则回退到 `{indir}/index/{genome}/genome.fa(.fai/.dict)` |
| `reference.interval` / `reference.known_sites` | 区间文件与已知位点列表(此处仅用于选择输入 BAM 形态) |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.raw.vcf.gz          # HaplotypeCaller 原始调用
{outdir}/{genome}/{sample_id}/{sample_id}.filtered.vcf.gz     # 过滤后 VCF
{outdir}/{genome}/{sample_id}/{sample_id}.filtered.vcf.gz.csi # bcftools index
{logdir}/{genome}/{sample_id}/haplotypeCaller.log
{logdir}/{genome}/{sample_id}/haplotypeCaller-filtered.log
```

收尾规则 `gatk_haplotypeCaller_result` 聚合 filtered VCF + CSI。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.gatk` | 可执行文件本身 | `gatk` | GATK 可执行路径 |
| `Procedure.bcftools` | — | `bcftools` | 声明的 bcftools 路径;过滤规则实际以字面量 `bcftools index` 调用 |
| `Params.gatk.javaOptions` | `--java-options` | `-Xmx30g` | JVM 参数 |
| `Params.gatk.tmp-dir` | `--tmp-dir` | 不传 | GATK 临时目录 |
| `Params.gatk.FS_threshold` | `--filter "FS > X"` | 20.0 | Fisher 链偏倚过滤阈值 |
| `Params.gatk.QD_threshold` | `--filter "QD < X"` | 2.0 | 质量/深度比过滤阈值 |
| `Params.gatk.DP_threshold` | `--filter "DP < X"` | 10.0 | 位点深度过滤阈值 |
| `Params.gatk.QUAL_threshold` | `--filter "QUAL < X"` | 20.0 | 位点质量过滤阈值 |
| `genome.<G>.fasta/.fai_index/.dict_index` | `-R` | — | 每基因组参考三件套 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

资源默认:HaplotypeCaller threads=10、mem_mb=32000;
filterHaplotypeCallerVcf threads=10、mem_mb=16000(均经 `rule_res()` 可
覆盖)。过滤固定附加 `-window 35 -cluster 3`(35 bp 窗口内 ≥3 个变异聚簇
判为伪迹)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
gatk --java-options "-Xmx30g" HaplotypeCaller \
  -R genome.fa -I s1.sorted_markdup.bam -O s1.raw.vcf.gz

gatk --java-options "-Xmx30g" VariantFiltration \
  -R genome.fa -V s1.raw.vcf.gz -O s1.filtered.vcf.gz \
  -window 35 -cluster 3 \
  --filter-name FS   --filter "FS > 20.0" \
  --filter-name QD   --filter "QD < 2.0" \
  --filter-name DP   --filter "DP < 10.0" \
  --filter-name QUAL --filter "QUAL < 20.0"
bcftools index s1.filtered.vcf.gz
```

## 2. 算法原理与生物学意义

- **HaplotypeCaller**:对活跃区域做局部 de novo 组装重建单倍型,再比对
  reads 计算基因型似然,输出 gVCF/VCF;对 indel 与复杂位点远优于纯 pileup
  调用;
- **VariantFiltration 硬过滤**:按注释字段打 FILTER 标签——`FS`(Fisher
  链偏倚,链伪迹)、`QD`(质量/深度,低深度虚高质量)、`DP`(深度不足)、
  `QUAL`(位点质量不足),外加 `-window 35 -cluster 3` 的聚簇过滤
  (比对困难区常成簇出现假阳性)。FILTER 列记录每个位点命中的过滤器名,
  `PASS` 为全部通过。

生物学意义:胚系变异(SNP/indel)是遗传分析、单倍型与疾病关联的基础;
硬过滤是在无 truth set 场景下压低假阳性的标准做法。polygenomes 场景下同
一样本对多参考基因组分别调用(如人/鼠跨物种样本),用于物种归属与污染
排查。结果阅读:`bcftools stats` 看 Ti/Tv 与 indel 比例;FILTER 分布评估
过滤松紧;raw 与 filtered 的差异即被过滤伪迹。

## 3. 参考文献

1. GATK 官方文档(HaplotypeCaller / VariantFiltration)——
   https://gatk.broadinstitute.org/
2. McKenna A, et al. *The Genome Analysis Toolkit: a MapReduce framework
   for analyzing next-generation DNA sequencing data*. Genome Research 2010,
   20(9):1297–1303. DOI: 10.1101/gr.107524.110
