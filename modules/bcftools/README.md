# bcftools — VCF 排序与索引模块

封装 bcftools 两个原子操作:`bcftools sort`(VCF 按坐标排序并 bgzip 压缩)
与 `bcftools index`(生成 CSI 索引),为下游按区间随机访问 VCF 做准备。
纯工具封装模块。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 上游 VCF 目录;`bcftools_sort` 消费 `{indir}/{sample_id}/{sample_id}.vcf`,`bcftools_index_csi` 消费 `{indir}/{sample_id}/{sample_id}.vcf.gz`(形态 A 布局拼接) |
| `samples` | 样本 ID 列表(经 subworkflow 通配符) |

### 1.2 输出

```
{outdir}/sort/{sample_id}/{sample_id}.sorted.vcf.gz      # 排序+压缩 VCF
{outdir}/{sample_id}/{sample_id}.vcf.gz.csi              # CSI 索引
```

收尾规则 `bcftools_result` 聚合上述两个输出。

### 1.3 参数

无 `Params` 配置段;工具参数固定(`sort -O z`)。

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.bcftools` | — | `bcftools` | bcftools 可执行路径 |
| mem_mb | — | 8000(每规则) | `rule_res("bcftools_sort"/"bcftools_index_csi"/"bcftools_result")` 配置 |

直接调用示例(脱离流程):

```bash
bcftools sort -O z -o S1.sorted.vcf.gz S1.vcf
bcftools index S1.vcf.gz     # 生成 S1.vcf.gz.csi
```

## 2. 算法原理与生物学意义

`bcftools sort` 把 VCF 记录按染色体/坐标排序(bcftools 内部先分块排序再
归并,内存受 `--max-mem` 控制),`-O z` 输出 bgzip 压缩的 VCF;`bcftools
index` 对 bgzip VCF 建 CSI 索引,支持 `bcftools view -R/-r`、IGV 等按
区间随机访问。生物学意义:变异检测器(如 GATK/manta)直出的 VCF 可能
未按坐标排序或缺索引,排序 + 索引是进入合并、注释、过滤等下游步骤的
标准前置。纯工具封装。

## 3. 参考文献

1. bcftools — https://samtools.github.io/bcftools/
2. Danecek P et al. *Twelve years of SAMtools and BCFtools*. GigaScience
   10(2):giab008, 2021. DOI: 10.1093/gigascience/giab008
