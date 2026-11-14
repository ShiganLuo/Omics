# tabix — VCF bgzip 压缩与 tabix 索引模块

封装 htslib 的 bgzip + tabix:把 VCF 压缩为 bgzip 并建 tabix(.tbi)
索引,支持下游按区间随机访问。单规则 `tabix_bgzip` + 收尾
`tabix_result`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | VCF 目录;input 函数拼 `{indir}/{sample_id}/{sample_id}.vcf`(未压缩 VCF) |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.vcf.gz      # bgzip 压缩 VCF
{outdir}/{sample_id}/{sample_id}.vcf.gz.tbi  # tabix 索引
```

收尾规则 `tabix_result` 聚合 vcf.gz + tbi。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.bgzip` | `bgzip <in> -o <out>` | `bgzip` | 压缩工具 |
| `Procedure.tabix` | `tabix -p vcf <out>` | `tabix` | 索引工具(按 VCF preset) |
| mem_mb | — | 8000 | `rule_res("tabix_bgzip")` 配置 |

注:.smk 回执中 tabix 可执行经 `params.tabix` 引用,而规则 params 仅
声明了 `bgzip`;接线时确保 `Procedure`/`params` 提供该键。

直接调用示例(脱离流程):

```bash
bgzip S1.vcf -o S1.vcf.gz
tabix -p vcf S1.vcf.gz
```

## 2. 算法原理与生物学意义

bgzip 是按 BGZF 块格式的 gzip 变体,支持随机访问;tabix 对排序后的
坐标型文件(VCF/BED/GFF)建索引,记录每块的虚拟文件偏移,使
`tabix file.vcf.gz chr1:1000-2000` 式区间查询无需全文件扫描。生物学
意义:VCF 进入 IGV、bedtools、bcftools 区间操作的前置格式要求。纯
工具封装。

## 3. 参考文献

1. htslib(tabix/bgzip)— https://samtools.github.io/htslib/
2. Li H. *Tabix: fast retrieval of sequence features from generic
   TAB-delimited files*. Bioinformatics 27(5):718-719, 2011.
   DOI: 10.1093/bioinformatics/btq671
