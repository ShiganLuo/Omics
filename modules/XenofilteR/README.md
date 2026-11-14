# XenofilteR — 异种移植(PDX)宿主/污染源交叉 reads 过滤模块

对人源肿瘤-小鼠基质(PDX)等"宿主基因组 + 污染源基因组"双比对数据,
用 XenofilteR 判定并剔除污染源(如小鼠)reads,输出仅保留宿主(人)
reads 的过滤 BAM,消除跨物种污染对变异检出/表达定量的干扰。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `Params.XenofilteR.sample_contamination`(config) | sample_id → `{host, contaminant}` 基因组名(第一列目标/宿主基因组,第二列污染源基因组) |
| `indir`(config) | 双比对 BAM 目录,input 函数按布局拼 `{indir}/{genome}/{sample_id}/{sample_id}.bam`(host 与 contaminant 各一份;host == contaminant 时两者同文件,规则跳过过滤) |
| `Parameters.XenofilteR.MM` | 错配阈值(注意 config 根键为 `Parameters` 而非 `Params`) |
| `Procedure.Rscript` | Rscript 可执行路径 |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.csv      # 输入样本表(host_bam,contaminant_bam 一行)
{outdir}/{genome}/{sample_id}/{sample_id}.bam      # 过滤后 BAM(temp)
{outdir}/{genome}/{sample_id}/{sample_id}.bam.bai  # 过滤后 BAM 索引(temp)
```

host == contaminant 时输出 BAM/BAI 为宿主 BAM 的符号链接(不过滤);
否则 XenofilteR 原始产出先落在
`{outdir}/{sample_id}/Filtered_bams/{sample_id}_Filtered.bam(.bai)`,
规则内 mv 到声明输出。日志:`logdir/{sample_id}/{genome}/XenofilteR.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Parameters.XenofilteR.MM` | `--MM` | 8(smk;R 脚本自身默认 4) | 污染判定的错配阈值;工具建议 150 bp reads 取 8 |
| — | `--inputFile` | — | 样本 CSV:每行 `host_bam,contaminant_bam`(规则内 echo 生成) |
| — | `--outputDir` | `{outdir}/{sample_id}` | XenofilteR 输出目录 |
| — | `--renameSamples` | `{sample_id}` | 输出 BAM 重命名 |
| — | `--workers` | 1(固定) | 并行样本数(BiocParallel SOCK) |

直接调用示例(脱离流程):

```bash
echo "host.bam,contaminant.bam" > S1.csv
Rscript bin/XenofilteR.r \
  --inputFile S1.csv --outputDir out/S1 --renameSamples S1 --MM 8 --workers 1
```

## 2. 算法原理与生物学意义

XenofilteR R 包把样本同时比对到宿主基因组与污染源基因组的两份 BAM 逐
read 比对:与污染源基因组匹配更好(错配数 ≤ `MM_threshold`)的 read 判为
污染并剔除,其余保留为宿主 reads(`bin/XenofilteR.r` 封装
`XenofilteR::XenofilteR()`,BiocParallel `SnowParam` 控制并行)。

生物学意义:PDX 模型中人源肿瘤长在小鼠基质里,测序 reads 会同时来自
两个物种;不剔除鼠源 reads 会造成假突变、表达定量失真。双基因组比对 +
错配判定比单纯比对质量过滤更可靠(污染 reads 在宿主基因组上往往也有
次优比对)。模块的 `MM` 随 read 长度调整(150 bp 建议 8);host 与
contaminant 同基因组时无污染问题,直接透传宿主 BAM。

## 3. 参考文献

1. Kluin RJC, Kemper K, Kuilman T, et al. *XenofilteR: computational
   deconvolution of mouse and human reads in tumor xenograft sequence
   data*. BMC Bioinformatics 19(1):366, 2018.
   DOI: 10.1186/s12859-018-2353-5
2. XenofilteR R 包(方法实现,见论文 Availability of data 一节)。
