# pbsv/polygenomes — 多基因组结构变异检测子模块(pbsv)

polygenomes(每基因组独立目录)布局下的 pbsv 长读长 SV 检测:每样本、
每基因组各做 discover → call,输出 bgzip 压缩并建索引的 SV VCF。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | 样本 BAM 目录(上游 config 的 `outdir`):`{indir}/{sample_id}/{sample_id}{.bam_substring}.bam` 与同名 `.bai` 索引;注意 input 拼接**不含** genome 段,同一 BAM 对每个 genome 各跑一次(布局由 subworkflow 保证) |
| `bam_substring` (config,顶层) | BAM 文件名后缀段(如 `sorted` → `sample.sorted.bam`);缺省空 |
| `samples` (config) | 样本 ID 列表 |
| `genome.<genome_id>.fasta` (config) | 每基因组参考 FASTA 字典(pbsv_call 按 wildcards.genome 取用) |
| wildcards `genome` / `sample_id` | 基因组 ID 与样本 ID |

环境共享父目录:`conda: "../pbsv.yaml"`、容器 `sif("../pbsv.yaml", ...)`
均指向 `modules/pbsv/pbsv.yaml`(与非 polygenomes 版共享)。

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.svsig.gz          # pbsv discover 证据信号
{outdir}/{genome}/{sample_id}/{sample_id}.sv.vcf.gz         # pbsv call SV VCF(bgzip)
{outdir}/{genome}/{sample_id}/{sample_id}.sv.vcf.gz.csi     # bcftools index 索引
{outdir}/{genome}/{sample_id}/pbsv_discover.sh              # 回执
{outdir}/{genome}/{sample_id}/pbsv_call.sh                  # 回执
```

日志:`{logdir}/{genome}/{sample_id}/pbsv_discover.log`、`pbsv_call.log`。
`pbsv_result` 为聚合规则(消费 VCF+索引)。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `bam_substring` | — | "" | 见输入 |
| `genome.<genome_id>.fasta` | `pbsv call` 第 1 位置参数 | — | 每基因组参考 FASTA |
| `Procedure.pbsv` | 引擎可执行 | `pbsv` | — |
| `Procedure.bgzip` | `bgzip` | `bgzip` | VCF 压缩 |
| `Procedure.bcftools` | `bcftools index` | `bcftools` | 建索引(.csi) |
| (threads) | `--num-threads` | 4(discover)/ 8(call) | `rule_res` 可覆盖 |
| (mem_mb) | — | 16000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
pbsv discover --num-threads 4 sample.bam sample.svsig.gz

pbsv call --num-threads 8 genome.fa sample.svsig.gz sample.sv.vcf
bgzip sample.sv.vcf -o sample.sv.vcf.gz
bcftools index sample.sv.vcf.gz
```

## 2. 算法原理与生物学意义

pbsv 是 PacBio 长读长(SMRT/HiFi)结构变异检测器:

- **discover**:扫描长读长比对中的 split-read 与大插入/缺失信号,
  逐样本汇总为 `.svsig.gz` 证据文件(每候选断点的读段证据);
- **call**:以参考 FASTA 校正断点序列,对证据聚类判为 DEL/INS/DUP/
  INV/BND 等 SV 类型,输出 VCF。

生物学意义:长读长能跨越重复区与复杂结构,检出短读长遗漏的 SV;
polygenomes 布局下每个基因组一套 VCF,可比较同一读段在不同参考
(如宿主/病原、不同组装版本)下的 SV 调用差异。VCF 经 bgzip+csi
索引后供 IGV 浏览与下游 sv 模块(MEI/TRA/OncoPrint 等分析)消费。

## 3. 参考文献

1. pbsv(Pacific Biosciences)—— https://github.com/PacificBiosciences/pbsv
