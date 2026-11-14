# pbsv — PacBio 长读长结构变异检测模块

用 pbsv 从比对 BAM 检出结构变异(SV):discover → call → result 三规则。

## 1. 模块输入输出与参数

### 1.1 输入

形态 A:`indir` + `bam_substring` 布局参数:

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}{<.bam_substring>}.bam` | 长读长比对 BAM;`bam_substring` 非空时为 `{sample_id}.{bam_substring}.bam` |
| 同名 `.bai` | BAM 索引(与 BAM 同目录、同 stem) |
| `genome.fasta` | 参考基因组 FASTA(pbsv call 需要) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.svsig.gz        # discover 结构变异证据
{outdir}/{sample_id}/{sample_id}.sv.vcf.gz       # SV 变异集(VCF,bgzip 压缩)
{outdir}/{sample_id}/{sample_id}.sv.vcf.gz.csi   # bcftools index 索引
{outdir}/{sample_id}/pbsv_discover.sh            # 回执(pbsv_call.sh 同目录)
logdir/{sample_id}/pbsv_discover.log             # 运行日志(pbsv_call.log 同目录)
```

中间未压缩 `{sample_id}.sv.vcf` 在 bgzip 后被 `rm` 清理。`pbsv_result`
规则无 output,只消费 vcf+索引作为 DAG 终点锚点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.pbsv` | 可执行文件 | null(回退 `pbsv`) | pbsv 可执行路径 |
| `Procedure.bgzip` | 可执行文件 | null(回退 `bgzip`) | bgzip 路径 |
| `Procedure.bcftools` | 可执行文件 | null(回退 `bcftools`) | bcftools 路径(建索引) |
| `genome.fasta` | pbsv call 位置参数 | — | 参考基因组 FASTA |
| `bam_substring` | — | "" | 输入 BAM 名中样本名后的后缀(如 `aligned` → `{s}.aligned.bam`) |
| `samples` | — | [] | 样本 ID 列表 |
| `indir` | — | `input` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

线程/内存经 `rule_res(...)` 配置:pbsv_discover threads 4 / mem 16000,
pbsv_call threads 8 / mem 16000。

引擎命令(回执脚本原文):

```bash
pbsv discover {in.bam} {outdir}/{s}/{s}.svsig.gz
pbsv call --num-threads 8 {fasta} {outdir}/{s}/{s}.svsig.gz {outdir}/{s}/{s}.sv.vcf
bgzip {outdir}/{s}/{s}.sv.vcf -o {outdir}/{s}/{s}.sv.vcf.gz
bcftools index {outdir}/{s}/{s}.sv.vcf.gz
```

## 2. 算法原理与生物学意义

纯工具封装:pbsv 分两步检出结构变异——`discover` 扫描 BAM 中长读长的
split/soft-clip/大插入片段信号,写成 svsig 证据文件;`call` 聚类证据并
结合参考序列调用 DEL/INS/DUP/INV/BND 等 SV,输出 VCF。长读长跨越重复
区与复杂结构的能力使其对中大片段 SV 的断点定位显著优于短读长。输出 VCF
可直接进 IGV 核查或下游群体 SV 分析(pbsv call 支持多样本联合,本模块为
单样本模式)。

## 3. 参考文献

1. pbsv(PacBio)— https://github.com/PacificBiosciences/pbsv
2. bgzip/bcftools(htslib)— https://www.htslib.org/
