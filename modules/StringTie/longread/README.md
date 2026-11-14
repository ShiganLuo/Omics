# StringTie/longread — 长读长转录组组装与定量子模块(StringTie -L)

long-read 转录组三段式:StringTie `-L` 长读长组装 → `--merge` 合并
多样本转录本 → `-e -B` 按合并注释定量,输出基因丰度表与定量 GTF。
父目录资源(环境 `../StringTie.yaml`、容器、`StringTie/README.md` 与
bin/ 嵌合脚本)与短读长主模块共享;本子模块不使用 TEChimeric 脚本。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `upstream_indir` / `indir` (config) | 长读长比对 BAM 目录(上游 config 的 `outdir`,优先取 `upstream_indir`):`{dir}/{sample_id}/{sample_id}.sorted.bam` |
| `samples` (config) | 样本 ID 列表 |
| `genome.gtf` (config) | 参考注释 GTF;缺省 `/dev/null`(无注释从头组装) |
| `merged_gtf` (config) | 定量用合并注释;缺省 `{outdir}/stringtie_merged.gtf`(即本模块 merge 产物) |
| `outdir` / `logdir` (config) | 输出/日志目录 |

### 1.2 输出

```
{outdir}/raw/{sample_id}/{sample_id}.gtf            # 每样本组装转录本
{outdir}/stringtie_merged.gtf                       # 多样本合并注释(stringtie --merge)
{outdir}/{sample_id}/{sample_id}.gene_abund.tab     # 基因丰度表(-A)
{outdir}/{sample_id}/{sample_id}.quant.gtf          # 定量转录本 GTF(-e -B)
```

回执:`stringtie_lr_assemble_{sample_id}.sh`、`stringtie_lr_merge.sh`、
`stringtie_lr_quant_{sample_id}.sh`(各自输出目录下);日志
`{logdir}/{sample_id}/stringtie_lr_*.log`、`{logdir}/stringtie_lr_merge.log`。
`stringtie_lr_result` 为聚合规则。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `upstream_indir`(回退 `indir`) | — | `input` | BAM 目录 |
| `genome.gtf` | `-G` | `/dev/null` | 参考注释(组装辅助/合并参照) |
| `merged_gtf` | `-G`(quant) | `{outdir}/stringtie_merged.gtf` | 定量注释 |
| `Procedure.stringtie` | 引擎可执行 | `stringtie` | — |
| (threads) | `-p` | 8(assemble)/ 4(quant) | `rule_res` 可覆盖 |
| (mem_mb) | — | 32000(assemble)/ 16000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
# 1) 每样本长读长组装
stringtie -L -p 8 -G annotation.gtf -o sample.gtf sample.sorted.bam

# 2) 多样本合并
stringtie --merge s1.gtf s2.gtf -o stringtie_merged.gtf -G annotation.gtf

# 3) 按合并注释定量(-e 只定量已知转录本,-B 输出 Ballgown 表)
stringtie -e -B -p 4 -G stringtie_merged.gtf \
    -o sample.quant.gtf -A sample.gene_abund.tab sample.sorted.bam
```

## 2. 算法原理与生物学意义

StringTie `-L` 针对长读长(ONT/PacBio)比对做转录本组装:长读长
直接跨越多个外显子,剪接结构由单条读段串起,组装歧义远低于短读长,
能恢复全长异构体与复杂剪接。三段式的语义:

- **assemble**(`-L` + `-G`):逐样本预测转录本结构(有注释则引导);
- **merge**(`--merge`):把各样本 GTF 合并为统一位点级注释,消除
  样本间命名/碎片差异,是跨样本可比定量的前提;
- **quant**(`-e -B`):只对合并注释中的转录本定量(不预测新结构),
  `-A` 输出基因丰度表(gene_id、参考/匹配长度、TPM、FPKM、cov),
  `-B` 输出 Ballgown 兼容表。

生物学意义:基因丰度表(TPM/FPKM)供差异表达与表达谱分析;合并
GTF 本身即长读长转录组的结构产物(新异构体、新剪接位点)。注意
`-e` 模式的定量强依赖 merged_gtf 的完整性。

## 3. 参考文献

1. StringTie — https://github.com/gpertea/stringtie
2. Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL.
   *StringTie enables improved reconstruction of a transcriptome from
   RNA-seq reads*. Nature Biotechnology 33(3):290-295, 2015.
   DOI: 10.1038/nbt.3122.
3. Kovaka S, Zimin AV, Pertea GM, Razavi R, Salzberg SL, Pertea M.
   *Transcriptome assembly from long-read RNA-seq alignments with
   StringTie2*. Genome Biology 20:275, 2019.
   DOI: 10.1186/s13059-019-1910-9.(`-L` 长读长模式方法出处)
