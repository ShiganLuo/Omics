# frip_score — ChIP-seq FRiP 质控指标模块

计算 FRiP(Fraction of Reads in Peaks):峰区内 reads 占比,ChIP-seq 实验
质量的核心指标。conda 环境 `frip_score.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 样本 ID 布局(形态 A) | `{indir}/{sample_id}/{sample_id}.sorted_markdup.bam`:markdup 后的比对 BAM |
| `peaks_indir` + 样本 ID 布局 | `{peaks_indir}/{sample_id}/{sample_id}_peaks.narrowPeak`:该样本的峰文件(MACS2 等) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.FRiP.txt    # 两列:sample_id<TAB>frip_score(6 位小数)
{logdir}/{sample_id}/frip_score.log
```

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.bedtools` | 可执行文件本身 | `bedtools` | 计算峰内 reads |
| `Procedure.samtools` | 可执行文件本身 | `samtools` | 统计总比对 reads |
| `indir` | — | `input` | BAM 根目录 |
| `peaks_indir` | — | `peaks` | 峰文件根目录 |
| `outdir` / `logdir` | — | `output` / `log` | 输出/日志根目录 |

资源默认:`mem_mb = 16000`(经 `rule_res("frip_score")` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则内嵌 shell 的等价命令):

```bash
READS_IN_PEAKS=$(bedtools intersect -a s1.sorted_markdup.bam -b s1_peaks.narrowPeak -u | samtools view -c -)
TOTAL_MAPPED=$(samtools flagstat s1.sorted_markdup.bam | grep 'mapped (' | grep -v primary | head -1 | awk '{print $1}')
awk "BEGIN {printf \"%.6f\", $READS_IN_PEAKS / $TOTAL_MAPPED}"
```

## 2. 算法原理与生物学意义

FRiP = 峰区内 reads 数 / 总比对 reads 数。规则实现:`bedtools intersect
-u` 计数与 narrowPeak 有重叠的 reads(每 read 只计一次),`samtools
flagstat` 的 mapped 行(排除 primary 统计行)取总比对 reads;总 reads 为
0 时输出 0.000000。生物学意义:FRiP 衡量 ChIP 富集效率——抗体好、富集强
的实验 reads 高度集中在峰区;ENCODE 建议 FRiP ≥ 1%(理想 > 10%~20%),
低 FRiP 提示抗体失败或起始量不足。结果阅读:`*.FRiP.txt` 的分数可直接
横向比较样本;与 peak 数、测序深度联看,深度饱和后 FRiP 趋稳。

## 3. 参考文献

1. bedtools intersect —— https://bedtools.readthedocs.io/;
   samtools flagstat —— https://www.htslib.org/doc/samtools-flagstat.html
2. Landt SG, et al. *ChIP-seq guidelines and practices of the ENCODE and
   modENCODE consortia*. Genome Research 2012, 22(9):1813–1831.
   DOI: 10.1101/gr.136184.111(FRiP 质控指标出处)
