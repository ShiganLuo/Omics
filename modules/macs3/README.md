# macs3 Module

MACS3 (Model-based Analysis of ChIP-Seq) peak calling for ChIP-seq/DIP-seq data.

## Output

每个 sample 同时输出 narrow peak 和 broad peak：

| 文件 | 说明 |
|------|------|
| `{sample}_peaks.narrowPeak` | Narrow peak (ENCODE narrowPeak 格式, 10列) |
| `{sample}_peaks.xls` | Narrow peak 详细信息 (含 p/q-value) |
| `{sample}_broad_peaks.broadPeak` | Broad peak (ENCODE broadPeak 格式, 9列) |
| `{sample}_broad_peaks.xls` | Broad peak 详细信息 |
| `{sample}_cutoff_analysis.txt` | Narrow peak 阈值分析报告 |
| `{sample}_broad_cutoff_analysis.txt` | Broad peak 阈值分析报告 |

## Narrow Peak vs Broad Peak

### Narrow Peak
- **适用**: 转录因子 (TF)、H3K4me3 等点状富集信号
- **峰宽**: 通常 100-500 bp
- **特点**: 精确定位峰的边界，适合 motif 分析
- **阈值**: `-p` (p-value, 默认 1e-5)

### Broad Peak
- **适用**: H3K27me3、H3K36me3、H3K9me3 等组蛋白修饰信号
- **峰宽**: 数 kb 到数十 kb，覆盖大片染色质区域
- **特点**: 峰边界模糊，定义染色质状态区域（活跃/抑制）
- **阈值**: `--broad-cutoff` (q-value, 默认 0.1)

### 选择依据

| 标记物类型 | 推荐 Peak 类型 | 示例 |
|-----------|---------------|------|
| 转录因子 | Narrow | CTCF, YY1, p53 |
| 活性启动子 | Narrow | H3K4me3 |
| 活性增强子 | Narrow | H3K4me1, H3K27ac |
| 抑制性修饰 | Broad | H3K27me3, H3K9me3 |
| 转录延伸 | Broad | H3K36me3 |

## 配置参数

```json
{
    "Params": {
        "macs3": {
            "bw": 200,
            "pvalue": "1e-5",
            "genome_size": "mm",
            "broad_cutoff": "0.1"
        }
    }
}
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `bw` | 200 | 带宽 (bandwidth), 用于 reads 扩展和 peak 平滑 |
| `pvalue` | 1e-5 | Narrow peak p-value 阈值 |
| `genome_size` | mm | 基因组大小 (hs/mm/ce/dm 或具体数值) |
| `broad_cutoff` | 0.1 | Broad peak q-value 阈值 |

## Cutoff Analysis

两条命令均无条件启用 `--cutoff-analysis`，生成阈值分析报告。该报告列出不同 p-value/q-value 阈值对应的 peak 数量，供人工评估最优阈值。MACS3 不支持自动调参，需根据报告手动调整 `pvalue` 或 `broad_cutoff` 后重跑。
