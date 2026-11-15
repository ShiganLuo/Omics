# scRNAseq_report — 单细胞分析脚本集(独立脚本,暂未接入 .smk)

两个独立 Python 脚本,读取 annotated h5ad 输出表格与图片;目前仓库内
没有任何 .smk/流程规则引用它们(手动调用)。

## 脚本与输出

| 脚本 | 输出(--output-dir 下) |
|---|---|
| `bin/cell_type_composition.py` | `{label}_counts.tsv`、`{label}_pct.tsv`、`{label}_sample_composition.<fmt>`(stacked bar) |
| `bin/pseudobulk_deg.py` | `{label}/{cell_type}_{comp}_deg.tsv`、`{label}/{cell_type}_{comp}_volcano.<fmt>`(火山图) |

## 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `--input` | `--input`(可重复/多值) | 必填 | annotated h5ad 文件 |
| `--output-dir` | `--output-dir` | 必填 | 表格与图片输出目录 |
| `formats` | `--format`(可重复) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba |

示例:

```bash
python3 modules/scRNAseq_report/bin/cell_type_composition.py \
  --input sample_annotated.h5ad --output-dir out \
  --format png --format pdf --format svg
```
