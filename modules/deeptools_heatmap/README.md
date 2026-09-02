# deeptools_heatmap 模块

computeMatrix + plotHeatmap 一步完成。

## Region 模式

### 1. TSS 模式

配置：`"regions": "tss"`

- BED：每个 transcript 的 TSS ± flank，strand 来自 GTF
- computeMatrix：reference-point + TSS
- TSS 标签：✓ 每行是一个 transcript 的 TSS，strand 可信

### 2. Peaks 模式

配置：`"regions": "peaks"`

- BED：MACS3 narrowPeak，无 strand 列
- computeMatrix：reference-point + center
- TSS/TES 标签：✗ peak 无方向性，标签应为 "Peak Center"

### 3. Genes 模式（合并热图）

配置：`"regions": {"genes": [...], ...}`，`per_gene: false`（默认）

所有基因/TE 合并在一张热图，每行一个区间。

#### 3a. gene GTF（`gtf` 不设或 `"gene"`）

- BED：每个基因一行，exon 合并成 gene body (min_start, max_end)
- strand：取第一个 exon 的 strand（同一基因所有 exon 同 strand）
- computeMatrix：reference-point 或 scale-regions
- TSS/TES 标签：✓ 基因 body 有明确的 TSS/TES

#### 3b. TE GTF（`gtf: "te"`）

- BED：每个 TE locus 独立一行（不 merge），strand 来自该 locus 自身
- 一个 subfamily 的所有 locus 分散在各染色体，每个 locus 有独立 strand
- computeMatrix：scale-regions（推荐）或 reference-point + TSS
- TSS/TES 标签：✓ 每个 locus 是独立的有方向区间，TSS = 5' 端，TES = 3' 端

### 4. Per-gene 模式（每基因独立热图）

配置：`"regions": {"genes": [...], "per_gene": true, ...}`

每个基因/TE 生成一张独立热图，由 `heatmap_gene` rule 处理。

#### 4a. gene GTF

- BED：该基因的所有 exon，每个 exon 独立一行
- computeMatrix：scale-regions
- TSS/TES 标签：✓ 同 3a

#### 4b. TE GTF

- BED：该 subfamily 的所有 locus，每个 locus 独立一行
- computeMatrix：scale-regions
- TSS/TES 标签：✓ 同 3b

---

## 汇总

| 场景 | BED 每行代表 | strand 来源 | TSS/TES 适用 |
|---|---|---|---|
| TSS mode | transcript TSS ± flank | GTF transcript | ✓ |
| Peaks mode | peak 区间 | 无 | ✗ |
| Genes + gene GTF | 基因 body (exons merged) | 第一个 exon | ✓ |
| Genes + TE GTF | 单个 TE locus | 该 locus 自身 | ✓ |
| Per-gene + gene GTF | 基因的所有 exons | 第一个 exon | ✓ |
| Per-gene + TE GTF | subfamily 的所有 loci | 各 locus 自身 | ✓ |

**原则：BED 每行是有方向的独立区间 → TSS/TES 适用；无 strand → 不适用。**

---

## merge 行为

`--merge` 仅影响 gene GTF，TE GTF 始终不 merge。

merge 把一个基因的所有 exon 合并成一个 BED 区间 (min_start, max_end)，区间内包含内含子。这样 scale-regions 模式下热图展示整个基因 body 的信号分布（exon 区域信号高，intron 区域信号低），一个基因一行。

不 merge 时每个 exon 独立一行，一个 10 exon 的基因变成 10 行热图，丢失基因 body 上下文。

TE 没有 intron/exon 结构，且不同 locus 分布在不同染色体，merge 无意义。

---

## 输出文件

| 模式 | 输出 | rule |
|---|---|---|
| TSS | `{sample}_tss_heatmap.png` | heatmap |
| Peaks | `{sample}_peaks_heatmap.png` | heatmap |
| Genes 合并 | `{sample}_genes_heatmap.png` | heatmap |
| Per-gene | `{sample}_{gene_name}_heatmap.png` | heatmap_gene |

---

## 节点控制

```json
{
    "Params": {
        "computeMatrix": {
            "samples": ["Rpp21IP", "Rpp14IP"],
            "regions": {
                "genes": ["L1MdTf", "GSAT_MM"],
                "gtf": "te",
                "match_by": "gene_id",
                "per_gene": true
            }
        }
    }
}
```

- `samples`: 控制哪些 IP 样本生成热图（空 = 全部）
- `per_gene`: true = 每个 TE 一张图，false = 合并一张
- 流程由 node.py 注册的输出文件决定，不由 input 函数决定
