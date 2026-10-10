# Fiber-seq 下游分析：需要补齐的文件清单

对应 `config/Fiberseq.json`。所有路径留空（null）的键按下表补齐后，把对应分析的
`Params.<x>.enabled` 设为 `true` 即可（已在默认配置中开启的项也需要文件）。
缺文件时 DAG 会报错并提示具体 JSON 键名。

## 1. genome.references.<组装>（按你选用的组装逐个填）

| JSON 键 | 文件 | 用途 | 建议来源 |
|---|---|---|---|
| `fasta` | 基因组 FASTA（已有） | 全部 | 已填 |
| `fai` | FASTA 索引（`.fai`） | ft call-peaks FDR 模式（shuffle null）、bigWig 轨道 | `samtools faidx <fasta>` 生成 |
| `gtf` | 基因注释 GTF | peak 注释、TSS 核小体定位 | GENCODE/Ensembl 对应组装 |
| `motif_file` | MEME 格式 motif 文件 | motif 扫描（任意物种） | JASPAR MEME 导出（如 CTCF、NFY 等目标 TF） |
| `repeat_bed` | 重复注释 BED（name 列=重复类型） | peak 重复注释（可选） | RepeatMasker 输出转 BED |
| `censat_bed` | 区域类别 BED（name 列=类别） | 着丝粒专题（可选） | T2T 组装用 chm13v2.0_censat_v2.1.bed；其他组装自定义区域集（如 CDR_core / alphaSat_HOR / flank_1kb） |
| `mcpg_pileup` | modkit pileup bedGraph | 着丝粒甲基化定量（可选，ONT 数据推荐） | `modkit pileup --cpg` 输出 |

## 2. Params

| JSON 键 | 文件 | 格式 |
|---|---|---|
| `diff.design_tsv` | 比较设计表 | 三列 TSV：`sample`、`group`、`contrast`（同一样本可出现在多个 contrast）。开启 `diff.enabled` 必填 |

`diff.annot_sample`（可选）：用哪个样本的 peak 注释做差异结果分层（默认第一个样本）。

## 3. Params.function（复用 modules/function 模块）

GO/KEGG + GSEA 走已有的 `function_go_kegg` / `function_gsea` 规则（R
clusterProfiler，离线，支持 Hs/Mm/Mmu）。diff_accessibility 自动生成
`{contrast}.TEcount_Gene.name.tsv` 喂给 go-kegg.r。

| JSON 键 | 文件 | 说明 |
|---|---|---|
| `Params.function.enabled` | — | true 时按 contrast 出 GO/KEGG 上下调结果 |
| `Params.function.species` | — | human / mouse / rhesus |
| `Params.function.gmt` + `genome.references.<asm>.geneIDAnno` | GMT、基因 ID 注释 | 两者都填才跑 GSEA |
| `group_pairs`（顶层，可选） | — | function 的对比标注；不填时从 design_tsv 自动派生（与 diff 的组顺序一致） |

## 4. Params.ldsc（保留接口）

| JSON 键 | 文件 | 说明 |
|---|---|---|
| `annot_bed` | 待检验 peak 集 BED | 如合并后的 FIRE peaks；`annot_cluster_from_name: true` 时按 name 列分 cluster |
| `gwas` | `{性状名: sumstats路径}` 字典 | 每个 GWAS 汇总统计一个条目（加条目即多跑一个性状） |
| `merge_alleles` | HapMap3 SNP 列表 | `munge_sumstats.py --merge-alleles` 用（如 w_hm3.snplist） |
| `ref_ld_chr` | LD score 参考目录前缀 | 如 `eur_w_ld_chr/`（LDSC 官方提供） |
| `w_ld_chr` | 回归权重目录前缀 | 同上 |

## 5. 单倍型分析（Params.haplotype）

无需外部文件——变异检测（DeepVariant）、SV（pbsv）、分相（hiphase）、按 HP 拆分
提取均在流程内完成。需要的 SIF：`deepvariant`、`pbsv`、`hiphase`（env 映射已列）。

## 6. 需要构建的 SIF 容器

| env 键 | 来源 yaml | 说明 |
|---|---|---|
| `fiberseq` | modules/fiberseq/fiberseq.yaml | 下游分析家族唯一环境（python+numpy/scipy/matplotlib/pysam+UCSC tools，布局 C 共享） |
| `function` | modules/function/function.yaml | GO/KEGG + GSEA（已有模块） |
| `ldsc` | modules/ldsc/ldsc.yaml | LDSC |
| `deepvariant`/`pbsv`/`hiphase` | 各模块 yaml | 单倍型链（已有模块，按现有方式构建） |

构建命令（同 FIBERSEQ_ANALYSIS.md）：
```bash
python workflow/Omics/src/common/util/EnvUtil.py apptainer all \
    -y workflow/Omics/modules/<module>/<module>.yaml -o /home/luosg/Database/env/<module> --force
```

## 7. 输出速查（每个分析一个模块，布局 `results/fiberseq_<模块>/<样本|对比>/`）

- `fiberseq_annotate/` peak 注释 · `fiberseq_motif/` motif 位点/谱 · `fiberseq_footprint/` 三分类 footprint
- `fiberseq_nucpos/` 核小体定位 · `fiberseq_quant/` 区域定量+matrix · `fiberseq_track/` bigWig+trackDb
- `fiberseq_coactuation/` 共活化 · `fiberseq_haplotype/` 单倍型差异 · `fiberseq_censat/` 着丝粒专题
- `fiberseq_diff/<contrast>/` 组间差异+基因表 · `function/<contrast>/` GO/KEGG · `output/ldsc/`
