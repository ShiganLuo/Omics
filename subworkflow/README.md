# RNA-SNP/subworkflow 目录说明

本目录包含多个 Snakemake 子工作流（subworkflow），用于不同类型的转录组/表观组数据分析。每个 .smk 文件为一个分析流程主入口，集成了多个标准化模块。

---

## 各工作流简介

### 1. CLIP.smk

- **用途**：iCLIP/CLIP-seq 数据分析全流程。
- **主要模块**：
  - fastqc_raw：原始数据质控
  - cutadapt：去接头/质控
  - fastqc_trimmed：修剪后质控
  - star/hisat2：比对
  - 后续分析（如 PureCLIP、UmiTools 等可扩展）
- **输入**：原始fastq，配置json
- **输出**：标准bam、质控报告等

igv模块准备:

```json
  "igv": {
        "js": "/data/pub/zhousha/Reference/igv.min.js",
        "id": "mm39",
        "name": "Mouse (GRCm39/mm39)",
        "publicPathMap": {
            "/data/pub/zhousha/": "/data/",
            "/data/pub/zhousha/Reference/": "/ref/"
        },
        "fastaURL": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa",
        "indexURL": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa.fai",
        "cytobandURL": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/mm39.cytoBand.txt",
        "tracks": [
            {
                "name": "Gencode vM38 genes",
                "format": "gtf",
                "url": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.basic.gene_exon.sorted.gtf.gz",
                "indexUrl": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.basic.gene_exon.sorted.gtf.gz.tbi"
            },
            {
                "name": "Gencode rmsk repeats",
                "format": "gtf",
                "url": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.sorted.gtf.gz",
                "indexUrl": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.sorted.gtf.gz.tbi"
            }
        ]
    }
```

1. 注释轨道

- 注释轨道需要自己建立索引，以防止浏览器全量加载注释卡死（F12观察返回码是否为200,200则是没有配置索引）
```sh
tabix -g gff [gft.gz]
```
- nginx需要确保sendfile是开启状态(on)

2. publicPathMap

为了让nginx能够获取流程生成html内的内容，特意配置了这个map，以生成相对路径让nginx能够读取

待做：
- [] 增加UMI提取方式字段

#### UMI-tools工作原理

Read: TAGCCGGCTTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTAG Barcode: NNNXXXXNN

Barcode: TAGCCGGCT
UMI: TAGCT
library barcode: CCGG
Processed data: CCGGTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTAG



#### UMI提取和trim顺序

UMI提取依赖序列不被破坏，建议先提取UMI，再做trim比较安全

#### 常见UMI模式

1. NNNXXXXNN

常见iCLIP实验

2. NNNNNXXXXXXNNNN

[iCLIP2 protocol](https://doi.org/10.1016/j.ymeth.2019.10.003)

解释：

- N为随机序列，即UMI
- XXXX为barcode，用于区分样本，后续demultiplex拆分样本(如果有多个样本的话)

注意：

如果是单个样本，X需要被替换为N，即全是UMI


### 2. CoCulture.smk
- **用途**：共培养体系转录组分析。
- **主要模块**：
  - SOAPnuke：原始数据过滤
  - hisat2：多物种比对
  - ngs_disambiguate：去除混合比对
  - TEtranscripts（TEcount + TElocal）：基因和转座子表达定量
- **输入**：fastq，物种基因组信息
- **输出**：区分物种的 BAM、TE 表达量矩阵

### 3. MERIP.smk
- **用途**：MeRIP-seq/m6A-seq 数据分析。
- **主要模块**：
  - cutadapt：去接头
  - hisat2：比对
  - igv：可视化
  - exomePeak：甲基化位点检测
- **输入**：fastq，设计信息
- **输出**：dedup bam、peak表等

### 4. RNAseq.smk

- **用途**：常规转录组综合分析（表达定量、差异分析、融合基因、变异检测）。

#### 管道流程

```
FASTQ
  ├─ cutadapt / trimmomatic / trim_galore ── trimmed FASTQ
  │    ├─ star / hisat2 ── BAM ── TEtranscripts (TEcount + TElocal) ── TE 矩阵
  │    │    └─ DESeq2 ── 差异表达
  │    │         └─ function ── GO/KEGG/GSEA
  │    ├─ hisat2 (dta) ── StringTie ── 转录本组装 + TE chimeric
  │    ├─ RmrRNA + bowtie2 ── rRNA去除 ── STAR (chimeric) ── gatk_prepare ── arriba
  │    ├─ STAR (SNP) ── XenofilteR ── gatk_prepare_SNP ── gatk_RNAseq ── RNA SNP VCF
  │    └─ RNAseq_report ── 汇总报告
  └─
```

#### 模块控制 (enabled)

每个模块通过 `Params.<module>.enabled` 控制是否产出最终文件。默认 `true`，向后兼容。

下游模块拥有最高权限——Snakemake DAG 自动解析依赖，下游 enabled 会自动拉起上游产出。

| Params key | 产出文件 | 说明 |
| --- | --- | --- |
| `Params.TEtranscripts.enabled` | `counts/{genome}/TEcount/all_TEcount.tsv`<br>`counts/{genome}/TElocal/all_TElocal.tsv` | TE 表达定量矩阵 |
| `Params.DESeq2.enabled` | `diff_expression/{genome}/{ctr}_vs_{exp}/DESeq2.done` | 差异表达分析；需要 TEcount 作为输入 |
| `Params.function.<organism>.enabled` | `function/{genome}/{pair}/go_*.png`<br>`function/{genome}/{pair}/kegg_*.png`<br>`function/{genome}/{pair}/GSEA/*.jpeg`<br>`function/{genome}/{pair}/*.csv` | GO/KEGG/GSEA，per-organism 控制 |
| `Params.StringTie.enabled` | `transcripts/{genome}/stringtie_merged.gtf`<br>`transcripts/{genome}/TE_chimeric/*.png`<br>`transcripts/{genome}/TE_chimeric/*.tsv` | 转录本组装 + TE chimeric 分析 |
| `Params.arriba.enabled` | `fusion/{genome}/arriba_report/arriba_fusion_report.html`<br>`fusion/{genome}/{sample}/{sample}_passed_fusions.tsv` | 融合基因检测 |
| `Params.gatk_RNAseq.enabled` | `variation/germline_snv_indel_RNAseq/{genome}/{sample}/{sample}.filtered.vcf.gz` | RNA 编辑 / 胚系 SNV/INDEL |
| `Params.report.enabled` | `results/{genome}/RNAseq_report.pptx` | 汇总报告 |

#### 配置示例

```json
{
    "Params": {
        "TEtranscripts": {"enabled": true},
        "DESeq2": {"enabled": true},
        "StringTie": {"enabled": true},
        "arriba": {"enabled": true},
        "gatk_RNAseq": {"enabled": true},
        "report": {"enabled": true},
        "function": {
            "GRCm39": {"enabled": true, "species": "mouse"},
            "GRCh38": {"enabled": true, "species": "human"}
        }
    }
}
```

**输入**：fastq，meta_input.tsv（sample_id / organism / group / data_id / layout / fastq_1 / fastq_2 / contaminated_organism）

**输出**：TE 表达矩阵、差异分析、融合基因、RNA SNP、功能富集、汇总报告


### 5. ncRNAseq.smk
- **用途**：非编码 RNA 分析流程（含 tailer 和小 RNA）。
- **主要模块**：
  - fastqc：原始数据和修剪后质控
  - demultiplexer：barcode 解复用、UMI 提取、去重
  - trim_galore：去接头
  - hisat2/star：比对
  - gatk_prepare：addReadsGroup + MarkDuplicates
  - star small RNA：小 RNA 提取 + 3pass 比对
  - tailer：3' tail 分析
- **输入**：fastq
- **输出**：ncRNA 表达量、tailer 分析结果

### 6. PeakCalling.smk
- **用途**：ChIP-seq / ChIRP-seq / DIP-seq peak calling 分析。
- **主要模块**：
  - fastqc：原始数据和修剪后质控
  - trim_galore：去接头
  - bowtie2：建立索引和比对
  - gatk_prepare：addReadsGroup + MarkDuplicates
  - igv：dedup BAM → bigWig 可视化
  - track：IGV/UCSC track 生成
  - macs3：peak calling + cutoff 图
  - frip_score：FRiP 质控
  - deeptools_heatmap：bigWig ratio + 热图
  - homer：annotatepeaks 峰注释
  - peak_te_overlap：peak-TE overlap 分析 + 可视化
  - PeakCalling_report：汇总报告
- **输入**：fastq，IP 样本和 Input 对照样本设计信息
- **输出**：比对 BAM、narrowPeak peak 文件、FRiP 分数、deeptools 热图、Homer 注释、peak-TE overlap、HTML 报告
- **特点**：
  - 支持有对照（IP vs Input）和无对照模式
  - 支持多个基因组
  - 适用于转录因子结合位点（ChIP-seq）和 DNA 修饰（DIP-seq）分析

### 7. Mutation.smk
- **用途**：体细胞突变分析（tumor vs normal）及cfDNA片段分析。
- **主要模块**：
  - fastqc：质控
  - cutadapt：去接头
  - bwa-mem2：比对
  - gatk：BQSR、Mutect2 体细胞突变检测、胚系突变检测
  - spectrum：突变频谱可视化
  - fragment_size：cfDNA片段长度分析（可选）
  - manta：结构变异检测（可选）
  - cnvkit：拷贝数变异检测（可选）
- **输入**：fastq，design pairs（tumor/normal配对信息）
- **输出**：体细胞突变 VCF、突变频谱图、片段长度分布图、SV VCF、CNV结果
- **特点**：
  - 支持跳过片段长度分析（`Params.skip_fragment_size=true`）
  - 支持跳过SV检测（`Params.skip_sv=true`）
  - 支持跳过CNV检测（`Params.skip_cnv=true`）
  - CNVkit支持对照样本构建参考（`control_samples`参数）
  - 片段长度分析适用于cfDNA/ctDNA液体活检样本

### 8. PacVar.smk
- **用途**：PacBio 长读长变异检测。
- **主要模块**：
  - pbmm2：PacBio 专用比对
  - deepvariant/gatk4：SNP/INDEL 检测
  - pbsv：结构变异检测
  - hiphase：单倍型 phasing
  - trgt：重复序列分析
  - telogator2：端粒长度分析（per-chromosome-arm）
  - telomere assembly scan：基于 assembly 的端粒扫描（推荐用于小鼠）
  - telomere read density：基于 read 的端粒 k-mer 密度估算
  - tidk：社区工具端粒扫描
  - centromere：着丝粒分析（hifiasm + RepeatMasker）
- **输入**：PacBio BAM/fastq
- **输出**：SNP VCF、SV VCF、phasing 结果、重复序列分析结果、端粒长度、着丝粒统计
- **特点**：
  - 支持跳过特定步骤（skip_snp/skip_sv/skip_phase/skip_repeat）
  - 支持多种 SNV caller（deepvariant/gatk4）
  - 支持端粒分析（4 种方法）和着丝粒分析
  - 端粒方法：telogator2（TL_p75）、assembly scan（推荐小鼠）、read density（全基因组平均）、tidk
- **问题**：
  - 如何从PacBio数据识别端粒长度超过hifi读长的物种或品种（比如小鼠，hifi平均reads长度大约在15kb，适合测量人类的端粒长度），流程目前所集成的几大方法都有很大缺陷。确定着丝粒长度也非常有难度。除非组装出T2T基因组

### 9. KARRseq.smk
- **用途**：Kethoxal-Assisted RNA-RNA interaction sequencing 分析。
- **主要模块**：
  - STAR：比对
  - 自定义脚本：提取 chimeric reads、去除重复、生成 ligation pairs
- **输入**：fastq
- **输出**：去重后的 ligation pairs 文件
- **特点**：
  - 用于研究 RNA-RNA 相互作用
  - 需要自定义 STAR 参数输出 chimeric reads

### 10. QuantMS.smk
- **用途**：定量蛋白质组学分析（TMT/LFQ/DIA）。
- **主要模块**：
  - DecoyDatabase：生成诱饵数据库
  - CometAdapter/MSGFPlusAdapter/SageAdapter：数据库搜索引擎
  - PercolatorAdapter：PSM 重评分
  - FalseDiscoveryRate：PSM FDR 控制
  - Epifany：蛋白质推断
  - ProteomicsLFQ/ProteinQuantifier：蛋白质定量
  - MSstatsConverter：统计分析
- **输入**：mzML 文件、蛋白质数据库（FASTA）
- **输出**：mzTab 定量结果、MSstats 统计分析结果
- **特点**：
  - 支持 TMT、LFQ、DIA 三种定量方法
  - 支持多种搜索引擎（Comet、MSGF+、Sage）
  - 支持跳过 MSstats 分析（skip_post_msstats）
  - 基于 OpenMS 工具集

### 11. tRNAseq.smk
- **用途**：tRNA 修饰诱导错配测序分析（mim-tRNAseq）。
- **模块架构**：6 个独立子模块，通过 pickle 状态文件传递中间数据
  - `tRNAtools`：tRNA 注释、聚类、SNP 索引构建 → `state/`
  - `align`：GSNAP 比对（SNP-tolerant，支持两轮 remap）→ `samples/{sample}/*.bam`
  - `clusters`：isodecoder 去卷积（簇拆分为唯一定位转录本）
  - `mods`：错配 / 修饰定量（min_cov + misinc_thresh 过滤）
  - `coverage`：覆盖度计算 + CCA 分析 → `coverage_byaa.txt`、`coverage_bygene.txt`
  - `deseq`：DESeq2 差异表达（需指定 control_cond）
- **输入**：trimmed FASTQ（sample data sheet, TSV 格式）、物种 tRNA 参考
- **输出**：覆盖度报告、修饰表格、CCA 统计、DESeq2 结果、可视化图
- **特点**：
  - 基于 mim-tRNAseq 工具（mimseq Python 库）
  - 支持内置物种（Hsap, Mmus, Scer 等）和自定义 tRNA 参考
  - 所有样本一起处理，通过 sample data sheet 驱动
  - 子模块详细文档见 `modules/mimseq/README.md`

### 12. Fiberseq.smk
- **用途**：Fiber-seq 表观遗传分析（PacBio Fiber-seq 数据）。
- **主要模块**：
  - fibertools_m6a：ft predict m6A 甲基化检测
  - fibertools_nuc：ft add nucleosomes 核小体定位
  - fibertools_fire：ft fire 染色质 FIRE 活性区域检测
  - fibertools_extract：ft extract 数据提取
- **输入**：PacBio HiFi BAM
- **输出**：m6A 甲基化 BAM、核小体 BAM、FIRE 区域、提取结果

### 13. Population_genomics.smk
- **用途**：群体基因组学分析。
- **主要模块**：
  - population_metadata：样本元数据生成
  - gatk_population：HaplotypeCaller 单样本变异检测 + joint genotyping
  - bcftools_population：变异过滤
  - plink2_population：群体统计分析
  - vcftools_population：VCF 统计
  - admixture：群体结构推断（可选）
  - easySFS：群体历史推断（可选）
  - PopLDdecay：LD 衰减分析（可选）
- **输入**：BAM 文件、样本元数据
- **输出**：joint-called VCF、过滤后 VCF、plink2 统计、群体结构、SFS、LD 衰减

## 14. scRNAseq.smk

- **用途**：标准化 Scanpy 单细胞 RNA-seq 分析。
- **输入**：由 `tissue_samples` 按组织分组，每个样本的 h5ad 路径由 `indir + sample_id` 拼接。
- **基础分析**：QC、线粒体比例过滤、归一化、log1p、高变基因、PCA、neighbors、UMAP、Leiden、marker 和差异表达。
- **AI 自动注释**：LLM-assisted 细胞类型自动注释（详见下方）。
- **高级分析**：可配置启用 DPT 拟时序、RNA velocity、LIANA 细胞通讯和 infercnvpy CNV。
- **架构**：规则位于 `modules/scanpy/scanpy.smk`，本 subworkflow 只负责配置和编排。

### aligner / counter 组合

| aligner | counter | 说明 |
| --- | --- | --- |
| star | scTE | STARsolo 比对 + scTE TE 定量 |
| cellranger | scTE | Cell Ranger 比对 + scTE TE 定量（推荐） |
| cellranger | cellranger | Cell Ranger 比对 + Cell Ranger 基因定量 |

不存在 `aligner: star, counter: star`，因为 Cell Ranger 本身就是对 STAR 的修改，且只接受 FASTQ 输入。
选择了 `cellranger` + `scTE` 时，虽然 Cell Ranger 自身也会做基因定量，但 scTE 的 TE 定量是独立的。

### AI 自动聚类注释

聚类后的细胞类型注释是最依赖经验的环节。CellTypist 等工具只覆盖人和小鼠，非模式物种（恒河猴、猪等）只能靠人工逐个查 marker 文献。AI 注释的目标是在聚类完成后秒级给出有据可查的初始注释，把"从聚类到初步结果"的时间从几小时压缩到几分钟。

#### 整体架构

```
聚类完成 (leiden)
      │
      ▼
┌─────────────────────────────┐
│  逐 cluster 注释循环         │
│  ┌───────────────────────┐  │
│  │ 1. 提取 cluster 统计   │  │
│  │    - top 50 DEG       │  │
│  │    - 细胞数/基因数/UMI │  │
│  │    - MT%              │  │
│  └──────────┬────────────┘  │
│             ▼               │
│  ┌───────────────────────┐  │
│  │ 2. 构建结构化 prompt   │  │
│  │    - 组织类型上下文    │  │
│  │    - 质量标记规则      │  │
│  │    - JSON 输出约束     │  │
│  └──────────┬────────────┘  │
│             ▼               │
│  ┌───────────────────────┐  │
│  │ 3. LLM 推断            │  │
│  │    - OpenAI API 兼容   │  │
│  │    - 或本地 Ollama     │  │
│  │    - temperature=0.1   │  │
│  │    - json_object 格式  │  │
│  └──────────┬────────────┘  │
│             ▼               │
│  ┌───────────────────────┐  │
│  │ 4. PubMed 文献验证     │  │
│  │    - 为 key_markers    │  │
│  │      检索 PMIDs        │  │
│  │    - 限速 0.35s/req    │  │
│  └──────────┬────────────┘  │
│             ▼               │
│  ┌───────────────────────┐  │
│  │ 5. 程序化质量校验      │  │
│  │    - 与 AI flag 交叉   │  │
│  │    - 标记待过滤 cluster │  │
│  └───────────────────────┘  │
└─────────────────────────────┘
      │
      ▼
  写入 adata.obs["cell_type"]
  写入 adata.uns["annotation"]
```

#### Prompt 注入设计

每个 cluster 的 prompt 包含两层信息：

**第一层：统计上下文**——把 cluster 的数值特征直接注入 prompt，让 LLM 判断质量：

```
## Cluster 3
- Cells: 1523
- Mean genes/cell: 2100
- Mean UMI/cell: 8500
- MT%: 3.2%
- Top 50 DEGs (ranked by Wilcoxon): COL1A1, COL3A1, DCN, LUM, ...
```

这些数字来自 `adata.obs` 的实时计算（`n_genes_by_counts`、`total_counts`、`pct_counts_mt`），不是配置里的静态值。

**第二层：质量标记规则**——硬编码在 prompt 末尾，作为 LLM 的判断依据：

| 规则 | 触发条件 | 输出 |
| --- | --- | --- |
| 未注释基因 | top DEG 以 ENSMMUG 开头 | `cell_type="Unannotated"`, `quality_flag="unannotated"` |
| MT 污染 | top DEG 以 MT-/mt- 开头 | `quality_flag="low_quality"` |
| 核糖体污染 | top DEG 为 RPS/RPL | `quality_flag="ribosomal"` |
| TE 富集 | top DEG 为 Alu/ERVK/LINE/SINE/LTR | `cell_type="ERVK_high"`, `quality_flag="te_dominated"` |
| 无法判断 | 以上都不匹配且无明确 marker | `cell_type="Unknown"`, `quality_flag="unknown"` |

这些规则写进 prompt 而不是后处理，是因为 LLM 能结合基因名称和生物学意义做更灵活的判断（比如同时有少量 MT 基因和明确的上皮 marker 时，不会误判为低质量）。

**输出格式约束**：`response_format={"type": "json_object"}` + prompt 末尾 "Output ONLY valid JSON"，强制 LLM 返回结构化结果，避免解析失败。

#### 错误兜底机制

**API 调用层**：`_call_openai()` 和 `_call_ollama()` 都包裹在 try/except 中。调用失败时返回空 dict `{}`，不会中断整个流程。后续代码通过 `annotation.get("cell_type", "Unknown")` 兜底，失败的 cluster 标记为 Unknown。

**JSON 解析层**：依赖 `response_format={"type": "json_object"}`（OpenAI）和 `"format": "json"`（Ollama）确保返回合法 JSON。如果解析仍然失败，except 捕获后返回空 dict。

**集群级兜底**：`_analyze_cluster_quality()` 在 AI 注释之后做程序化二次校验——独立计算 mean_genes、mean_counts、pct_mt，与 AI 返回的 `quality_flag` 交叉验证。即使 AI 漏标，程序化检查也能兜住。

**PubMed 检索兜底**：`_search_pubmed()` 失败时返回空列表，不影响注释结果，只丢失文献引用。

**逐 cluster 隔离**：每个 cluster 独立调用 LLM，单个 cluster 的 API 失败不会影响其他 cluster。cluster 之间有 `time.sleep(0.5)` 限速，避免触发 API 频率限制。

#### 上下文维护

注释结果写入两个位置：

- `adata.obs["cell_type"]`：每个 cell 的细胞类型标签，用于下游可视化和差异分析
- `adata.uns["annotation"]`：完整的注释详情，结构为：

```python
{
    "3": {
        "cell_type": "Stromal",
        "key_markers": ["COL1A1", "COL3A1", "DCN", "LUM"],
        "reasoning": "High expression of extracellular matrix genes...",
        "confidence": "high",
        "quality_flag": null,
        "references": {
            "COL1A1": [{"pmid": "12345678", "title": "...", "year": "2020"}],
            "DCN": [{"pmid": "87654321", "title": "...", "year": "2019"}]
        }
    }
}
```

`references` 字段将 PubMed 检索结果嵌入注释，后续人工审核时可以直接点开 PMID 查阅原文，不需要再单独查文献。

#### LLM 后端配置 (mode auto)

```json
{
    "llm_method": "anthropic",
    "llm_model": "MiniMax-M3",
    "llm_api_key": "",
    "llm_base_url": "https://api.minimax.cn/anthropic"
}
```

- `llm_method`：`"anthropic"`（anthropic 协议）/ `"openai"`（OpenAI API 兼容）/ `"ollama"`（本地部署）
- `llm_api_key` 和 `llm_base_url` 支持通过环境变量 `LLM_API_KEY` 和 `LLM_BASE_URL` 注入，config 中留空时自动读取环境变量
- 上述字段属于 `mode auto`，配置在 `Params.scanpy.{counter}.auto.*` 下；`mode annotate` 不再支持 LLM 字段（已合并到 `mode auto`）
- 官方推荐 base_url: `https://api.minimax.cn/anthropic`


示例：

```json
{
  "input_h5ad": "/path/input.h5ad",
  "Params": {
    "scanpy": {
      "advanced": {
        "trajectory": true,
        "velocity": false,
        "communication": false,
        "cnv": false
      }
    }
  }
}
```

### 15. spatial_transcriptomics.smk

- **用途**：10x Visium 和自定义 spot×gene h5ad 空间转录组分析。
- **Visium 输入**：`visium_h5` 和可选 `spatial_dir`。
- **自定义输入**：`input_h5ad`，必须包含 `adata.obsm["spatial"]` 空间坐标。
- **基础分析**：空间 QC、归一化、高变基因、PCA、空间聚类和空间 marker。
- **高级分析**：Squidpy 空间邻域和 Moran's I 空间自相关。
- **架构**：规则位于 `modules/spatial_scanpy/spatial_scanpy.smk`，本 subworkflow 只负责配置和编排。

cell2location 依赖已在环境模板中预留，但只有在明确参考 scRNA 表达矩阵和细胞类型 signature 输入后才应启用。

Visium 输入示例：

```json
{
  "visium_h5": "/path/filtered_feature_bc_matrix.h5",
  "spatial_dir": "/path/spatial"
}
```

自定义 h5ad 输入示例：

```json
{
  "input_h5ad": "/path/spatial.h5ad"
}
```

高级分析示例：

```json
{
  "visium_h5": "/path/filtered_feature_bc_matrix.h5",
  "spatial_dir": "/path/spatial",
  "Params": {
    "spatial": {
      "advanced": {
        "spatial_autocorrelation": true
      }
    }
  }
}
```

## 使用说明

- 每个 .smk 文件可作为 Snakemake 主入口，需配合对应 config json/yaml。
- 支持 conda 环境自动管理。
- 具体参数和模块细节见各模块目录及主 config。
