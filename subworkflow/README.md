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

待做：
- [] 增加UMI提取方式字段

### UMI-tools工作原理

Read: TAGCCGGCTTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTAG Barcode: NNNXXXXNN

Barcode: TAGCCGGCT
UMI: TAGCT
library barcode: CCGG
Processed data: CCGGTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTAG



### UMI提取和trim顺序

UMI提取依赖序列不被破坏，建议先提取UMI，再做trim比较安全

### 常见UMI模式

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
- **输入**：fastq，物种基因组信息
- **输出**：区分物种的bam、表达量等

### 3. MERIP.smk
- **用途**：MeRIP-seq/m6A-seq 数据分析。
- **主要模块**：
  - cutadapt：去接头
  - hisat2：比对
  - igv：可视化
  - exomePeak：甲基化位点检测
- **输入**：fastq，设计信息
- **输出**：dedup bam、peak表等

### 4. RNA_SNP.smk
- **用途**：RNA变异检测流程。
- **主要模块**：
  - cutadapt：去接头
  - hisat2/star：比对
  - XenofilterR/gatk：变异检测
- **输入**：fastq
- **输出**：SNP/INDEL结果

### 5. RNAseq.smk
- **用途**：常规转录组分析。
- **主要模块**：
  - cutadapt：去接头
  - hisat2/star：比对
  - TEtranscripts：转座子表达
- **输入**：fastq
- **输出**：表达量矩阵


### 6. ncRNAseq.smk
- **用途**：非编码RNA分析流程。
- **主要模块**：
  - cutadapt
  - hisat2/star
  - featureCounts
- **输入**：fastq
- **输出**：ncRNA表达量

### 7. PeakCalling.smk
- **用途**：ChIP-seq / ChIRP-seq / DIP-seq peak calling 分析。
- **主要模块**：
  - cutadapt：去接头/质控
  - bowtie2：建立索引和比对
  - macs3：peak calling
- **输入**：fastq，IP样本和Input对照样本设计信息
- **输出**：比对BAM、narrowPeak peak文件
- **特点**：
  - 支持有对照（IP vs Input）和无对照模式
  - 支持多个基因组
  - 适用于转录因子结合位点（ChIP-seq）和DNA修饰（DIP-seq）分析

### 8. Mutation.smk
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

### 9. PacVar.smk
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

### 10. KARRseq.smk
- **用途**：Kethoxal-Assisted RNA-RNA interaction sequencing 分析。
- **主要模块**：
  - STAR：比对
  - 自定义脚本：提取 chimeric reads、去除重复、生成 ligation pairs
- **输入**：fastq
- **输出**：去重后的 ligation pairs 文件
- **特点**：
  - 用于研究 RNA-RNA 相互作用
  - 需要自定义 STAR 参数输出 chimeric reads

### 11. QuantMS.smk
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

### 12. tRNAseq.smk
- **用途**：tRNA 修饰诱导错配测序分析（mim-tRNAseq）。
- **主要步骤**：
  - tRNA 聚类 + SNP 索引构建
  - GSNAP 比对（SNP-tolerant）
  - 簇拆分 / 反卷积为唯一 tRNA 转录本
  - 覆盖度计算 + QC 图
  - 错配 / 修饰定量
  - 3'-CCA 分析（tRNA 完整性）
  - DESeq2 差异表达
  - （可选）SLAC 修饰-氨基酸酰化串扰分析
- **输入**：trimmed FASTQ（sample data sheet, TSV 格式）、物种 tRNA 参考
- **输出**：覆盖度报告、修饰表格、CCA 统计、DESeq2 结果、可视化图
- **特点**：
  - 基于 mim-tRNAseq 一体化工具（mimseq 命令）
  - 支持内置物种（Hsap, Mmus, Scer 等）和自定义 tRNA 参考
  - 所有样本一起处理，通过 sample data sheet 驱动

---

## 使用说明
- 每个 .smk 文件可作为 Snakemake 主入口，需配合对应 config json/yaml。
- 支持 conda 环境自动管理。
- 具体参数和模块细节见各模块目录及主 config。
