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

### 6. TEtranscipts.smk
- **用途**：转座子表达分析专用流程。
- **主要模块**：
  - TEtranscripts
- **输入**：bam
- **输出**：TE表达量

### 7. featureCounts.smk
- **用途**：基因/转座子计数。
- **主要模块**：
  - featureCounts
- **输入**：bam
- **输出**：count矩阵

### 8. ncRNAseq.smk
- **用途**：非编码RNA分析流程。
- **主要模块**：
  - cutadapt
  - hisat2/star
  - featureCounts
- **输入**：fastq
- **输出**：ncRNA表达量

---

## 使用说明
- 每个 .smk 文件可作为 Snakemake 主入口，需配合对应 config json/yaml。
- 支持 conda 环境自动管理。
- 具体参数和模块细节见各模块目录及主 config。
