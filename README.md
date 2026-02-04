# RNASNP

## 快速开始

1. 创建envs中的conda环境，对于各个subworflow，配置好同名yaml内容

注意：
- trim_galore只是打包程序，需要确保cutadapt存在
- 物种名是动态变化的，在yaml中，相关文件的键需要根据metadata改变，如果物种名为有空格隔开，需要改空格为_，如Mus musculus -> Mus_musculus

2. 修改scripts下RNA-SNP_prepare.sh脚本，执行

3. 运行

```sh
snakemake -s workflow/RNA-SNP/run.smk --config indir=data/fq outdir=output metadata=data/target_fq.tsv --cores 45 --use-conda
```
4. 亮点

- 支持多个物种
- 支持单端测序和双端测序文件
- 统一的日志管理

注意：

- 双端测序文件：
  `_R1*fastq` , `_R2*fastq` or `_R1*fq` , `_R2*fq`
  `_1*fastq` , `_2*fastq` or `_1*fq` , `_2*fq`
- 单端测序文件：非双端测序文件模式
- 具体查看utils/fastq_utils.py

## subworkflow详解

### Align

1. **流程控制变量**
- single_samples: 单端样本id列表
- paired_samples: 双端样本id列表
- single_sample_genome_pairs: 单端样本id，物种配对列表，如：[(SE_sample_id， organsim)……]
- paired_sample_genome_pairs: 双端样本id，物种配对列表，如：[(SE_sample_id， organsim)……]

### Annovar

### SNP
- XenofilterR_target_genome: 代表可能含有其它物种序列污染的样本列表，在示例流程中为用小鼠饲养层细胞培育的人类干细胞

### StringTie

### TEtranscripts

1. **流程控制变量**

- single_sample_genome_pairs: 单端样本id，物种配对列表，如：[(SE_sample_id， organsim)……]
- paired_sample_genome_pairs: 双端样本id，物种配对列表，如：[(SE_sample_id， organsim)……]

### XenofilterR

1. **流程控制变量**
- XenofilterR_target_samples: 代表可能含有其它物种序列污染的样本列表，在示例流程中为用小鼠饲养层细胞培育的人类干细胞
- XenofilterR_target_genome: 代表目标物种，在示例流程中为人类
- XenofilterR_pollution_source_genome: 代表污染来源基因组字符串，在示例流程中为小鼠

解释：
- 各subworkflow变量需要在run.smk中定义，各subworkflow同名变量是同一个变量，均来自run.smk

## 参数详解

- `--config`
  - `indir=value`：指定 fastq.gz 输入路径
  - `outdir=value`：指定输出文件目录
  - `metadata=value`: 指定fastaq元信息，必须包含["Data_id","Sample_id","Organism"],详见utils/fastq_utils.py如何处理metadata


## 待改进之处
- 目前只支持一对污染物种运行，可以针对XenofilterR和SNP规则进行修改使其支持

## 设计规范

1. capability 以 最终结果 为中心设计; 每个 capability 只有一个 facade rule

2. 可替换实现用 include + config dispatch

3. schema / defaults 只声明，不参与运行

4. 文件命名通过 contract 函数统一管理

5. work/ ≠ results/，下游不碰中间态

6. 下游通过 rules.xxx.input 耦合能力，而不是路径

## 下游分析

### annotation

  - gene_id2name: 映射基因id到基因名称
  - geneIDAnnotation: 注释基因id
  - vcf_annovar: annovar注释vcf

### count

  - count: RNAseq count标准化
  - normalization: r版本，待完善

### download

  - ascp_download: 通过ascp下载对应 sra id(run number)代表的fastq序列
  - GSE_runinfo: 获取GSE对应的run信息
  - GSM_metadata: 获取GSM对应信息
  - GSM_resolver: 解析多种输入的  GSM id
  - old/*: 一些旧有脚本，已通过python重构

### function

  - DESeq2: 批量差异分析
  - gmt.py: GSEA富集分析
  - go-kegg_back: GO,KEGG富集分析，绘制背对背柱状图
  - go-kegg: GO,KEGG富集分析，绘制散点图
  - gsea: GSEA富集分析
  - TEsite_subfamily: 统计差异分析中TE subfamily情况

#### DESeq2

举例

```sh
Rscript workflow/scripts/DESeq2/DESeq2.r --mode TEcount \
    --matrix output/${outdir}/counts/mouseTEcount.cntTable \
    --group ${group} \
    --pattern  ${control} ${experiment} \
    --outdir output/${outdir}/ \
    --annotation /ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv \
    --figure pca heatmap volcano \
    --TEcountMode all
```


### fusion

  - fusionGeneAnnotation: 融合基因名称注释
  - volcanoImportant: 融合基因火山图

### gatk

  - gatkPrepare: gatk上游分析shell脚本

### map

  - star: star比对shell脚本
  - starPrepare: star建立索引shell脚本

### plot

  - heatmap: 热图
  - pie: 饼状图
  - violin: 小提琴图 

### SNP

  - bash/*: 一些shell脚本
  - count/*: count分析脚本
  - python/*: 一些老旧Python脚本
  - QC/*: 质量控制
  - Rscript/*: 一些老旧R脚本
  - Terr_RNA/*: 端粒转录序列分析
  - vcf/*: vcf分析脚本

### SV
  - pbsv_sv_diff_analysis: 比较对照组和实验组结构变异类型差异情况
  - PlaB_only: 比较实验组特有alt序列与对照组对应alt序列注释情况差异，绘制富集图，以及一些统计
  - run_circos: 绘制vcf的结构变异circos图，注意vcf文件寻找模式替换
  - utils/repeatmasker_analysis: 注释序列重复元件情况，比较不同序列注释结果
  - utils/repeatmasker_plot: 绘制重复序列注释结果差异富集柱状图
  - utils/SV_TYPE_plot: 绘制结构变异类型图
  - utils/SV_TYPE: 统计vcf中结构变异情况

### train

潜在因子分析