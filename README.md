# Omics

记录各种组学分析流程，实践最佳流程管理


## 快速开始

1. 创建envs中的conda环境，对于各个subworflow，配置好同名yaml内容

注意：
- trim_galore只是打包程序，需要确保cutadapt存在
- 物种名是动态变化的，在yaml中，相关文件的键需要根据metadata改变，如果物种名为有空格隔开，需要改空格为_，如Mus musculus -> Mus_musculus

2. 修改scripts下RNA-SNP_prepare.sh脚本，执行

3. 运行

--config参数可以在main.json对应字段中填写

```sh
snakemake -s workflow/RNA-SNP/main.smk --config indir=data/fq outdir=output metadata=data/target_fq.tsv --cores 45 --use-conda --conda-prefix /path/to/enviroment

```
好用参数：--rerun-triggers input

- --conda-prefix 制定conda包下载地址，比如：/data/pub/zhousha/env/mutation_0.1

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

## run.py

读取模板json，调用指定分析流程:
  - 解析meta信息，覆盖相关字段
  - 模板json默认配置保留
  - 支持命令行参数更新指定字段

## modules

见各模块下json文件
问题：
- 哪些字段可以为空（初步想法是，标注null的字段可以为空，不标注的不可以不传递）
- 字段值约束如何做（初步想法是，增加字段名_comment作为字段约束，后续可以在流程内加约束条件）
- 模块接口（初步想法是通过indir,outdir；拿indir举个例子，有内层结构就input/内层路径，无内层结构就input）
- 是否需要给文件输入输出加后缀，还是同样的文件格式都是sample_id为文件名

## subworkflow

无

## 参数详解

- `--config`
  - `indir=value`：指定 fastq.gz 输入路径
  - `outdir=value`：指定输出文件目录
  - `metadata=value`: 指定fastaq元信息，必须包含["Data_id","Sample_id","Organism"],详见utils/fastq_utils.py如何处理metadata


## 设计规范

规则尽量不硬编码路径，只接受indir, outdir和logdir，还有生成的wildcards
规则应该尽可能少涉及与执行无关的信息

1. 规则可复用，灵活性高
2. 可指定分析终点
3. 具备强大的元信息处理
4. 尽量在流程开始前规范化输入信息，流程本身不承担规范责任


## 待做
- [x] 元信息控制CoCulture流程，记得完善json
- [x] RNAseq流程重写
- [ ] 多模块并行，（指定.snakemake生成与不同位置，或者拼接规则）

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
Rscript workflow/scripts/function/DESeq2.r --mode TEcount \
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