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
  `_R1*fastq` , `_R2*fastq` or `_R1*fq` , `_R2*fq` or `_1*fastq` , `_2*fastq` or `_1*fq` , `_2*fq`
- 单端测序文件：非双端测序文件模式
- 具体查看utils/fastq_utils.py

## subworkflow像解

### Align

1. **流程控制变量**
- single_samples: 单端样本id列表
- paired_samples: 双端样本id列表
- genomes: 物种列表

### Annovar

### StringTie

### TEtranscripts

1. **流程控制变量**
- all_samples: 样本id列表
- genomes: 物种列表

### XenofilterR

1. **流程控制变量**
- XenofilterR_target_samples: 代表可能含有其它物种序列污染的样本列表，在示例流程中为用小鼠饲养层细胞培育的干细胞
- XenofilterR_pollution_source_genome: 代表污染来源基因组字符串，在示例流程中为小鼠

解释：
- 各subworkflow变量需要在run.smk中定义，各subworkflow同名变量是同一个变量，均来自run.smk

## 参数详解

- `--config`
  - `indir=value`：指定 fastq.gz 输入路径
  - `outdir=value`：指定输出文件目录
  - `metadata=value`: 指定fastaq元信息，必须包含["data_id","sample_id","organism"],详见utils/fastq_utils.py如何处理metadata


## 流程脆弱之处
1. XenofilterR与SNP衔接
为了支持去除bam中的污染序列，采用了XenofilterR,这个R脚本输入是个文件，没有提供命令行参数。需要使用临时目录作为input，同时对addReadsGroup规则
的params区分了污染物种和非污染物种的bam，以正确生成vcf。之所以说是脆弱，要想跑通它，需要理解并修改这两个规则。其它只要按照yaml填空就行
