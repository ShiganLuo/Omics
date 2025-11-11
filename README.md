# RNASNP

## 快速开始

1. 创建envs中的conda环境，对于各个subworflow，配置好同名yaml内容

注意：trim_galore只是打包程序，需要确保cutadapt存在

2. 运行

```sh
snakemake -s workflow/RNA-SNP/run.smk --config indir=data/fq outdir=output metadata=data/target_fq.tsv--cores 45 --use-conda
```
3. 亮点

- 支持多个物种
- 支持单端测序和双端测序文件
- 统一的日志管理

注意：

- 双端测序文件：
  _R1*fastq , _R2*fastq or _R1*fq , _R2*fq
  _1*fastq , _2*fastq or _1*fq , _2*fq
- 单端测序文件：非双端测序文件
- 建议用函数控制输入文件

## 参数详解

- `--config`
  - `indir=value`：指定 fastq.gz 输入路径
  - `outdir=value`：指定输出文件目录
  - `metadata=value`: 指定fastaq元信息，必须包含["data_id","sample_id","organism"],详见utils/fastq_utils.py如何处理metadata


## changelog
1. 20251022-fix(all):将原流程拆分，统一用include导入，旧有流程丢入example，包括其configfile。(snakemake对docker的原生支持不如nextflow，还是使用本地环境，在configfile中配置)
2. 20251111-feat(run,utils):引入metadata，添加工具类统一处理metadata，在之前支持单双端文件同时执行的基础上支持多物种同时运行