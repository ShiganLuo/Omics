# RNASNP

## 快速开始

1. 创建envs中的conda环境，对于各个subworflow，配置好同名yaml内容

注意：
- trim_galore只是打包程序，需要确保cutadapt存在
- 物种名是动态变化的，在yaml中，相关文件的键需要根据metadata改变，如果物种名为有空格隔开，需要改空格为_，如Mus musculus -> Mus_musculus

2. 修改scripts下RNA-SNP_prepare.sh脚本，执行

3. 运行

```sh
snakemake -s workflow/RNA-SNP/run.smk --config indir=data/fq outdir=output metadata=data/target_fq.tsv--cores 45 --use-conda
```
4. 亮点

- 支持多个物种
- 支持单端测序和双端测序文件
- 统一的日志管理

注意：

双端测序文件：
  _R1*fastq , _R2*fastq or _R1*fq , _R2*fq
  _1*fastq , _2*fastq or _1*fq , _2*fq
单端测序文件：非双端测序文件

## 参数详解

- `--config`
  - `indir=value`：指定 fastq.gz 输入路径
  - `outdir=value`：指定输出文件目录
  - `metadata=value`: 指定fastaq元信息，必须包含["data_id","sample_id","organism"],详见utils/fastq_utils.py如何处理metadata

