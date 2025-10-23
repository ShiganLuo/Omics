# RNASNP

## 快速开始

1. 按照configfile的Procedure安装程序

注意：trim_galore只是打包程序，需要确保cutadapt存在

2. 运行

```sh
snakemake -s workflow/RNA-SNP/run.smk --config indir=data/fq outdir=output --cores 45
```

## 参数详解

- `--config`
  - `indir=value`：指定 fastq.gz 输入路径
  - `outdir=value`：指定输出文件目录
  - `genome=value`：指定分析物种，需提前在 configfile 中配置，默认 `human`


## changelog
1. 20251022-fix(all):将原流程拆分，统一用include导入，旧有流程丢入example，包括其configfile。(snakemake对docker的原生支持不如nextflow，还是使用本地环境，在configfile中配置)