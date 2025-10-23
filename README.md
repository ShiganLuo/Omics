# RNASNP

## 快速开始

1. 按照configfile的Procedure安装程序

注意：trim_galore只是打包程序，需要却cutadapt存在

## 参数详解

- --containerize启用docker
- --config input=value output=value 指定输入和输出文件路径

## changelog
1. fix(all):20251022:将原流程拆分，统一用include导入，旧有流程丢入example，包括其configfile。(snakemake对docker的原生支持不如nextflow，还是使用本地环境，在configfile中配置)