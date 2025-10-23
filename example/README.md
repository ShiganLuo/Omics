# Description

fastq文件必须在fq/下，而且后缀名为.fq.gz，前缀为SRR。

单端为：SRR[number*n].fq.gz

双端为: SRR[number*n]_1.fq.gz,SRR[number*n]_2.fq.gz

说明：

- 本流程在鉴定人类样本SNP的去除了小鼠基因组污染，在进行转录组分析时没必要去除小鼠基因组污染
- 在鉴定小鼠样本SNP的时候，由于是同物种，没必要有额外步骤；转录组分析更是如此

## workflow

### **运行**

#### 如果在snakefile所在目录运行

`time snakemake -s snakefile --cores {cores>25}} --use-conda`

#### 如果在非snakefile目录运行

`time snakemake -s /your/snakefile/path --cores 25 --directory /the/directroy/of/snakefile --use-conda`

1. --config indir="data/other" outdir="results/other"
    - 可自定义输入输出目录，默认输入目录为README.md同目录下的data,输出目录为同目录下的output
注意：目录不一定是相对于工作目录而言，指定了--directory，是相对该目录而言
举例：

```shell
snakemake -s workflow/snakemake/RNASNP-noRemove.smk --config indir="../../data/GSE166216" outdir="../../output/GSE166216" --cores 25 --directory workflow/snakemake/ --use-conda --dry-run
```

注意事项：
    - fq.gz 不允许fastaq.gz由于rule triming
    - 最好不要改变workflow层级关系

### RNASNP-human.smk

- 适合人类（小鼠饲养层细胞培养，污染了小鼠基因组序列）

### RNASNP-noRemove.smk

适合小鼠等没有其它基因组污染的样本

### prepare and downstream

#### prepare

运行前需要做的准备工作

#### downstream

在运行之前，你需要根据实际情况更改assets中的group.csv分组信息或修改downstrea.sh以适应新的分组信息

## todo

1.合并单双端测序文件处理流程

    - 20250329: 暂时分开单双端测序流程，再原双端测序中加入single_samples变量，使得不至于影响原来双端测序流程，
    又能兼容先运行单端测序流程，再运行双端测序流程
    - 20250329: 合并单双端测序流程，并完善人类流程
    - 20251021: 拆分模块，将成型流程置于example，根目录为当前运行流程run.smk
