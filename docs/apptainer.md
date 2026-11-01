# Apptainer 容器模式详解

[← 返回 README](../README.md)

`run.py --sdm` 用预构建的 Apptainer SIF 镜像运行所有工具，软件不需要手动安装。本文是容器模式的完整说明：原理、镜像构建、bind 路径推导、手动调用 Snakemake、常见问题，以及脱离 Snakemake 直接用 `apptainer exec` 调工具的速查表。

## 原理

每个模块的 conda 环境 YAML（`modules/<module>/<name>.yaml`）可构建为 Apptainer SIF 镜像。Snakemake 通过 `container:` 指令在运行时自动调用对应的 SIF。

## 构建 SIF 镜像

使用 `EnvUtil.py` 从 conda YAML 生成 `.def` 文件并构建 `.sif`：

```bash
# 生成 .def 文件（不构建）
python workflow/Omics/src/common/util/EnvUtil.py gen \
  --modules-dir workflow/Omics/modules \
  --outdir /path/to/env

# 构建已有 .def 文件为 .sif
python workflow/Omics/src/common/util/EnvUtil.py build \
  --modules-dir workflow/Omics/modules \
  --outdir /path/to/env

# 完整流程：YAML → .def → .sif（推荐）
python workflow/Omics/src/common/util/EnvUtil.py all \
  --modules-dir workflow/Omics/modules \
  --outdir /path/to/env
```

构建后的 SIF 存放结构：

```
/path/to/env/
  bowtie2/
    bowtie2.sif
  macs3/
    macs3.sif
  fastqc/
    fastqc.sif
  ...
```

## SIF 路径解析

`common.smk` 中的 `sif()` 函数负责路径解析：

```python
# modules/common/common.smk
def sif(yaml_filename: str) -> str:
    """根据 config["env"] 映射查找 SIF 路径"""
```

配置 JSON 中的 `env` 字段建立 YAML stem → SIF 路径的映射：

```json
{
  "env": {
    "env_dir": "/path/to/env",
    "bowtie2": "/path/to/env/bowtie2/bowtie2.sif",
    "macs3": "/path/to/env/macs3/macs3.sif"
  }
}
```

未显式映射的模块会自动回退到 `<env_dir>/<module>/<yaml_stem>.sif`。

## Bind 路径自动推导

`run.py --sdm` 模式下，`_collect_bind_paths()` 会递归扫描配置 JSON 中的所有路径字符串，自动推导需要 bind-mount 的目录列表，并合并父目录以最小化 bind 数量。

自动 bind 包含：

- 输入数据目录（raw FASTQ 路径）
- 输出目录
- 参考基因组目录
- 软件/env 目录
- `/tmp`（临时文件）

如需额外 bind 路径：

```bash
python workflow/Omics/run.py \
  -m data/meta.tsv -w PeakCalling -o output \
  --sdm \
  --singularity-args '--bind /extra/data,/extra/ref'
```

## 手动运行 Snakemake + Apptainer

不通过 `run.py`，直接调用 Snakemake：

```bash
snakemake -s workflow/Omics/subworkflow/PeakCalling.smk \
  --configfile output/<project>/PeakCalling/raw.json \
  --cores 48 \
  --sdm apptainer \
  --singularity-args '--bind /data,/home/user/Database,/home/user/Data,/tmp' \
  --rerun-triggers code input mtime params software-env
```

## 常见问题

- **`command not found`**：SIF 内的工具不在 PATH 中，检查 `.def` 文件的 `%environment` 或 `%post` 段是否正确设置了 PATH。
- **权限错误**：确保 `--bind` 的目标目录对当前用户可读写。
- **SIF 过期**：conda YAML 更新后需重新构建 SIF（`EnvUtil.py all`）。
- **启动容器报 `Could not write info to setgroups: Permission denied` / `Error while waiting event for user namespace mappings: no event received`**：
  conda 装的非特权版 apptainer 依赖 user namespace，Ubuntu 24.04 默认 AppArmor
  禁止非特权用户命名空间（`kernel.apparmor_restrict_unprivileged_userns=1`）导致
  启动失败。解决：改用 setuid 版 `apptainer-suid`（见 README「环境依赖」）。
  修复后若 conda 环境里的非特权版 apptainer 在 PATH 中优先，会复现同样报错。
  同时出现的 `fuse2fs not found` 只是无关提示，可忽略。

## Apptainer 直接执行命令

SIF 镜像内的工具可直接通过 `apptainer exec` 调用，无需经过 Snakemake。每个 SIF 的 `%runscript` 设置为 `exec "$@"`，因此传入的命令会直接执行。

基本语法：

```bash
apptainer exec \
  --bind /data,/home/user/Database,/tmp \
  /path/to/env/<module>/<module>.sif \
  <tool> [arguments...]
```

常用工具示例：

```bash
BIND="--bind /home/user/Data,/home/user/Database,/tmp"
ENV="/home/user/Database/env"

# FastQC 质控
apptainer exec $BIND $ENV/fastqc/fastqc.sif \
  fastqc --threads 8 -o output/QC/ sample_1.fq.gz sample_2.fq.gz

# Trim Galore 接头修剪
apptainer exec $BIND $ENV/trim-galore/trim-galore.sif \
  trim_galore --quality 25 --paired -o output/trimmed/ sample_1.fq.gz sample_2.fq.gz

# Bowtie2 比对
apptainer exec $BIND $ENV/bowtie2/bowtie2.sif \
  bowtie2 -x /path/to/index -1 sample_1.fq.gz -2 sample_2.fq.gz --very-sensitive -S output.sam

# SAMtools 处理
apptainer exec $BIND $ENV/samtools/samtools.sif \
  samtools sort -@ 8 -o output.bam input.sam

# GATK MarkDuplicates
apptainer exec $BIND $ENV/gatk/gatk4.sif \
  gatk MarkDuplicates -I input.bam -O output.bam -M metrics.txt

# MACS3 Peak Calling
apptainer exec $BIND $ENV/macs3/macs3.sif \
  macs3 callpeak -t IP.bam -c Input.bam -f BAM -g mm -n sample --outdir peaks/

# HOMER 注释
apptainer exec $BIND $ENV/homer/homer.sif \
  annotatePeaks.pl peaks.narrowPeak mm39 -gtf annotation.gtf > annotation.txt

# BEDTools 交集
apptainer exec $BIND $ENV/bedtools/bedtools.sif \
  bedtools intersect -a peaks.bed -b TE.gtf -wa -wb > overlap.bed

# deeptools computeMatrix
apptainer exec $BIND $ENV/deeptools_heatmap/deeptools_heatmap.sif \
  computeMatrix reference-point -S signal.bigwig -R peaks.bed -a 3000 -b 3000 -o matrix.gz

# deeptools plotHeatmap
apptainer exec $BIND $ENV/deeptools_heatmap/deeptools_heatmap.sif \
  plotHeatmap -m matrix.gz -o heatmap.png --whatToShow 'heatmap and colorbar'

# FRiP score 计算（bedtools + samtools 组合）
apptainer exec $BIND $ENV/samtools/samtools.sif \
  samtools view -c -F 4 aligned.bam  # 总 mapped reads
apptainer exec $BIND $ENV/bedtools/bedtools.sif \
  bedtools intersect -a aligned.bam -b peaks.bed -u | samtools view -c  # peaks 内 reads
```

调试容器内部环境：

```bash
# 进入容器交互式 shell
apptainer shell --bind /home/user/Data,/home/user/Database,/tmp \
  /home/user/Database/env/bowtie2/bowtie2.sif

# 查看容器内可用命令
apptainer exec $ENV/bowtie2/bowtie2.sif which bowtie2
apptainer exec $ENV/bowtie2/bowtie2.sif bowtie2 --version
```
