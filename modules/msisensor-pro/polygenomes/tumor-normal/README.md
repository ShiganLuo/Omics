# msisensor-pro/polygenomes/tumor-normal — 多基因组 MSI 检测子模块(MSIsensor-pro,配对模式)

polygenomes 布局下的 tumor-vs-normal 微卫星不稳定(MSI)检测:对每个
基因组参考先 `scan` 建微卫星位点表,再逐配对 case 跑 `msi`,最后合并
为一张 MSI 分数表。独立 `all` 规则驱动,输入走样本表
(`sample_path` TSV)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `sample_path` (config) | 样本表 TSV(首列作 index=case 名),须含 `tumor_path` / `normal_path` 两列(BAM 路径),由 pandas 读入 |
| `reference` (config) | 当前基因组参考 FASTA(scan/msi 的 `-g`/`-d` 输入) |
| `genome_version` (config) | 基因组 ID(进输出路径与文件名;wildcard `{genome_version}`) |
| `outdir` (config) | 输出目录(末尾自动补 `/`) |
| `genome.<g>.fasta` (config) | 示例 config 里列出的多基因组参考字典;本 .smk 实际只消费 `reference` 键 |

环境文件:规则声明 `conda: "../../msisensor-pro.yaml"`(共享祖父目录的
msisensor-pro 环境)。

### 1.2 输出

```
{outdir}/reference/{genome_version}.msisensor.scan.list              # 微卫星位点表(scan)
{outdir}/tumor_normal_output/{case}/{case}.{genome_version}.msisensor # 每 case MSI 结果
{outdir}/tumor_normal_output.{genome_version}.merge.tsv               # 全 case 合并表
```

日志:`{outdir}/reference/{genome_version}.msisensor.scan.log` 等(日志
与输出同目录,无独立 logdir)。`merge_msi_result` 合并表列:
`case_name`、`Total_number_of_sites`、`Number_of_unstable_sites`、
`MSI_score`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `sample_path` | — | — | 样本表 TSV,必需(顶层 config 键) |
| `reference` | `scan -d` / `msi -g` | — | 参考 FASTA,必需 |
| `genome_version` | — | — | 基因组 ID,必需(进文件名) |
| `msisensor-pro` | 引擎可执行 | `msisensor-pro` | 顶层 config 键(非 `Procedure` 嵌套) |
| (mem_mb) | — | 8000(scan/merge)/ 16000(msi) | `rule_res` 可覆盖 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
msisensor-pro scan -d genome1.fa \
    -o out/reference/genome1.msisensor.scan.list

msisensor-pro msi \
    -d out/reference/genome1.msisensor.scan.list \
    -n normal.bam -t tumor.bam -g genome1.fa \
    -o out/tumor_normal_output/case1/case1.genome1.msisensor
```

## 2. 算法原理与生物学意义

MSIsensor-pro 对参考基因组 `scan` 找出微卫星(短串联重复)位点并记录
各等位长度;`msi` 阶段比较 tumor 与 matched normal 在这些位点上的
读段长度分布,长度分布显著漂移(插入/缺失滑移)的位点判为不稳定,
输出三列:总位点数、不稳定位点数、MSI_score(不稳定占比)。

- **MSI_score** ≥ 2 常用作 MSI-H(微卫星高度不稳定)阈值(MSIsensor
  系惯例;本仓库 qc_summary/metrics.md 同此约定),MSI-H 是 dMMR/
  错配修复缺陷的标志,提示免疫治疗(如 PD-1 抗体)获益可能;
- 配对 normal 校正个体胚系的重复长度多态,tumor-only 阈值可靠性下降
  (对应子模块 `tumor-only/`);
- polygenomes 变体对每个 `genome_version` 独立 scan+msi:同一对样本在
  不同参考下的位点表不同,MSI 分数需按基因组分别解读,合并表文件名
  带 `.{genome_version}` 区分。

## 3. 参考文献

1. MSIsensor-pro — https://github.com/bioinform/msisensor-pro
   (模块文件内未附论文引文;配套指标解读见本仓库
   `modules/qc_summary/metrics.md` §7)
