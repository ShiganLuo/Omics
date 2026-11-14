# manta/polygenomes — 多参考基因组 Manta 结构变异检测子模块

Manta(结构变异 calling)的多参考基因组(polygenomes)变体:按 `genome`
wildcard 为每样本 × 每基因组配置并运行 Manta,输出候选 SV 的 VCF。conda
环境共享父目录 `../manta.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 样本 ID 布局(形态 A) | `{indir}/{sample_id}/{sample_id}.bam`:坐标排序 BAM(比对结果) |
| `genome.<G>.fasta` | 每基因组参考 FASTA;`<G>` 为 genome wildcard |
| `samples` | 样本 ID 列表(供 subworkflow 注册驱动 DAG) |
| `ROOT_DIR` | 项目根目录(config 声明,当前规则未直接消费) |

### 1.2 输出

Manta 的 `--runDir` 固定为 `{outdir}/{genome}/{sample_id}/`,产物按 Manta
目录约定落位:

```
{outdir}/{genome}/{sample_id}/runWorkflow.py.config.pickle   # manta_config 产出的工作流配置
{outdir}/{genome}/{sample_id}/runWorkflow.py                 # manta_config 生成的运行脚本(manta_run 的输入)
{outdir}/{genome}/{sample_id}/results/variants/candidateSV.vcf.gz      # 候选结构变异
{outdir}/{genome}/{sample_id}/results/variants/candidateSV.vcf.gz.tbi  # tabix 索引
{outdir}/{genome}/{sample_id}/manta_config.sh / manta_run.sh # 回执
{logdir}/{genome}/{sample_id}/manta_config.log / manta_run.log
```

收尾规则 `manta_result` 聚合 candidateSV VCF + TBI。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.manta.threads` | `runWorkflow.py -j` | 8 | Manta 工作流并行 jobs(亦为规则 threads) |
| `genome.<G>.fasta` | `configManta.py --referenceFasta` | — | 每基因组参考 FASTA |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

模块无 `Procedure` 工具路径配置:`configManta.py` 与 `runWorkflow.py` 直接
调用(来自 conda 环境/容器内的 Manta 安装)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
configManta.py --bam s1.bam --referenceFasta GRCh38.fa \
  --runDir output/GRCh38/s1
python output/GRCh38/s1/runWorkflow.py -j 8
```

## 2. 算法原理与生物学意义

Manta 通过 read-pair 与 split-read 证据发现候选断点,对每个候选区间做局部
de novo 组装,识别缺失(DEL)、重复(DUP)、倒位(INV)、插入(INS)与
断点重排(BND),并对胚系样本给出基因型。本模块输出 `candidateSV.vcf.gz`
(候选集,未做进一步打分/过滤)。生物学意义:结构变异影响基因剂量与调控
结构,是遗传病、肿瘤与演化研究的关键变异类型;长/短读长比对的成对证据均可
使用。polygenomes 场景下同一样本对多参考基因组分别检测(如人/鼠跨物种
样本)。结果阅读:VCF 记录类型与断点 CIPOS/CIEND 置信区间;过滤 PASS
候选;`manta_result` 产物可直接送注释/可视化。

## 3. 参考文献

1. Manta —— https://github.com/Illumina/manta(源码与文档)
2. Chen X, et al. *Manta: rapid detection of structural variants and indels
   for germline and cancer sequencing applications*. Bioinformatics 2016,
   32(8):1220–1222. DOI: 10.1093/bioinformatics/btv710
